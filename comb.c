
/**
 * @file comb.c
 *
 * @brief comb aligner impl.
 *
 * @author Hajime Suzuki
 * @date 2016/4/10
 * @license Apache v2.
 */
#define UNITTEST_UNIQUE_ID		5
#include "unittest.h"

#include <ctype.h>
#include <stdint.h>
#include "ptask.h"			/* parallel task dispatcher */
#include "sr.h"				/* sequence reader */
#include "gref.h"			/* graphical sequence indexer */
#include "ggsea.h"			/* graph-to-graph seed-and-extend alignment */
#include "aw.h"				/* alignment writer */
#include "mem.h"
#include "sassert.h"
#include "log.h"
#include "lmm.h"


/* inline directive */
#define _force_inline				inline

/* max and min */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )


/* constants */
#ifndef COMB_VERSION_STRING
#define COMB_VERSION_STRING			"0.0.1"
#endif


/* alignment core functions */
/**
 * @struct comb_align_params_s
 */
struct comb_align_params_s {
	/* global */
	uint8_t help;
	uint8_t version;
	uint8_t pad1[5];

	uint8_t message_level;
	int (*message_printer)(void *ctx, char const *fmt, ...);
	void *message_context;

	int64_t num_threads;
	int64_t mem_size;
	int64_t pool_size;

	char *command;
	char *command_base;
	char *program_name;
	int64_t program_id;

	char *ref_name;
	char *query_name;
	char *out_name;

	uint8_t ref_format;
	uint8_t query_format;
	uint8_t out_format;

	/* indexing parameters */
	int64_t k;

	/* filtering parameters */
	int64_t kmer_cnt_thresh;
	int64_t overlap_thresh;
	int64_t gapless_thresh;

	/* scoring parameters */
	int64_t xdrop;
	int8_t m, x, gi, ge;
	char clip;
	int8_t pad2[2];

	/* reporting parameters */
	uint8_t include_unmapped;
	int64_t score_thresh;
};

/**
 * @struct comb_align_worker_args_s
 */
struct comb_align_worker_args_s {
	struct comb_align_params_s const *params;
	struct sr_gref_s *r;
	ggsea_ctx_t *ctx;
	sr_t *ref;
	sr_t *query;
	aw_t *aw;
};

/**
 * @struct comb_align_worker_item_s
 */
struct comb_align_worker_item_s {
	lmm_t *lmm;
	struct sr_gref_s *q;
	struct ggsea_result_s *res;
};

/* multithread workers */
/**
 * @fn comb_align_worker_init
 */
static _force_inline
struct comb_align_worker_args_s **comb_align_worker_init(
	struct comb_align_params_s const *params,
	ggsea_conf_t const *conf,
	sr_t *ref,
	sr_t *query,
	aw_t *aw)
{
	int64_t num_worker = MAX2(1, params->num_threads);
	struct comb_align_worker_args_s **w = (struct comb_align_worker_args_s **)malloc(
		  (num_worker + 1) * sizeof(struct comb_align_worker_args_s *)
		+ (num_worker + 1) * sizeof(struct comb_align_worker_args_s));
	struct comb_align_worker_args_s *base = (struct comb_align_worker_args_s *)(w + num_worker + 1);

	/* set pointers */
	for(int64_t i = 0; i < num_worker + 1; i++) {
		w[i] = &base[i];
	}

	/* init worker object */
	for(int64_t i = 0; i < num_worker; i++) {
		struct sr_gref_s *r = sr_get_index(ref);
		base[i] = (struct comb_align_worker_args_s){
			.params = params,
			.r = r,
			.ctx = ggsea_ctx_init(conf, r->gref),
			.ref = ref,
			.query = query,
			.aw = aw
		};
	}
	memset(&base[num_worker], 0, sizeof(struct comb_align_worker_args_s));
	return(w);
}

/**
 * @fn comb_align_worker_clean
 */
static _force_inline
void comb_align_worker_clean(
	struct comb_align_worker_args_s **w)
{
	if(w == NULL) { return; }

	for(struct comb_align_worker_args_s **p = w; (*p)->params != NULL; p++) {
		ggsea_ctx_clean((*p)->ctx);
		sr_gref_free((*p)->r);
		memset(*p, 0, sizeof(struct comb_align_worker_args_s));
	}
	free(w);
	return;
}

/**
 * @fn comb_align_source
 */
static
void *comb_align_source(
	void *arg)
{
	struct comb_align_worker_args_s *a = (struct comb_align_worker_args_s *)arg;

	/* get the next iterator */
	struct sr_gref_s *q = NULL;
	while((q = sr_get_iter(a->query)) != NULL) {
		debug("iterator fetched a(%p), q(%p), iter(%p)", a, q, q->iter);

		if(q->iter != NULL) {
			/* valid iterator was returned, create worker_item object */
			struct comb_align_worker_item_s *i = (struct comb_align_worker_item_s *)lmm_malloc(
				q->lmm, sizeof(struct comb_align_worker_item_s));

			/* init item */
			*i = (struct comb_align_worker_item_s){
				.lmm = q->lmm,
				.q = q
			};

			debug("worker created, ptr(%p)", i);
			return((void *)i);
		}
	}

	/* reached the end */
	return(NULL);
}

/**
 * @fn comb_align_worker
 */
static
void *comb_align_worker(
	void *arg,
	void *item)
{
	struct comb_align_worker_args_s *a = (struct comb_align_worker_args_s *)arg;
	struct comb_align_worker_item_s *i = (struct comb_align_worker_item_s *)item;

	/* do alignment */
	debug("align a(%p), i(%p), iter(%p)", a, i, i->q->iter);
	i->res = ggsea_align(a->ctx, i->q->gref, i->q->iter, i->lmm);
	return((void *)i);
}

/**
 * @fn comb_align_drain
 */
static
void comb_align_drain(
	void *arg,
	void *item)
{
	struct comb_align_worker_args_s *a = (struct comb_align_worker_args_s *)arg;
	struct comb_align_worker_item_s *i = (struct comb_align_worker_item_s *)item;

	/* append to result queue */
	if(i->res->cnt == 0 && a->params->include_unmapped != 0) {
		aw_append_unmapped(a->aw, i->res->ref, i->res->query);
	} else {
		aw_append_alignment(a->aw, i->res->ref, i->res->query, i->res->aln, i->res->cnt);
	}

	/* cleanup */
	ggsea_aln_free(i->res);

	/* lmm must be freed before gref_free */
	debug("worker destroyed, ptr(%p)", i);
	struct sr_gref_s *q = i->q;
	lmm_free(i->lmm, (void *)i);
	sr_gref_free(q);
	return;
}


/* align main functions */
/**
 * @fn comb_align_print_option_summary
 */
static _force_inline
void comb_align_print_option_summary(
	struct comb_align_params_s const *params)
{
	char *s = (char *)malloc(4096);
	char *p = s;

	#define _p(...)		p += sprintf(p, __VA_ARGS__);
	_p("%s align ", params->command_base);
	_p("-t%" PRId64 " ", params->num_threads);
	_p("-k%" PRId64 " ", params->k);
	_p("-r%" PRId64 " ", params->kmer_cnt_thresh);
	_p("-d%" PRId64 " ", params->overlap_thresh);
	_p("-f%" PRId64 " ", params->gapless_thresh);
	_p("-a%d ", params->m);
	_p("-b%d ", params->x);
	_p("-p%d ", params->gi);
	_p("-q%d ", params->ge);
	_p("-x%" PRId64 " ", params->xdrop);
	_p("-m%" PRId64 " ", params->score_thresh);
	_p("-c%c", params->clip);
	#undef _p

	params->message_printer(params->message_context, "%s\n", s);
	free(s);
	return;
}

/**
 * @fn comb_align
 */
static
int comb_align(
	struct comb_align_params_s const *params)
{
	int ret = 1;

	ggsea_conf_t *conf = NULL;
	sr_t *ref = NULL;
	sr_t *query = NULL;
	aw_t *aw = NULL;
	struct comb_align_worker_args_s **w = NULL;
	ptask_t *pt = NULL;

	#define comb_align_error(expr, ...) { \
		if(!(expr)) { \
			if(params->message_level != 0) { \
				params->message_printer(params->message_context, "[ERROR] " __VA_ARGS__); \
			} \
			goto _comb_align_error_handler; \
		} \
	}

	/* print option summary */
	if(params->message_level != 0) {
		comb_align_print_option_summary(params);
	}

	/* build ggsea configuration object */
	debug("build conf");
	conf = ggsea_conf_init(GGSEA_PARAMS(
		.xdrop = params->xdrop,
		.score_matrix = GABA_SCORE_SIMPLE(params->m, params->x, params->gi, params->ge),
		.k = params->k,
		.kmer_cnt_thresh = params->kmer_cnt_thresh,
		.overlap_thresh = params->overlap_thresh,
		.gapless_thresh = params->gapless_thresh,
		.score_thresh = params->score_thresh));
	comb_align_error(conf != NULL, "Failed to create alignment configuration.");

	/* build reference sequence index */
	ref = sr_init(params->ref_name,
		SR_PARAMS(
			.format = params->ref_format,
			.k = params->k,
			.seq_direction = SR_FW_ONLY,
			.num_threads = params->num_threads
		));
	comb_align_error(ref != NULL, "Failed to open reference file `%s'.", params->ref_name);

	/* build read pool */
	query = sr_init(params->query_name,
		SR_PARAMS(
			.format = params->query_format,
			.k = params->k,
			.seq_direction = SR_FW_ONLY,
			.pool_size = params->pool_size,
		));
	comb_align_error(query != NULL, "Failed to open query file `%s'.", params->query_name);

	/* build alignment writer */
	struct sr_gref_s *r = sr_get_index(ref);
	comb_align_error(r != NULL, "Failed to build reference index.");
	aw = aw_init(params->out_name, r->gref,
		AW_PARAMS(
			.format = params->out_format,
			.clip = params->clip,
			.program_id = params->program_id,
			.program_name = params->program_name,
			.command = params->command
		));
	sr_gref_free(r);
	comb_align_error(aw != NULL, "Failed to open output file `%s'.", params->out_name);

	/* initialize parallel task dispatcher */
	w = comb_align_worker_init(params, conf, ref, query, aw);
	pt = ptask_init(comb_align_worker, (void **)w, params->num_threads, params->pool_size);
	comb_align_error(w != NULL && pt != NULL, "Failed to initialize parallel worker threads.");

	/* run tasks */
	ptask_stream(pt, comb_align_source, (void *)*w, comb_align_drain, (void *)*w, params->pool_size / 4);

	/* destroy objects */
	ret = 0;
_comb_align_error_handler:;
	if(ref != query) {
		sr_clean(query); query = NULL;
	}
	sr_clean(ref); ref = NULL;
	aw_clean(aw); aw = NULL;
	comb_align_worker_clean(w); w = NULL;
	ggsea_conf_clean(conf); conf = NULL;
	ptask_clean(pt); pt = NULL;
	return(ret);
}

/**
 * @fn comb_align_clean
 */
static _force_inline
void comb_align_clean(
	struct comb_align_params_s *params)
{
	if(params == NULL) {
		return;
	}

	free(params->command); params->command = NULL;
	free(params->command_base); params->command_base = NULL;
	free(params->program_name); params->program_name = NULL;
	free(params->ref_name); params->ref_name = NULL;
	free(params->query_name); params->query_name = NULL;
	free(params->out_name); params->out_name = NULL;
	free(params);
	return;
}



/* index core functions (not yet implemented) */
/**
 * @struct comb_index_params_s
 */
struct comb_index_params_s {
	/* global */
	uint8_t help;
	uint8_t version;
	uint8_t pad1[5];

	uint8_t message_level;
	int (*message_printer)(void *ctx, char const *fmt, ...);
	void *message_context;

	int64_t num_threads;
	int64_t mem_size;

	char *command;
	char *command_base;
	char *program_name;
	int64_t program_id;

	char *ref_name;
	char *prefix;

	/* indexing parameters */
	int64_t k;
};

/**
 * @fn comb_index_print_option_summary
 */
static _force_inline
void comb_index_print_option_summary(
	struct comb_index_params_s const *params)
{
	char *s = (char *)malloc(4096);
	char *p = s;

	#define _p(...)		p += sprintf(p, __VA_ARGS__);
	_p("%s index ", params->command_base);
	_p("-t%" PRId64 " ", params->num_threads);
	_p("-k%" PRId64 " ", params->k);
	#undef _p

	params->message_printer(params->message_context, "%s\n", s);
	free(s);
	return;
}

/**
 * @fn comb_index
 */
static
int comb_index(
	struct comb_index_params_s const *params)
{
	if(params->message_level != 0) {
		comb_index_print_option_summary(params);
	}
	return(1);
}

/**
 * @fn comb_index_clean
 */
static _force_inline
void comb_index_clean(
	struct comb_index_params_s *params)
{
	if(params == NULL) {
		return;
	}

	free(params->command); params->command = NULL;
	free(params->command_base); params->command_base = NULL;
	free(params->program_name); params->program_name = NULL;
	free(params->ref_name); params->ref_name = NULL;
	free(params);
	return;
}


/* subcommand (compatibility layer) implementations */
/**
 * @struct comb_inst_s
 */
struct comb_inst_s {
	void *params;
	int (*process)(void *);
	void (*clean)(void *);
};

/**
 * @struct comb_args_s
 */
struct comb_args_s {
	int argc;
	char **argv;
};


/* message printers */
/**
 * @fn comb_print_version
 */
static _force_inline
void comb_print_version(
	void)
{
	fprintf(stderr, "comb aligner (%s)\n", COMB_VERSION_STRING);
	return;
}

/**
 * @fn comb_print_unknown_option
 */
static _force_inline
void comb_print_unknown_option(
	char c)
{
	fprintf(stderr, "[WARNING] Unknown option `%c'.\n", c);
	return;
}

/**
 * @fn comb_print_invalid_args
 */
static _force_inline
void comb_print_invalid_args(
	char *argv[],
	int from,
	int to)
{
	fprintf(stderr, "[ERROR] Invalid number of arguments.\n");
	return;
}

/**
 * @fn comb_print_error
 */
static _force_inline
void comb_print_error(
	char const *msg,
	...)
{
	va_list l;
	va_start(l, msg);
	fprintf(stderr, "[ERROR] ");
	vfprintf(stderr, msg, l);
	fprintf(stderr, "\n");
	va_end(l);
	return;
}


/* option parset helper functions */
/**
 * @fn comb_build_short_option_string
 */
static _force_inline
char *comb_build_short_option_string(
	struct option const *opts)
{
	char *str = NULL, *ps;
	int len = 0;
	struct option const *po = NULL;

	for(po = opts; po->val != 0; po++) { len++; }
	str = ps = (char *)malloc(2 * len + 1);
	for(po = opts; po->val != 0; po++) {
		if(!isalnum(po->val)) { continue; }
		*ps++ = (char)po->val;
		if(po->has_arg != no_argument) {
			*ps++ = ':';
		}
	}
	*ps = '\0';
	return(str);
}

/**
 * @fn comb_build_command_string
 */
static _force_inline
char *comb_build_command_string(
	int argc,
	char *argv[])
{
	lmm_kvec_t(char) s;

	/* pass NULL to force use libc malloc (to be freed with libc free) */
	lmm_kv_init(NULL, s);

	for(int i = 0; i < argc; i++) {
		char const *p = argv[i];
		while(*p != '\0') {
			lmm_kv_push(NULL, s, *p);
			p++;
		}
		lmm_kv_push(NULL, s, ' ');
	}
	lmm_kv_at(s, lmm_kv_size(s) - 1) = '\0';
	return(lmm_kv_ptr(s));
}

/**
 * @fn comb_atoi
 */
static _force_inline
int64_t comb_atoi_prefix(
	char const *val)
{
	if(*val == 'x') {
		return(strtoll(++val, NULL, 16));
	} else if(*val == 'd') {
		return(strtoll(++val, NULL, 10));
	} else if(*val == 'b') {
		return(strtoll(++val, NULL, 2));
	}
	return(strtoll(val, NULL, 8));
}
static _force_inline
int64_t comb_atoi_dec(
	char const *val)
{
	int64_t len = strlen(val);
	int64_t mul = 1, base = 1000;
	if(isdigit(val[len - 1])) {
		return(strtoll(val, NULL, 10));
	}
	do {
		switch(val[len - 1]) {
			case 'i': base = 1024; len--; continue;
			case 'T': mul *= base;
			case 'G': mul *= base;
			case 'M': mul *= base;
			case 'K':
			case 'k': mul *= base;
			default:;
				char tmp[len];
				memcpy(tmp, val, len - 1);
				tmp[len - 1] = '\0';
				return(strtoll(tmp, NULL, 10));
		}
	} while(0);
	return(0);
}
static _force_inline
int64_t comb_atoi(
	char const *val)
{
	if(*val == '0') {
		return(comb_atoi_prefix(val + 1));
	} else {
		return(comb_atoi_dec(val));
	}
}

/**
 * @fn comb_parse_format
 */
static _force_inline
int comb_parse_format(
	char const *str)
{
	struct format_map_s {
		char const *str;
		int format;
	};

	static struct format_map_s const map[] = {
		{ "fasta", SR_FASTA },
		{ "fa", SR_FASTA },
		{ "fastq", SR_FASTQ },
		{ "fq", SR_FASTQ },
		{ "fast5", SR_FAST5 },
		{ "f5", SR_FAST5 },
		{ "sam", AW_SAM },
		{ "bam", AW_BAM },
		{ "maf", AW_MAF },
		{ "gpa", AW_GPA },
		{ NULL, 0 }
	};
	for(uint64_t i = 0; i < sizeof(map) / sizeof(struct format_map_s); i++) {
		if(map[i].str == NULL) { continue; }
		if(strcmp(str, map[i].str) == 0) {
			return(map[i].format);
		}

	}
	return(0);
}


/* `unittest' subcommand implementation */
/**
 * @fn comb_process_unittest
 */
static
int comb_process_unittest(void *params)
{
	struct comb_args_s *args = (struct comb_args_s *)params;
	return(unittest_main(args->argc, args->argv));
}

/**
 * @fn comb_init_unittest
 */
static
struct comb_inst_s *comb_init_unittest(
	int argc,
	char *argv[])
{
	struct comb_inst_s *o = (struct comb_inst_s *)malloc(
		sizeof(struct comb_inst_s) + sizeof(struct comb_args_s));
	if(o == NULL) { return(NULL); }

	struct comb_args_s *args = (struct comb_args_s *)(o + 1);
	*args = (struct comb_args_s){
		.argc = argc,
		.argv = argv
	};

	*o = (struct comb_inst_s){
		.params = (void *)args,
		.process = comb_process_unittest,
		.clean = NULL
	};
	return(o);
}


/* `align' subcommand option parser implementation */
static
char const *const comb_align_help_message =
	"\n"
	"    comb aligner (%s)\n"
	"\n"
	"  Comb aligner is a prototype implementation of a seed-and-extend alignment\n"
	"on two string graphs. The aligner accept FASTA / FASTQ and GFA formats for the\n"
	"input files (reference and query) and handle SAM and GPA (Graphical Pairwise\n"
	"Alignment format) for the output file.\n"
	"\n"
	"  Usage\n"
	"\n"
	"    $ comb align [options] <reference> <query> <output>\n"
	"\n"
	"  Options and defaults\n"
	"    Global option\n"
	"      -t<int>  [0]  Number of threads.\n"
	"\n"
	"    Seeding option\n"
	"      -k<int>  [14] k-mer length in indexing and matching.\n"
	"\n"
	"    Filtering options\n"
	"      -r<int>  [30] Repetitive k-mer filter threshold.\n"
	"      -d<int>  [3]  Overlap filter threshold.\n"
	"      -f<int>  [10] Gapless alignment filter threshold.\n"
	"\n"
	"    Extension options\n"
	"      -a<int>  [1]  Match award (in positive integer)\n"
	"      -b<int>  [1]  Mismatch penalty (in positive integer)\n"
	"      -p<int>  [1]  Gap-open penalty (pos. int. or 0 (=linear-gap penalty))\n"
	"      -q<int>  [1]  Gap-extension penalty (positive integer)\n"
	"      -x<int>  [60] X-drop threshold\n"
	"\n"
	"    Reporting options\n"
	"      -m<int>  [10] Minimum score for reporting.\n"
	"      -c<char> [S]  Clip operation in CIGAR string. (H (hard) or S (soft))\n"
	"\n"
	"    Miscellaneous options\n"
	"      -h       Print help (this) message.\n"
	"      -v       Print version information.\n"
	"\n";

/**
 * @fn comb_align_print_help
 */
static _force_inline
void comb_align_print_help(
	void)
{
	fprintf(stderr, comb_align_help_message, COMB_VERSION_STRING);
	return;
}

/* default scoring parameters */
/**
 * @fn comb_init_align_default_gapless_thresh
 */
static _force_inline
int64_t comb_init_align_default_gapless_thresh(
	struct comb_align_params_s const *params)
{
	return((int64_t)(15.0 * params->x / (params->m + params->x)));
}

/**
 * @fn comb_init_align_default_xdrop_thresh
 */
static _force_inline
int64_t comb_init_align_default_xdrop_thresh(
	struct comb_align_params_s const *params)
{
	return(15 * (params->m + 2 * params->ge) + params->gi);
}

/**
 * @fn comb_init_align_default_score_thresh
 */
static _force_inline
int64_t comb_init_align_default_score_thresh(
	struct comb_align_params_s const *params)
{
	return(100 * params->m);
}

/**
 * @fn comb_init_align
 */
#define ID_BASE					( 256 )
#define ID_REF_FORMAT			( ID_BASE + 1 )
#define ID_QUERY_FORMAT			( ID_BASE + 2 )
#define ID_OUT_FORMAT			( ID_BASE + 3 )
#define ID_INCLUDE_UNMAPPED		( ID_BASE + 4 )
#define ID_OMIT_UNMAPPED		( ID_BASE + 5 )
static
struct comb_align_params_s *comb_init_align(
	char const *base,
	int argc,
	char *argv[])
{
	if(argc == 1) {
		comb_align_print_help();
		return(NULL);
	}

	/* build params object with default params */
	struct comb_align_params_s *params = (struct comb_align_params_s *)malloc(
		sizeof(struct comb_align_params_s));
	*params = (struct comb_align_params_s){
		.message_level = 1,
		.message_printer = (int (*)(void *, char const *, ...))fprintf,
		.message_context = (void *)stderr,
		.num_threads = 0,
		.mem_size = mem_estimate_free_size(),
		.pool_size = 256,
		.command = comb_build_command_string(argc, argv),
		.command_base = strdup(base),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.k = 14,
		.kmer_cnt_thresh = 30,
		.overlap_thresh = 3,
		.gapless_thresh = 0,
		.xdrop = 0,		/* default xdrop threshold is derived from scoring parameters */
		.m = 1, .x = 2, .gi = 2, .ge = 1,
		.clip = 'S',
		.include_unmapped = 1,
		.score_thresh = 0
	};

	static struct option const opts_long[] = {
		/* global */
		{ "help", no_argument, NULL, 'h' },
		{ "version", no_argument, NULL, 'v' },
		{ "verbose", no_argument, NULL, 'V' },
		{ "threads", required_argument, NULL, 't' },
		{ "memory", required_argument, NULL, 'M' },
		{ "out", required_argument, NULL, 'o' },

		/* file formats */
		{ "ref-format", required_argument, NULL, ID_REF_FORMAT },
		{ "query-format", required_argument, NULL, ID_QUERY_FORMAT },
		{ "output-format", required_argument, NULL, ID_OUT_FORMAT },
		{ "rf", required_argument, NULL, ID_REF_FORMAT },
		{ "qf", required_argument, NULL, ID_QUERY_FORMAT },
		{ "of", required_argument, NULL, ID_OUT_FORMAT },

		/* indexing params */
		{ "seed-length", required_argument, NULL, 'k' },

		/* filtering params */
		{ "repcnt", required_argument, NULL, 'r' },
		{ "depth", required_argument, NULL, 'd' },
		{ "popcnt", required_argument, NULL, 'f' },

		/* scoring params */
		{ "match", required_argument, NULL, 'a' },
		{ "mismatch", required_argument, NULL, 'b' },
		{ "gap-open", required_argument, NULL, 'p' },
		{ "gap-extend", required_argument, NULL, 'q' },
		{ "xdrop", required_argument, NULL, 'x' },
		{ "clip-penalty", required_argument, NULL, 'C' },

		/* reporting params */
		{ "min", required_argument, NULL, 'm' },
		{ "clip", required_argument, NULL, 'c' },
		{ "include-unmapped", no_argument, NULL, ID_INCLUDE_UNMAPPED },
		{ "omit-unmapped", no_argument, NULL, ID_OMIT_UNMAPPED },
		{ 0 }
	};
	char *opts_short = comb_build_short_option_string(opts_long);

	/* optional arguments */
	int c, idx;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			/* global */
			case 'h': comb_align_print_help(); goto _comb_init_align_error_handler;
			case 'v': comb_print_version(); goto _comb_init_align_error_handler;
			case 'V': params->message_level = 3; break;
			case 't': params->num_threads = comb_atoi(optarg); break;
			case 'M': params->mem_size = comb_atoi(optarg); break;
			case 'o': params->out_name = strdup(optarg); break;

			/* formats */
			case ID_REF_FORMAT: params->ref_format = comb_parse_format(optarg); break;
			case ID_QUERY_FORMAT: params->query_format = comb_parse_format(optarg); break;
			case ID_OUT_FORMAT: params->out_format = comb_parse_format(optarg); break;

			/* params */
			case 'k': params->k = comb_atoi(optarg); break;
			case 'r': params->kmer_cnt_thresh = comb_atoi(optarg); break;
			case 'd': params->overlap_thresh = comb_atoi(optarg); break;
			case 'f': params->gapless_thresh = comb_atoi(optarg); break;
			case 'a': params->m = comb_atoi(optarg); break;
			case 'b': params->x = comb_atoi(optarg); break;
			case 'p': params->gi = comb_atoi(optarg); break;
			case 'q': params->ge = comb_atoi(optarg); break;
			case 'x': params->xdrop = comb_atoi(optarg); break;
			case 'm': params->score_thresh = comb_atoi(optarg); break;
			case 'c': params->clip = optarg[0]; break;
			case ID_INCLUDE_UNMAPPED: params->include_unmapped = 1; break;
			case ID_OMIT_UNMAPPED: params->include_unmapped = 0; break;
			
			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
	}

	/* restore default params */
	if(params->gapless_thresh == 0) {
		params->gapless_thresh = comb_init_align_default_gapless_thresh(params);
	}
	if(params->xdrop == 0) {
		params->xdrop = comb_init_align_default_xdrop_thresh(params);
	}
	if(params->score_thresh == 0) {
		params->score_thresh = comb_init_align_default_score_thresh(params);
	}

	/* positional arguments */
	int argcnt = argc - optind;
	switch(argcnt) {
		case 1:
			/* all-to-all alignment mode */
			params->ref_name = strdup(argv[optind]);
			params->query_name = strdup(argv[optind++]);
			break;
		case 2:
			/* mapping mode */
			params->ref_name = strdup(argv[optind++]);
			params->query_name = strdup(argv[optind++]);
			break;
		case 3:
			/* paired-end mapping mode */
			comb_print_error("paired-end mapping mode is not implemented.");
			goto _comb_init_align_error_handler;
		default:
			comb_print_invalid_args(argv, optind, argc);
			goto _comb_init_align_error_handler;
	}

	if(params->out_name == NULL) {
		params->out_name = strdup("-");
	}

	free(opts_short);
	return(params);

_comb_init_align_error_handler:;
	free(opts_short);
	free(params);
	return(NULL);
}


/* `index' subcommand option parser implementation */
static
char const *const comb_index_help_message =
	"\n"
	"    comb aligner (%s) index subcommand\n"
	"\n";

/**
 * @fn comb_index_print_help
 */
static _force_inline
void comb_index_print_help(
	void)
{
	fprintf(stderr, comb_index_help_message, COMB_VERSION_STRING);
	return;
}

/**
 * @fn comb_init_index
 */
static
struct comb_index_params_s *comb_init_index(
	char const *base,
	int argc,
	char *argv[])
{
	if(argc == 1) {
		comb_index_print_help();
		return(NULL);
	}

	/* build params object with default params */
	struct comb_index_params_s *params = (struct comb_index_params_s *)malloc(
		sizeof(struct comb_index_params_s));
	*params = (struct comb_index_params_s){
		.message_level = 1,
		.message_printer = (int (*)(void *, char const *, ...))fprintf,
		.message_context = (void *)stderr,
		.num_threads = 0,
		.mem_size = mem_estimate_free_size(),
		.command = comb_build_command_string(argc, argv),
		.command_base = strdup(base),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.k = 14
	};

	static struct option const opts_long[] = {
		/* global */
		{ "help", no_argument, NULL, 'h' },
		{ "version", no_argument, NULL, 'v' },
		{ "verbose", no_argument, NULL, 'V' },
		{ "threads", required_argument, NULL, 't' },
		{ "memory", required_argument, NULL, 'M' },
		{ "prefix", required_argument, NULL, 'p' },

		/* indexing params */
		{ "seedlength", required_argument, NULL, 'k' },
		{ 0 }

	};
	char *opts_short = comb_build_short_option_string(opts_long);

	/* optional arguments */
	int c, idx;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			/* global */
			case 'h': comb_index_print_help(); goto _comb_init_index_error_handler;
			case 'v': comb_print_version(); goto _comb_init_index_error_handler;
			case 'V': params->message_level = 3; break;
			case 't': params->num_threads = comb_atoi(optarg); break;
			case 'M': params->mem_size = mem_estimate_free_size(); break;
			case 'p': params->prefix = strdup(optarg); break;

			/* params */
			case 'k': params->k = comb_atoi(optarg); break;

			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
	}

	/* positional arguments */
	if(argc - optind == 1) {
		params->ref_name = strdup(argv[optind++]);
	} else {
		comb_print_invalid_args(argv, optind, argc);
		goto _comb_init_index_error_handler;
	}

	/* use ref file name as prefix if not specified */
	if(params->prefix == NULL) {
		params->prefix = strdup(params->ref_name);
	}

	free(opts_short);
	return(params);

_comb_init_index_error_handler:;
	free(opts_short);
	free(params);
	return(NULL);
}


/* `comb align' and `comb index' dispatcher */
/**
 * @fn comb_init
 */
static
struct comb_inst_s *comb_init(
	int argc,
	char *argv[])
{
	struct comb_inst_s *o = (struct comb_inst_s *)malloc(
		sizeof(struct comb_inst_s));
	if(o == NULL) { return(NULL); }

	/* argc must be greater than or equal to 2 */
	char const *base = argv[0];
	if(strcmp(argv[1], "index") == 0) {
		*o = (struct comb_inst_s){
			.params = (void *)comb_init_index(base, argc - 1, argv + 1),
			.process = (int (*)(void *))comb_index,
			.clean = (void (*)(void *))comb_index_clean
		};

	} else {
		if(strcmp(argv[1], "align") == 0) {
			argc--; argv++;
		}

		*o = (struct comb_inst_s){
			.params = (void *)comb_init_align(base, argc, argv),
			.process = (int (*)(void *))comb_align,
			.clean = (void (*)(void *))comb_align_clean
		};
	}

	if(o->params == NULL) {
		o->clean(o->params);
		free(o); o = NULL;
	}
	return(o);
}


/* `bwa' subcommand dispatcher */
static
char const *const comb_bwa_help_message =
	"\n"
	"    comb aligner (%s) bwa compatibility layer\n"
	"\n";

/**
 * @fn comb_bwa_print_help
 */
static _force_inline
void comb_bwa_print_help(
	void)
{
	fprintf(stderr, comb_bwa_help_message, COMB_VERSION_STRING);
	return;
}

/**
 * @fn comb_init_bwa_index
 */
static
struct comb_index_params_s *comb_init_bwa_index(
	char const *base,
	int argc,
	char *argv[])
{
	if(argc == 1) {
		comb_bwa_print_help();
		return(NULL);
	}

	/* build params object with default params */
	struct comb_index_params_s *params = (struct comb_index_params_s *)malloc(
		sizeof(struct comb_index_params_s));
	*params = (struct comb_index_params_s){
		.message_level = 1,
		.message_printer = (int (*)(void *, char const *, ...))fprintf,
		.message_context = (void *)stderr,
		.num_threads = 0,
		.mem_size = 0,
		.command = comb_build_command_string(argc, argv),
		.command_base = strdup(base),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.k = 14			/* the default k-mer length (14) is used */
	};

	static struct option const opts_long[] = {
		/* global params */
		{ NULL, required_argument, NULL, 'p' },

		/* indexing params is available */
		{ NULL, required_argument, NULL, 'k' },

		/* bwa-specific indexing params are ignored */
		{ NULL, required_argument, NULL, 'a' },
		{ NULL, required_argument, NULL, 'b' },
		{ NULL, no_argument, NULL, '6' }
	};
	char *opts_short = comb_build_short_option_string(opts_long);

	/* optional arguments */
	int c, idx;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			/* global */
			case 'p': params->prefix = strdup(optarg); break;

			/* params */
			case 'k': params->k = comb_atoi(optarg); break;

			/* ignored options */
			case 'a': break;
			case 'b': break;
			case '6': break;

			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
	}

	/* positional arguments */
	if(argc - optind == 1) {
		params->ref_name = strdup(argv[optind++]);
	} else {
		comb_print_invalid_args(argv, optind, argc);
		goto _comb_init_bwa_index_error_handler;
	}

	/* use ref file name as prefix if not specified */
	if(params->prefix == NULL) {
		params->prefix = strdup(params->ref_name);
	}

	free(opts_short);
	return(params);

_comb_init_bwa_index_error_handler:;
	free(params);
	free(opts_short);
	return(NULL);
}

/**
 * @fn comb_init_bwa_mem
 */
static
struct comb_align_params_s *comb_init_bwa_mem(
	char const *base,
	int argc,
	char *argv[])
{
	if(argc == 1) {
		comb_bwa_print_help();
		return(NULL);
	}

	/* build params object with default params */
	struct comb_align_params_s *params = (struct comb_align_params_s *)malloc(
		sizeof(struct comb_align_params_s));
	*params = (struct comb_align_params_s){
		.message_level = 3,
		.message_printer = (int (*)(void *, char const *, ...))fprintf,
		.message_context = (void *)stderr,

		.num_threads = 0,
		.mem_size = 0,
		.pool_size = 256,
		.command = comb_build_command_string(argc, argv),
		.command_base = strdup(base),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.out_name = strdup("-"),
		.out_format = AW_SAM,

		.k = 14,
		.kmer_cnt_thresh = 500,
		.overlap_thresh = 3,
		.gapless_thresh = 0,
		.xdrop = 100,
		.m = 1, .x = 4, .gi = 6, .ge = 1,
		.clip = 'S',
		.score_thresh = 30,
		.include_unmapped = 1
	};

	static struct option const opts_long[] = {
		/* global */
		{ NULL, required_argument, NULL, 't' },		/* -t INT        number of threads [1] */
		{ NULL, required_argument, NULL, 'v' },		/* -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3] */

		/* indexing params */
		{ NULL, required_argument, NULL, 'k' },		/* -k INT        minimum seed length [19] */

		/* filtering params */
		{ NULL, required_argument, NULL, 'c' },		/* -c INT        skip seeds with more than INT occurrences [500] */

		/* scoring params */
		{ NULL, required_argument, NULL, 'A' },		/* -A INT        score for a sequence match */
		{ NULL, required_argument, NULL, 'B' },		/* -B INT        penalty for a mismatch [4] */
		{ NULL, required_argument, NULL, 'O' },		/* -O INT[,INT]  gap open penalties for deletions and insertions [6,6] */
		{ NULL, required_argument, NULL, 'E' },		/* -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1] */
		{ NULL, required_argument, NULL, 'd' },		/* -d INT        off-diagonal X-dropoff [100] */

		/* reporting params */
		{ NULL, required_argument, NULL, 'T' },		/* -T INT        minimum score to output [30] */

		/* unimplemented options */
		{ NULL, required_argument, NULL, 'm' },
		{ NULL, no_argument, NULL, 'S' },
		{ NULL, no_argument, NULL, 'P' },
		{ NULL, required_argument, NULL, 'L' },
		{ NULL, required_argument, NULL, 'U' },
		{ NULL, required_argument, NULL, 'x' },
		{ NULL, no_argument, NULL, 'p' },
		{ NULL, required_argument, NULL, 'R' },
		{ NULL, required_argument, NULL, 'H' },
		{ NULL, no_argument, NULL, 'j' },
		{ NULL, required_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 'a' },
		{ NULL, no_argument, NULL, 'C' },
		{ NULL, no_argument, NULL, 'V' },
		{ NULL, no_argument, NULL, 'Y' },
		{ NULL, no_argument, NULL, 'M' },
		{ NULL, required_argument, NULL, 'I' },
		{ 0 }
	};
	char *opts_short = comb_build_short_option_string(opts_long);

	/* optional arguments */
	int c, idx;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			/* global */
			case 't': params->num_threads = comb_atoi(optarg); break;
			case 'v': params->message_level = comb_atoi(optarg); break;

			/* params */
			case 'k': params->k = comb_atoi(optarg); break;
			case 'c': params->kmer_cnt_thresh = comb_atoi(optarg); break;
			case 'A': params->m = comb_atoi(optarg); break;
			case 'B': params->x = comb_atoi(optarg); break;
			case 'O': params->gi = comb_atoi(optarg); break;
			case 'E': params->ge = comb_atoi(optarg); break;
			case 'd': params->xdrop = comb_atoi(optarg); break;
			case 'T': params->score_thresh = comb_atoi(optarg); break;

			/* unimplemented options */
			case 'm':
			case 'S':
			case 'P':
			case 'L':
			case 'U':
			case 'x':
			case 'p':
			case 'R':
			case 'H':
			case 'j':
			case 'h':
			case 'a':
			case 'C':
			case 'V':
			case 'Y':
			case 'M':
			case 'I':
				comb_print_error("unimplemented option `%c'.", c);
				break;
			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
	}

	/* restore default params */
	params->gapless_thresh = comb_init_align_default_gapless_thresh(params);

	/* positional arguments */
	int argcnt = argc - optind;
	switch(argcnt) {
		case 1:
			/* all-to-all alignment mode */
			params->ref_name = strdup(argv[optind]);
			params->query_name = strdup(argv[optind++]);
			break;
		case 2:
			/* mapping mode */
			params->ref_name = strdup(argv[optind++]);
			params->query_name = strdup(argv[optind++]);
			break;
		case 3:
			/* paired-end mapping mode */
			comb_print_error("paired-end mapping mode is not implemented.");
			goto _comb_init_bwa_mem_error_handler;
		default:
			comb_print_invalid_args(argv, optind, argc);
			goto _comb_init_bwa_mem_error_handler;
	}

	free(opts_short);
	return(params);

_comb_init_bwa_mem_error_handler:;
	free(params);
	free(opts_short);
	return(NULL);
}

/**
 * @fn comb_init_bwa_bt
 * @brief bwa-backtrack option parser
 */
static
struct comb_index_params_s *comb_init_bwa_bt(
	char const *base,
	int argc,
	char *argv[])
{
	if(argc == 1) {
		comb_bwa_print_help();
		return(NULL);
	}
	return(NULL);
}

/**
 * @fn comb_init_bwa
 */
static
struct comb_inst_s *comb_init_bwa(
	int argc,
	char *argv[])
{
	/* argc must be greater than or equal to 2 */
	if(argc == 2) {
		comb_bwa_print_help();
		return(NULL);
	}

	struct comb_inst_s *o = (struct comb_inst_s *)malloc(
		sizeof(struct comb_inst_s));
	if(o == NULL) { return(NULL); }

	/* here argc >= 3 is guaranteed */
	char const *base = argv[0];
	if(strcmp(argv[2], "index") == 0) {
		*o = (struct comb_inst_s){
			.params = (void *)comb_init_bwa_index(base, argc - 2, argv + 2),
			.process = (int (*)(void *))comb_index,
			.clean = (void (*)(void *))comb_index_clean
		};
	} else if(strcmp(argv[2], "mem") == 0) {
		*o = (struct comb_inst_s){
			.params = (void *)comb_init_bwa_mem(base, argc - 2, argv + 2),
			.process = (int (*)(void *))comb_align,
			.clean = (void (*)(void *))comb_align_clean
		};
	} else {
		*o = (struct comb_inst_s){
			.params = (void *)comb_init_bwa_bt(base, argc - 1, argv + 1),
			.process = (int (*)(void *))comb_align,
			.clean = (void (*)(void *))comb_align_clean
		};
	}

	if(o->params == NULL) {
		o->clean(o->params);
		free(o); o = NULL;
	}
	return(o);
}

#ifdef MAIN
/**
 * @fn main
 */
int main(int argc, char *argv[])
{
	/* print help when argc == 1 */
	if(argc == 1) {
		comb_align_print_help();
		return(1);
	}

	/* recurse to main if `comb comb' */
	if(strcmp(argv[1], "comb") == 0) {
		return(main(argc - 1, argv + 1));
	}

	/* here is guaranteed argc >= 2 */
	struct subc_map_s {
		char const *sub;
		struct comb_inst_s *(*init)(int argc, char *argv[]);
	};
	struct subc_map_s map[] = {
		{ "unittest", comb_init_unittest },
		{ "align", comb_init },
		{ "index", comb_init },
		{ "bwa", comb_init_bwa }
	};

	struct comb_inst_s *(*init)(int argc, char *argv[]) = comb_init;
	for(uint64_t i = 0; i < sizeof(map) / sizeof(struct subc_map_s); i++) {
		if(strcmp(argv[1], map[i].sub) == 0) {
			init = map[i].init; break;
		}
	}

	struct comb_inst_s *o = init(argc, argv);
	if(o == NULL) { return(1); }

	int ret = (o->process != NULL) ? o->process(o->params) : 1;
	if(o->clean != NULL) { o->clean(o->params); }
	free(o);
	return(ret);
}
#endif /* MAIN */

/**
 * end of comb.c
 */
