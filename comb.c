
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
#include "sassert.h"
#include "kvec.h"
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

char const *comb_help_message =
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
	"    $ comb [options] <reference> <query> <output>\n"
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
 * @struct comb_params_s
 */
struct comb_params_s {
	/* global */
	uint8_t help;
	uint8_t version;
	uint8_t pad1[2];
	int64_t num_threads;

	char *command;
	char *command_base;
	char *program_name;
	int64_t program_id;

	char *ref_name;
	char *query_name;
	char *out_name;

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
	int8_t pad2[3];
	int64_t score_thresh;
};

/**
 * @struct comb_worker_args_s
 */
struct comb_worker_args_s {
	struct comb_params_s const *params;
	struct sr_gref_s *r;
	ggsea_ctx_t *ctx;
	sr_t *ref;
	sr_t *query;
	aw_t *aw;
};

/**
 * @struct comb_worker_item_s
 */
struct comb_worker_item_s {
	lmm_t *lmm;
	struct sr_gref_s *q;
	struct ggsea_result_s *res;
};


/**
 * @fn comb_print_help
 */
#define COMB_PRINT_VERSION			0
#define COMB_PRINT_HELP 			1
#define COMB_PRINT_INVALID			2
static _force_inline
void comb_print_help(
	int level)
{
	if(level == COMB_PRINT_VERSION) {
		fprintf(stderr, "comb aligner (%s)\n", COMB_VERSION_STRING);
	} else if(level == COMB_PRINT_HELP) {
		fprintf(stderr, comb_help_message, COMB_VERSION_STRING);
	}
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
 * @fn comb_print_option_summary
 */
static _force_inline
void comb_print_option_summary(
	struct comb_params_s const *params)
{
	#define _p(...)		fprintf(stderr, __VA_ARGS__);
	_p("%s ", params->command_base);
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
	_p("\n");
	#undef _p
	return;
}

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

	for(po = opts; po->name != NULL; po++) { len++; }
	str = ps = (char *)malloc(2 * len + 1);
	for(po = opts; po->name != NULL; po++) {
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
	kvec_t(char) s;

	kv_init(s);
	for(int i = 0; i < argc; i++) {
		char const *p = argv[i];
		while(*p != '\0') {
			kv_push(s, *p);
			p++;
		}
		kv_push(s, ' ');
	}
	kv_at(s, kv_size(s) - 1) = '\0';
	return(kv_ptr(s));
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
 * @fn comb_calc_gapless_thresh
 */
static _force_inline
int64_t comb_calc_gapless_thresh(
	struct comb_params_s const *params)
{
	return((int64_t)(15.0 * params->x / (params->m + params->x)));
}

/**
 * @fn comb_calc_xdrop_thresh
 */
static _force_inline
int64_t comb_calc_xdrop_thresh(
	struct comb_params_s const *params)
{
	return(15 * (params->m + 2 * params->ge) + params->gi);
}

/**
 * @fn comb_calc_score_thresh
 */
static _force_inline
int64_t comb_calc_score_thresh(
	struct comb_params_s const *params)
{
	return(100 * params->m);
}

/**
 * @fn comb_parse_args
 */
static
struct comb_params_s *comb_parse_args(
	int argc,
	char *argv[])
{
	/* build params object with default params */
	struct comb_params_s *params = (struct comb_params_s *)malloc(
		sizeof(struct comb_params_s));
	*params = (struct comb_params_s){
		.num_threads = 0,
		.command = comb_build_command_string(argc, argv),
		.command_base = strdup(argv[0]),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.k = 14,
		.kmer_cnt_thresh = 30,
		.overlap_thresh = 3,
		.gapless_thresh = 0,
		.xdrop = 0,		/* default xdrop threshold is derived from scoring parameters */
		.m = 1, .x = 2, .gi = 2, .ge = 1,
		.clip = 'S',
		.score_thresh = 0
	};

	static struct option const opts_long[] = {
		/* global */
		{ "help", no_argument, NULL, 'h' },
		{ "version", no_argument, NULL, 'v' },
		{ "threads", required_argument, NULL, 't' },

		/* indexing params */
		{ "seedlength", required_argument, NULL, 'k' },

		/* filtering params */
		{ "repcnt", required_argument, NULL, 'r' },
		{ "depth", required_argument, NULL, 'd' },
		{ "popcnt", required_argument, NULL, 'f' },

		/* scoring params */
		{ "match", required_argument, NULL, 'a' },
		{ "mismatch", required_argument, NULL, 'b' },
		{ "gapopen", required_argument, NULL, 'p' },
		{ "gapextend", required_argument, NULL, 'q' },
		{ "xdrop", required_argument, NULL, 'x' },

		/* reporting params */
		{ "min", required_argument, NULL, 'm' },
		{ "clip", required_argument, NULL, 'c' }
	};
	char *opts_short = comb_build_short_option_string(opts_long);

	/* optional arguments */
	int c, idx;
	while((c = getopt_long(argc, argv, opts_short, opts_long, &idx)) != -1) {
		switch(c) {
			/* global */
			case 'h': comb_print_help(COMB_PRINT_HELP); exit(1);
			case 'v': comb_print_help(COMB_PRINT_VERSION); exit(1);
			case 't': params->num_threads = comb_atoi(optarg); break;

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
			
			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
	}

	/* restore default params */
	if(params->gapless_thresh == 0) {
		params->gapless_thresh = comb_calc_gapless_thresh(params);
	}
	if(params->xdrop == 0) {
		params->xdrop = comb_calc_xdrop_thresh(params);
	}
	if(params->score_thresh == 0) {
		params->score_thresh = comb_calc_score_thresh(params);
	}

	/* positional arguments */
	int argcnt = argc - optind;
	switch(argcnt) {
		case 1:
			/* all-to-all alignment mode, dump to stdout */
			params->ref_name = strdup(argv[optind]);
			params->query_name = strdup(argv[optind++]);
			params->out_name = strdup("-");
			break;
		case 2:
			/* all-to-all alignment mode, dump to file */
			params->ref_name = strdup(argv[optind]);
			params->query_name = strdup(argv[optind++]);
			params->out_name = strdup(argv[optind++]);
			break;
		case 3:
			/* mapping mode, dump to file */
			params->ref_name = strdup(argv[optind++]);
			params->query_name = strdup(argv[optind++]);
			params->out_name = strdup(argv[optind++]);
			break;
		default:
			comb_print_invalid_args(argv, optind, argc);
			break;
	}

	free(opts_short);
	return(params);
}

/**
 * @fn comb_clean_params
 */
static _force_inline
void comb_clean_params(
	struct comb_params_s *params)
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


/**
 * @fn comb_worker_init
 */
static _force_inline
struct comb_worker_args_s **comb_worker_init(
	struct comb_params_s const *params,
	ggsea_conf_t const *conf,
	sr_t *ref,
	sr_t *query,
	aw_t *aw)
{
	int64_t num_worker = MAX2(1, params->num_threads);
	struct comb_worker_args_s **w = (struct comb_worker_args_s **)malloc(
		  (num_worker + 1) * sizeof(struct comb_worker_args_s *)
		+ (num_worker + 1) * sizeof(struct comb_worker_args_s));
	struct comb_worker_args_s *base = (struct comb_worker_args_s *)(w + num_worker + 1);

	/* set pointers */
	for(int64_t i = 0; i < num_worker + 1; i++) {
		w[i] = &base[i];
	}

	/* init worker object */
	for(int64_t i = 0; i < num_worker; i++) {
		struct sr_gref_s *r = sr_get_index(ref);
		base[i] = (struct comb_worker_args_s){
			.params = params,
			.r = r,
			.ctx = ggsea_ctx_init(conf, r->gref),
			.ref = ref,
			.query = query,
			.aw = aw
		};
	}
	memset(&base[num_worker], 0, sizeof(struct comb_worker_args_s));
	return(w);
}

/**
 * @fn comb_worker_clean
 */
static _force_inline
void comb_worker_clean(
	struct comb_worker_args_s **w)
{
	if(w == NULL) { return; }

	for(struct comb_worker_args_s **p = w; (*p)->params != NULL; p++) {
		ggsea_ctx_clean((*p)->ctx);
		sr_gref_free((*p)->r);
		memset(*p, 0, sizeof(struct comb_worker_args_s));
	}
	free(w);
	return;
}

/**
 * @fn comb_source
 */
static
void *comb_source(
	void *arg)
{
	struct comb_worker_args_s *a = (struct comb_worker_args_s *)arg;

	/* get the next iterator */
	struct sr_gref_s *q = NULL;
	while((q = sr_get_iter(a->query)) != NULL) {
		debug("iterator fetched a(%p), q(%p), iter(%p)", a, q, q->iter);

		if(q->iter != NULL) {
			/* valid iterator was returned, create worker_item object */
			struct comb_worker_item_s *i = (struct comb_worker_item_s *)lmm_malloc(
				q->lmm, sizeof(struct comb_worker_item_s));

			/* init item */
			*i = (struct comb_worker_item_s){
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
 * @fn comb_worker
 */
static
void *comb_worker(
	void *arg,
	void *item)
{
	struct comb_worker_args_s *a = (struct comb_worker_args_s *)arg;
	struct comb_worker_item_s *i = (struct comb_worker_item_s *)item;

	/* do alignment */
	debug("align a(%p), i(%p), iter(%p)", a, i, i->q->iter);
	i->res = ggsea_align(a->ctx, i->q->gref, i->q->iter, i->lmm);
	return((void *)i);
}

/**
 * @fn comb_drain
 */
static
void comb_drain(
	void *arg,
	void *item)
{
	struct comb_worker_args_s *a = (struct comb_worker_args_s *)arg;
	struct comb_worker_item_s *i = (struct comb_worker_item_s *)item;

	/* append to result queue */
	aw_append_alignment(a->aw, i->res->ref, i->res->query, i->res->aln, i->res->cnt);

	/* cleanup */
	ggsea_aln_free(i->res);

	/* lmm must be freed before gref_free */
	debug("worker destroyed, ptr(%p)", i);
	struct sr_gref_s *q = i->q;
	lmm_free(i->lmm, (void *)i);
	sr_gref_free(q);
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

/**
 * @fn comb_align
 */
static _force_inline
int comb_align(
	struct comb_params_s const *params)
{
	int ret = 1;

	ggsea_conf_t *conf = NULL;
	sr_t *ref = NULL;
	sr_t *query = NULL;
	aw_t *aw = NULL;
	struct comb_worker_args_s **w = NULL;
	ptask_t *pt = NULL;

	#define comb_align_error(expr, ...) { \
		if(!(expr)) { \
			comb_print_error("" __VA_ARGS__); \
			goto _comb_align_error_handler; \
		} \
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
			.k = params->k,
			.seq_direction = SR_FW_ONLY,
			.num_threads = params->num_threads
		));
	comb_align_error(ref != NULL, "Failed to open reference file `%s'.", params->ref_name);

	/* build read pool */
	query = sr_init(params->query_name,
		SR_PARAMS(
			.k = params->k,
			.seq_direction = SR_FW_ONLY,
			// .num_threads = params->num_threads
		));
	comb_align_error(query != NULL, "Failed to open query file `%s'.", params->query_name);

	/* build alignment writer */
	struct sr_gref_s *r = sr_get_index(ref);
	comb_align_error(r != NULL, "Failed to build reference index.");
	aw = aw_init(params->out_name, r->gref,
		AW_PARAMS(
			.clip = params->clip,
			.program_id = params->program_id,
			.program_name = params->program_name,
			.command = params->command
		));
	sr_gref_free(r);
	comb_align_error(aw != NULL, "Failed to open output file `%s'.", params->out_name);

	/* initialize parallel task dispatcher */
	w = comb_worker_init(params, conf, ref, query, aw);
	pt = ptask_init(comb_worker, (void **)w, params->num_threads, 2048);
	comb_align_error(w != NULL && pt != NULL, "Failed to initialize parallel worker threads.");

	/* run tasks */
	ptask_stream(pt, comb_source, (void *)*w, comb_drain, (void *)*w, 512);

	/* destroy objects */
	ret = 0;
_comb_align_error_handler:;
	if(ref != query) {
		sr_clean(query); query = NULL;
	}
	sr_clean(ref); ref = NULL;
	aw_clean(aw); aw = NULL;
	comb_worker_clean(w); w = NULL;
	ggsea_conf_clean(conf); conf = NULL;
	ptask_clean(pt); pt = NULL;
	return(ret);
}

#ifdef MAIN
/**
 * @fn main
 */
int main(int argc, char *argv[])
{
	/* print help */
	if(argc == 1) {
		comb_print_help(COMB_PRINT_HELP);
		return(1);
	}

	/* option parsers */
	struct comb_params_s *(*parse_args)(int argc, char *argv[]) = comb_parse_args;
	if(strcmp(argv[1], "unittest") == 0) {
		return(unittest_main(argc, argv));
	}

	struct comb_params_s *params = parse_args(argc, argv);

	comb_print_option_summary(params);
	int ret = comb_align(params);

	comb_clean_params(params);
	return(ret);		
}
#endif /* MAIN */

/**
 * end of comb.c
 */
