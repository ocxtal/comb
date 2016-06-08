
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
#include "ptask/ptask.h"	/* parallel task dispatcher */
#include "sr/sr.h"			/* sequence reader */
#include "gref/gref.h"		/* graphical sequence indexer */
#include "ggsea/ggsea.h"	/* graph-to-graph seed-and-extend alignment */
#include "aw/aw.h"			/* alignment writer */
#include "sassert.h"
#include "kvec.h"
#include "log.h"
#include "lmm.h"


/* inline directive */
#define _force_inline				inline


/* constants */
#define COMB_VERSION_STRING			"0.0.1"

char const *comb_help_message =
	"\n"
	"    comb aligner (version %s)\n"
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
	int64_t popcnt_thresh;

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
	struct ggsea_result_s aln;
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
		fprintf(stderr, "comb aligner version %s\n", COMB_VERSION_STRING);
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
		.num_threads = 1,
		.command = comb_build_command_string(argc, argv),
		.program_name = strdup("comb"),
		.program_id = UNITTEST_UNIQUE_ID,
		.k = 14,
		.kmer_cnt_thresh = 30,
		.overlap_thresh = 3,
		.popcnt_thresh = 10,
		.xdrop = 60,
		.m = 1, .x = 1, .gi = 1, .ge = 1,
		.clip = 'S',
		.score_thresh = 10
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
			case 'f': params->popcnt_thresh = comb_atoi(optarg); break;
			case 'a': params->m = comb_atoi(optarg); break;
			case 'b': params->x = comb_atoi(optarg); break;
			case 'p': params->gi = comb_atoi(optarg); break;
			case 'q': params->ge = comb_atoi(optarg); break;
			case 'x': params->xdrop = comb_atoi(optarg); break;
			case 'm': params->score_thresh = comb_atoi(optarg); break;
			
			/* unknown option */
			default: comb_print_unknown_option(c); break;
		}
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
struct comb_worker_args_s *comb_worker_init(
	struct comb_params_s const *params,
	ggsea_conf_t const *conf,
	sr_t *ref,
	sr_t *query,
	aw_t *aw)
{
	struct comb_worker_args_s *w = (struct comb_worker_args_s *)malloc(
		(params->num_threads + 1) * sizeof(struct comb_worker_args_s));

	for(int64_t i = 0; i < params->num_threads; i++) {
		w[i] = (struct comb_worker_args_s){
			.params = params,
			.ctx = ggsea_ctx_init(conf, sr_get_index(ref)->gref),
			.ref = ref,
			.query = query,
			.aw = aw
		};
	}
	memset(&w[params->num_threads], 0, sizeof(struct comb_worker_args_s));
	return(w);
}

/**
 * @fn comb_worker_clean
 */
static _force_inline
void comb_worker_clean(
	struct comb_worker_args_s *w)
{
	if(w == NULL) { return; }

	for(struct comb_worker_args_s *p = w; p->params != NULL; p++) {
		ggsea_ctx_clean(p->ctx);
		memset(p, 0, sizeof(struct comb_worker_args_s));
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
	struct sr_gref_s *q = sr_get_iter(a->query);
	if(q == NULL) {
		/* reached the end */
		return(NULL);
	}

	/* create worker_item object */
	struct comb_worker_item_s *i = (struct comb_worker_item_s *)lmm_malloc(
		q->lmm, sizeof(struct comb_worker_item_s));

	/* init item */
	*i = (struct comb_worker_item_s){
		.lmm = q->lmm,
		.q = q
	};
	return((void *)i);
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
	i->aln = ggsea_align(a->ctx, i->q->gref);
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
	aw_append_alignment(a->aw, i->aln.ref, i->aln.query, i->aln.aln, i->aln.cnt);

	/* cleanup */
	ggsea_aln_free(i->aln);
	lmm_free(i->lmm, item);
	sr_gref_free(i->q);
	return;
}


/**
 * @fn comb_run
 */
static _force_inline
int comb_run(
	struct comb_params_s const *params)
{
	int ret = 1;

	ggsea_conf_t *conf = NULL;
	sr_t *ref = NULL;
	sr_t *query = NULL;
	aw_t *aw = NULL;
	struct comb_worker_args_s *w = NULL;
	ptask_t *pt = NULL;

	/* build ggsea configuration object */
	conf = ggsea_conf_init(GGSEA_PARAMS(
		.xdrop = params->xdrop,
		.score_matrix = GABA_SCORE_SIMPLE(params->m, params->x, params->gi, params->ge),
		.kmer_cnt_thresh = params->kmer_cnt_thresh,
		.overlap_thresh = params->overlap_thresh,
		.popcnt_thresh = params->popcnt_thresh,
		.score_thresh = params->score_thresh));
	if(conf == NULL) {
		fprintf(stderr, "[ERROR] Failed to create alignment configuration.\n");
		goto _comb_run_error_handler;
	}


	/* build reference sequence index */
	ref = sr_init(params->ref_name,
		SR_PARAMS(
			.k = params->k,
			.seq_direction = SR_FW_ONLY,
			.num_threads = params->num_threads
		));
	if(ref == NULL) {
		fprintf(stderr, "[ERROR] Failed to open reference file `%s'.\n", params->ref_name);
		goto _comb_run_error_handler;
	}


	/* build read pool */
	query = (strcmp(params->ref_name, params->query_name) == 0)
		? ref
		: sr_init(params->query_name,
			SR_PARAMS(
				.k = params->k,
				.seq_direction = SR_FW_ONLY,
				.num_threads = params->num_threads
			));
	if(query == NULL) {
		fprintf(stderr, "[ERROR] Failed to open query file `%s'.\n", params->query_name);
		goto _comb_run_error_handler;
	}


	/* build alignment writer */
	struct sr_gref_s *r = sr_get_index(ref);
	aw = aw_init(params->out_name, r->gref,
		AW_PARAMS(
			.clip = params->clip,
			.program_id = params->program_id,
			.program_name = params->program_name,
			.command = params->command
		));
	sr_gref_free(r);
	if(aw == NULL) {
		fprintf(stderr, "[ERROR] Failed to open output file `%s'.\n", params->out_name);
		goto _comb_run_error_handler;
	}


	/* initialize parallel task dispatcher */
	w = comb_worker_init(params, conf, ref, query, aw);
	pt = ptask_init(comb_worker, (void **)&w, params->num_threads, 2048);
	if(w == NULL || pt == NULL) {
		fprintf(stderr, "[ERROR] Failed to initialize parallel worker threads.\n");
		goto _comb_run_error_handler;
	}

	/* run tasks */
	ptask_stream(pt, comb_source, (void *)w, comb_drain, (void *)w, 512);

	/* destroy objects */
	ret = 0;
_comb_run_error_handler:;
	ggsea_conf_clean(conf); conf = NULL;
	if(ref != query) {
		sr_clean(query); query = NULL;
	}
	sr_clean(ref); ref = NULL;
	aw_clean(aw); aw = NULL;
	comb_worker_clean(w); w = NULL;
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
		argc--; argv++;
		return(unittest_main(argc, argv));
	}

	struct comb_params_s *params = parse_args(argc, argv);
	
	comb_print_option_summary(params);
	int ret = comb_run(params);

	comb_clean_params(params);
	return(ret);
}
#endif /* MAIN */

/**
 * end of comb.c
 */
