
/**
 * @file ggsea.h
 *
 * @brief Graph-to-Graph Seed-and-Extend Alignment
 *
 * @author Hajime Suzuki
 * @date 2016/4/12
 */
#ifndef _GGSEA_H_INCLUDED
#define _GGSEA_H_INCLUDED

#include <stdint.h>
#include "gref.h"
#include "gaba.h"


/* types */
/**
 * @type ggsea_conf_t
 */
typedef struct ggsea_conf_s ggsea_conf_t;

/**
 * @type ggsea_ctx_t
 */
typedef struct ggsea_ctx_s ggsea_ctx_t;

/**
 * @struct ggsea_params_s
 */
struct ggsea_params_s {
	/* local memory manager */
	void *lmm;

	/* score parameters */
	int16_t xdrop;
	gaba_score_t const *score_matrix;

	/* repetitive kmer filter */
	int64_t k;
	int64_t kmer_cnt_thresh;		/* kmer count threshold */

	/* overlap filter thresh */
	int64_t overlap_thresh;			/* depth */

	/* popcnt filter thresh */
	int64_t gapless_thresh;			/* threshold */

	/* score thresh */
	int64_t score_thresh;
};
typedef struct ggsea_params_s ggsea_params_t;

/**
 * @macro GGSEA_PARAMS
 * @brief utility macro for gaba_init, see example on header.
 */
#define GGSEA_PARAMS(...)			( &((struct ggsea_params_s const) { __VA_ARGS__ }) )

/**
 * @struct ggsea_result_s
 */
struct ggsea_result_s {
	void *reserved1;
	gref_idx_t const *ref;
	gref_acv_t const *query;
	struct gaba_alignment_s const *const *aln;
	uint32_t cnt;
	uint32_t reserved2;
};
typedef struct ggsea_result_s ggsea_result_t;


/* functions */

/**
 * @fn ggsea_conf_init
 * @brief create configuration object
 */
ggsea_conf_t *ggsea_conf_init(
	ggsea_params_t const *params);

/**
 * @fn ggsea_conf_clean
 * @brief cleanup configuration object
 */
void ggsea_conf_clean(
	ggsea_conf_t *conf);

/**
 * @fn ggsea_ctx_init
 * @brief initialize thread-local context with const reference index object
 */
ggsea_ctx_t *ggsea_ctx_init(
	ggsea_conf_t const *conf,
	gref_idx_t const *ref);

/**
 * @fn ggsea_ctx_clean
 * @brief cleanup thread-local context
 */
void ggsea_ctx_clean(
	ggsea_ctx_t *ctx);

/**
 * @fn ggsea_align
 * @brief do pairwise local alignment between reference in the context and given query
 */
ggsea_result_t *ggsea_align(
	ggsea_ctx_t *_ctx,
	gref_acv_t const *query,
	gref_iter_t *iter,
	void *lmm);

/**
 * @fn ggsea_aln_free
 */
void ggsea_aln_free(
	ggsea_result_t *aln);


#endif /* _GGSEA_H_INCLUDED */
/**
 * end of ggsea.h
 */
