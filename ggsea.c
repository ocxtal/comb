
/**
 * @file ggsea.c
 *
 * @brief Graph-to-Graph Seed-and-Extend Alignment
 *
 * @author Hajime Suzuki
 * @date 2016/4/12
 */

#define UNITTEST_UNIQUE_ID			10
#include "unittest.h"

#include <stdint.h>
#include "ggsea.h"
#include "hmap.h"
#include "psort.h"
#include "tree.h"
#include "gref.h"
#include "gaba.h"
#include "arch/arch.h"
#include "kvec.h"
#include "sassert.h"
#include "lmm.h"
#include "log.h"

#include "ngx_rbtree.h"

/* inline directive */
#define _force_inline				inline

/* constants */
#define MARGIN_SEQ_SIZE				( 64 )
#define MARGIN_SEQ_OFFSET			( 16 )
#define MARGIN_SEQ_LEN				( 32 )

/* max and min */
#define MAX2(x,y) 		( (x) > (y) ? (x) : (y) )
#define MAX3(x,y,z) 	( MAX2(x, MAX2(y, z)) )
#define MAX4(w,x,y,z) 	( MAX2(MAX2(w, x), MAX2(y, z)) )

#define MIN2(x,y) 		( (x) < (y) ? (x) : (y) )
#define MIN3(x,y,z) 	( MIN2(x, MIN2(y, z)) )
#define MIN4(w,x,y,z) 	( MIN2(MIN2(w, x), MIN2(y, z)) )


/* static assert */
_static_assert(offsetof(struct gref_section_s, gid) == 0);

/* structs */
/**
 * @struct ggsea_conf_s
 * @brief configuration container
 */
struct ggsea_conf_s {
	/* alignment context */
	gaba_t *gaba;

	/* params */
	uint64_t init_rep_hash_size;
	uint64_t max_rep_vec_size;		/* max kmer vector size */
	uint64_t overlap_width;			/* overlap filter width */
	uint64_t res_lmm_size;			/* result memory manager size */
	struct ggsea_params_s params;
};

/**
 * @struct rep_seed_s
 */
struct rep_seed_s {
	hmap_header_t header;
	uint64_t vec_size;
	kvec_t(struct gref_gid_pos_s) rv;
	kvec_t(struct gref_gid_pos_s) qv;
};

/**
 * @struct rtree_node_s
 */
#if 0
struct rtree_node_s {
	rbtree_node_t h;		/* (40) */
	uint32_t prev_qpos;		/* q-coordinate at the previous r-index update */
	uint32_t path_qpos;		/* q-coordinate at the previous path adjustment */
	uint32_t path_ridx;		/* reverse index of the path string after the previous path adjustment */
	uint32_t path_rem;		/* rem length from the tail */
	uint64_t const *ptail;	/* path tail of the current section */

	struct gaba_alignment_s const *aln;
	uint32_t sidx;
	uint32_t pad1;

	struct qtree_node_s *qhead;
	uint32_t pad2[2];
};
#endif
struct rtree_node_s {
	rbtree_node_t h;		/* (40) */
	uint32_t prev_qpos;		/* q-coordinate on the previous r-index update */
	uint32_t path_qpos;		/* q-coordinate on the previous path adjustment */
	uint32_t qlim;			/* tail q-coordinate */
	uint32_t pad1;

	int64_t path_ridx;		/* reverse index of the path string after the previous path adjustment */
	uint64_t const *ptail;	/* path tail of the current section */

	struct gaba_alignment_s const *aln;
	uint32_t sidx;
	uint32_t pad2;

	struct qtree_node_s *qhead;
};
_static_assert(sizeof(struct rtree_node_s) == 96);

/**
 * @struct qtree_node_s
 */
struct qtree_node_s {
	rbtree_node_t h;		/* key: (id, pos) pair on query side */	
	struct gaba_alignment_s const *aln;
	uint32_t sidx;			/* section array index */
	uint32_t res_id;

	struct qtree_node_s *next;
};
_static_assert(sizeof(struct qtree_node_s) == 64);

/**
 * @struct dp_fill_pair_s
 */
struct dp_fill_pair_s {
	gaba_fill_t const *fw;
	gaba_fill_t const *rv;
};
_static_assert(sizeof(struct dp_fill_pair_s) == 16);

/**
 * @union ggsea_gid_pos_u
 */
union ggsea_gid_pos_u {
	struct gref_gid_pos_s p;
	uint64_t u;
};
_static_assert(sizeof(struct gref_gid_pos_s) == sizeof(uint64_t));
#define _cast_u(x)	( ((union ggsea_gid_pos_u){ .p = (x) }).u )
#define _cast_p(x)	( ((union ggsea_gid_pos_u){ .u = (x) }).p )

/**
 * @struct dp_front_s
 * @brief segment info container (to push into heapqueue)
 */
struct dp_front_s {
	int64_t psum;
	gaba_fill_t const *fill;
	uint32_t rgid;
	uint32_t qgid;
};
_static_assert(sizeof(struct dp_front_s) == 24);

/**
 * @struct ggsea_ctx_s
 */
struct ggsea_ctx_s {
	/* memory manager */
	lmm_t *lmm;

	/* constants */
	struct ggsea_conf_s conf;	/* copy of conf */

	/* current seqeunce info */
	gref_idx_t const *r;
	gref_acv_t const *q;

	/* repetitive kmer container */
	hmap_t *rep;

	/* seed filters */
	rbtree_t *rtree;
	rbtree_t *qtree;

	/* dp context */
	gaba_dp_t *dp;
	kvec_t(struct dp_front_s) queue;	/* segment queue */
	uint8_t *margin;
	struct gref_section_s fw_margin, rv_margin;

	/* result vector */
	lmm_t *res_lmm;
	kvec_t(struct gaba_alignment_s const *) aln;
};


/* conf init / clean */
/**
 * @fn ggsea_conf_clean
 */
void ggsea_conf_clean(
	ggsea_conf_t *_conf)
{
	struct ggsea_conf_s *conf = (struct ggsea_conf_s *)_conf;

	if(conf != NULL) {
		gaba_clean(conf->gaba); conf->gaba = NULL;
		free(conf);
	}
	return;
}

/**
 * @fn ggsea_conf_init
 */
ggsea_conf_t *ggsea_conf_init(
	ggsea_params_t const *params)
{
	/* restore default params */
	struct ggsea_params_s const default_params = { 0 };
	struct ggsea_params_s p = (params == NULL) ? default_params : *params;

	/* restore defaults */
	#define restore(param, def)		{ (param) = ((uint64_t)(param) == 0) ? (def) : (param); }

	restore(p.k, 4);
	restore(p.kmer_cnt_thresh, 100);
	restore(p.overlap_thresh, 3);
	restore(p.gapless_thresh, 10);
	restore(p.score_thresh, 0);			/* score positive */

	#undef restore

	/* malloc mem (params->lmm is ignored) */
	struct ggsea_conf_s *conf = (struct ggsea_conf_s *)malloc(
		sizeof(struct ggsea_conf_s));
	if(conf == NULL) {
		return(NULL);
	}
	memset(conf, 0, sizeof(struct ggsea_conf_s));

	/* init alignment context */
	conf->gaba = gaba_init(GABA_PARAMS(
		.filter_thresh = p.gapless_thresh,
		.xdrop = p.xdrop,
		.score_matrix = p.score_matrix));
	if(conf->gaba == NULL) {
		free(conf);
		return(NULL);
	}

	/* store constants */
	conf->init_rep_hash_size = 1024;
	conf->max_rep_vec_size = 128;
	conf->overlap_width = 48;
	conf->res_lmm_size = 16 * 1024 * 1024;		/* 16MB */
	conf->params = p;

	return((ggsea_conf_t *)conf);
}


/* context init / clean */
/**
 * @fn ggsea_ctx_clean
 */
void ggsea_ctx_clean(
	ggsea_ctx_t *ctx)
{
	if(ctx != NULL) {
		/* destroy repetitive kmer vectors */
		if(ctx->rep != NULL) {
			int64_t hcnt = hmap_get_count(ctx->rep);
			for(int64_t i = 0; i < hcnt; i++) {
				struct rep_seed_s *c = hmap_get_object(ctx->rep, i);
				debug("free rv(%p), qv(%p)", kv_ptr(c->rv), kv_ptr(c->qv));
				kv_destroy(c->rv);
				kv_destroy(c->qv);
			}
		}
		hmap_clean(ctx->rep); ctx->rep = NULL;

		/* destroy seed filter tree */
		rbtree_clean(ctx->rtree); ctx->rtree = NULL;
		rbtree_clean(ctx->qtree); ctx->qtree = NULL;

		/* destroy tree traversing queue */
		kv_hq_destroy(ctx->queue);

		/* margin sequence */
		free(ctx->margin); ctx->margin = NULL;

		/* dp context */
		gaba_dp_clean(ctx->dp); ctx->dp = NULL;

		/* ggsea context */
		free(ctx);
	}
	return;
}

/**
 * @fn ggsea_ctx_init
 */
ggsea_ctx_t *ggsea_ctx_init(
	ggsea_conf_t const *conf,
	gref_idx_t const *ref)
{
	/* malloc mem */
	struct ggsea_ctx_s *ctx = (struct ggsea_ctx_s *)malloc(
		sizeof(struct ggsea_ctx_s));
	if(ctx == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	memset(ctx, 0, sizeof(struct ggsea_ctx_s));

	/* copy conf */
	ctx->conf = *conf;

	/* clear ref and query pointers */
	if(ref == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	ctx->r = ref;
	ctx->q = NULL;

	/* init repetitive kmer filter */
	ctx->rep = hmap_init(
		sizeof(struct rep_seed_s),
		HMAP_PARAMS(.hmap_size = conf->init_rep_hash_size));
	if(ctx->rep == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}

	/* init overlap tree */
	ctx->rtree = rbtree_init(sizeof(struct rtree_node_s), NULL);
	ctx->qtree = rbtree_init(sizeof(struct qtree_node_s), NULL);

	/* init queue */
	kv_hq_init(ctx->queue);
	if(kv_ptr(ctx->queue) == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	debug("init, hq_size(%llu)", kv_hq_size(ctx->queue));

	/* init margin seq */
	uint64_t margin_size = 2 * sizeof(uint8_t) * (MARGIN_SEQ_SIZE + 32);
	if((ctx->margin = (uint8_t *)malloc(margin_size)) == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	memset(ctx->margin, 0, margin_size);

	/* init forward margin section */
	ctx->fw_margin = (struct gref_section_s){
		.gid = 0xfffc,
		.len = MARGIN_SEQ_LEN,
		.base = &ctx->margin[32]
	};

	/* init reverse margin section */
	ctx->rv_margin = (struct gref_section_s){
		.gid = 0xfffd,
		.len = MARGIN_SEQ_LEN,
		.base = &ctx->margin[MARGIN_SEQ_SIZE + 32]
	};

	/* init dp context (with invalid seq pair) */
	ctx->dp = gaba_dp_init(conf->gaba, NULL, NULL);
	if(ctx->dp == NULL) {
		goto _ggsea_ctx_init_error_handler;
	}
	return(ctx);

_ggsea_ctx_init_error_handler:;
	if(ctx != NULL) {
		ggsea_ctx_clean(ctx);
	}
	return(NULL);
}

/**
 * @fn ggsea_flush
 */
static _force_inline
void ggsea_flush(
	struct ggsea_ctx_s *ctx,
	gref_acv_t const *query,
	lmm_t *lmm)
{
	debug("flush called");

	/* set sequence info */
	ctx->q = query;

	/* flush hashmap */
	int64_t hcnt = hmap_get_count(ctx->rep);
	for(int64_t i = 0; i < hcnt; i++) {
		struct rep_seed_s *c = hmap_get_object(ctx->rep, i);
		debug("free rv(%p), qv(%p)", kv_ptr(c->rv), kv_ptr(c->qv));

		kv_destroy(c->rv);
		kv_destroy(c->qv);
	}
	hmap_flush(ctx->rep);

	/* flush tree */
	rbtree_flush(ctx->rtree);
	rbtree_flush(ctx->qtree);

	/* flush queues */
	kv_hq_clear(ctx->queue);

	/* flush dp context for the new read */
	debug("rlim(%p), qlim(%p)", gref_get_lim(ctx->r), gref_get_lim(ctx->q));
	gaba_dp_flush(ctx->dp, gref_get_lim(ctx->r), gref_get_lim(ctx->q));

	/* flush result vector */
	ctx->res_lmm = lmm;
	lmm_kv_init(ctx->res_lmm, ctx->aln);
	debug("flushed, aln(%p), lmm(%p), lim(%p)", lmm_kv_ptr(ctx->aln), ctx->res_lmm,
		(ctx->res_lmm != NULL) ? ctx->res_lmm->lim : NULL);
	return;
}


/* repetitive kmer filters */
/**
 * @fn rep_dedup_pos
 */
static _force_inline
int64_t rep_dedup_pos(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s *ptr,
	int64_t size)
{
	/* sort */
	psort_full((void *)ptr, size, sizeof(struct gref_gid_pos_s), 0);

	/* dedup */
	uint64_t *p = (uint64_t *)ptr;
	int64_t cnt = 0;
	for(int64_t i = 0; i < size; i++) {
		if(p[cnt] == p[i]) { continue; }
		p[++cnt] = p[i];
	}
	return(cnt + 1);
}

/**
 * @fn rep_save_pos
 */
static _force_inline
void rep_save_pos(
	struct ggsea_ctx_s *ctx,
	uint64_t kmer,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* add kmer to the hashmap */
	uint32_t id = hmap_get_id(ctx->rep,
		(char const *)&kmer, sizeof(uint64_t));
	struct rep_seed_s *c = hmap_get_object(ctx->rep, id);

	debug("id(%u), count(%u), diff(%d)", id, hmap_get_count(ctx->rep), id - hmap_get_count(ctx->rep) + 1);
	if(id == (hmap_get_count(ctx->rep) - 1)) {
		c->vec_size = ctx->conf.max_rep_vec_size;
		if(kv_ptr(c->rv) == NULL) { kv_init(c->rv); }
		if(kv_ptr(c->qv) == NULL) { kv_init(c->qv); }

		debug("malloc rv(%p), qv(%p)", kv_ptr(c->rv), kv_ptr(c->qv));
	}

	debug("save repetitive kmer(%llx), r(%u, %u), q(%u, %u), rv(%p, %llu), qv(%p, %llu)",
		kmer, rpos.gid, rpos.pos, qpos.gid, qpos.pos,
		kv_ptr(c->rv), kv_size(c->rv), kv_ptr(c->qv), kv_size(c->qv));

	/* add rpos */
	kv_push(c->rv, rpos);
	if(kv_size(c->rv) > c->vec_size) {
		kv_size(c->rv) = rep_dedup_pos(ctx, kv_ptr(c->rv), kv_size(c->rv));
	}

	kv_push(c->qv, qpos);
	if(kv_size(c->qv) > c->vec_size) {
		kv_size(c->qv) = rep_dedup_pos(ctx, kv_ptr(c->qv), kv_size(c->qv));
	}

	/* expand vector if spilled */
	if(kv_size(c->rv) > c->vec_size && kv_size(c->qv) > c->vec_size) {
		c->vec_size *= 2;
	}
	return;
}


/* graph traverse functions */
/**
 * @macro _rup, _qup
 * @brief extract update flag
 */
#define _rup(_fill)			( (_fill)->status & GABA_STATUS_UPDATE_A )
#define _qup(_fill)			( (_fill)->status & GABA_STATUS_UPDATE_B )
#define _rqup(_fill)		( (_fill)->status & (GABA_STATUS_UPDATE_A | GABA_STATUS_UPDATE_B) )
#define _term(_fill)		( (_fill)->status & GABA_STATUS_TERM )

/**
 * @fn dp_extend_leaf
 */
gaba_fill_t const *dp_extend_leaf(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec,
	uint32_t trigger_mask)
{
	debug("fill leaf");

	/* fill loop */
	while(1) {
		fill = gaba_dp_fill(ctx->dp, fill,
			(struct gaba_section_s *)rsec,
			(struct gaba_section_s *)qsec);
		debug("status(%x), max(%lld), r(%u), q(%u)",
			fill->status, fill->max, rsec->gid, qsec->gid);
		max = (fill->max > max->max) ? fill : max;

		if(fill->status & trigger_mask) { break; }
		if(_rup(fill) != 0) {
			trigger_mask |= GABA_STATUS_UPDATE_A;
			rsec = &ctx->fw_margin;
		}
		if(_qup(fill) != 0) {
			trigger_mask |= GABA_STATUS_UPDATE_B;
			qsec = &ctx->rv_margin;
		}
	}
	return(max);
}

/**
 * @fn dp_extend_update_queue
 */
gaba_fill_t const *dp_extend_update_queue(
	struct ggsea_ctx_s *ctx,
	gaba_fill_t const *fill,
	gaba_fill_t const *max,
	struct gref_section_s const *rsec,
	struct gref_section_s const *qsec)
{
	/* retrieve link info if update flag is set */
	struct gref_link_s rlink = { (uint32_t const *)rsec, 1 };
	struct gref_link_s qlink = { (uint32_t const *)qsec, 1 };
	uint32_t trigger_mask = GABA_STATUS_TERM;

	if(_rup(fill) != 0) {
		rlink = gref_get_link(ctx->r, rsec->gid);
		
		/* replace leaf link if no successors found */
		if(rlink.len == 0) {
			trigger_mask |= GABA_STATUS_UPDATE_A;
			rsec = &ctx->fw_margin;
		}
	}
	if(_qup(fill) != 0) {
		qlink = gref_get_link(ctx->q, qsec->gid);

		/* replace leaf link if no successors found */
		if(qlink.len == 0) {
			trigger_mask |= GABA_STATUS_UPDATE_B;
			qsec = &ctx->rv_margin;
		}
	}

	if(trigger_mask != 0) {
		return(dp_extend_leaf(ctx, fill, max, rsec, qsec, trigger_mask));
	}

	/* push section pairs */
	for(int64_t i = 0; i < rlink.len; i++) {
		for(int64_t j = 0; j < qlink.len; j++) {
			debug("push queue, fill(%p), psum(%lld), r(%u), q(%u)",
				fill, fill->psum, rlink.gid_arr[i], qlink.gid_arr[j]);

			kv_hq_push(ctx->queue, ((struct dp_front_s){
				.psum = (int64_t)fill->psum,
				.fill = fill,
				.rgid = rlink.gid_arr[i],
				.qgid = qlink.gid_arr[j]
			}));
		}
	}
	return(max);
}

/**
 * @fn dp_extend_intl
 */
static _force_inline
gaba_fill_t const *dp_extend_intl(
	struct ggsea_ctx_s *ctx,
	struct gref_section_s const *rsec,
	uint32_t rpos,
	struct gref_section_s const *qsec,
	uint32_t qpos)
{
	debug("seed: r(%u, %u), q(%u, %u)", rsec->gid, rpos, qsec->gid, qpos);

	/* flush queue */
	kv_hq_clear(ctx->queue);

	/* fill the first section */
	gaba_fill_t const *max = NULL;
	gaba_fill_t const *fill = max = gaba_dp_fill_root(ctx->dp,
		(struct gaba_section_s *)rsec, rpos,
		(struct gaba_section_s *)qsec, qpos);

	debug("root: status(%x), max(%lld), r(%u, %u), q(%u, %u)",
		fill->status, fill->max, rsec->gid, rpos, qsec->gid, qpos);

	/* check xdrop term */
	if((fill->status & GABA_STATUS_TERM) != 0) {
		return(max);
	}

	/* update first joint */
	max = dp_extend_update_queue(ctx, fill, max, rsec, qsec);

	/* loop */
	while(kv_hq_size(ctx->queue) > 0) {
		struct dp_front_s seg = kv_hq_pop(ctx->queue);
		debug("pop queue, fill(%p), psum(%lld), r(%u), q(%u)",
			seg.fill, seg.psum, seg.rgid, seg.qgid);

		if(seg.fill == NULL) {
			debug("fill == NULL (unexpected NULL pointer detected)");
			break;
		}

		/**
		 * Lazy merge detection comes here.
		 * tails at the head of the queue is scaned and overlapping pairs are
		 * detected in the lazy strategy. the detected pairs will be merged
		 * by ggsea_merge_tails function.
		 */

		/* extend */
		rsec = gref_get_section(ctx->r, seg.rgid);
		qsec = gref_get_section(ctx->q, seg.qgid);
		gaba_fill_t *fill = gaba_dp_fill(ctx->dp, seg.fill,
			(struct gaba_section_s const *)rsec,
			(struct gaba_section_s const *)qsec);

		/* update max */
		debug("check max, max(%lld), prev_max(%lld)", fill->max, max->max);
		max = (fill->max > max->max) ? fill : max;

		/* check xdrop term */
		if((fill->status & GABA_STATUS_TERM) == 0) {
			max = dp_extend_update_queue(ctx, fill, max, rsec, qsec);
			debug("queue updated, max(%lld)", max->max);
		}
	}

	debug("extend finished, max(%lld)", max->max);
	return(max);
}

/**
 * @fn dp_extend
 */
static _force_inline
struct dp_fill_pair_s dp_extend(
	struct ggsea_ctx_s *ctx,
	struct gaba_path_section_s *sec,
	int64_t len)
{
	/* forward section */
	debug("forward extend");
	struct gref_section_s const *rfsec = gref_get_section(ctx->r, sec[len - 1].aid);
	struct gref_section_s const *qfsec = gref_get_section(ctx->q, sec[len - 1].bid);
	gaba_fill_t const *fw_max = dp_extend_intl(ctx,
		rfsec, sec[len - 1].apos + sec[len - 1].alen,
		qfsec, sec[len - 1].bpos + sec[len - 1].blen);

	/* reverse section */
	debug("reverse extend");
	struct gref_section_s const *rrsec = gref_get_section(ctx->r, gref_rev_gid(sec[0].aid));
	struct gref_section_s const *qrsec = gref_get_section(ctx->q, gref_rev_gid(sec[0].bid));
	gaba_fill_t const *rv_max = dp_extend_intl(ctx,
		rrsec, rrsec->len - sec[0].apos,
		qrsec, qrsec->len - sec[0].bpos);

	debug("fw_max(%lld), rv_max(%lld), max(%lld)",
		fw_max->max, rv_max->max, fw_max->max + rv_max->max);
	/* return max pair */
	return((struct dp_fill_pair_s){
		.fw = fw_max,
		.rv = rv_max
	});
}

/**
 * @fn dp_expand_pos
 */
static _force_inline
int64_t dp_expand_pos(
	struct ggsea_ctx_s *ctx,
	struct gaba_path_section_s *sec,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* init rem length */
	int64_t const mask = GREF_K_MAX - 1;
	int64_t const ofs = ctx->conf.params.k - 1;
	int64_t rem = ofs;

	/* save gids of the first section */
	uint32_t rgid = sec->aid = rpos.gid;
	uint32_t qgid = sec->bid = qpos.gid;

	debug("rid(%x), qid(%x), k(%lld)", rgid, qgid, ctx->conf.params.k);

	/* calc remaining length */
	int64_t rlen = gref_get_section(ctx->r, rgid)->len;
	int64_t qlen = gref_get_section(ctx->q, qgid)->len;
	int64_t rrem = (int64_t)rpos.pos - rlen + ofs;
	int64_t qrem = (int64_t)qpos.pos - qlen + ofs;

	/* save apos and bpos */
	int64_t ridx = sec->apos = rlen - ofs + (((rrem < 0) ? -1 : mask) & rrem);
	int64_t qidx = sec->bpos = qlen - ofs + (((qrem < 0) ? -1 : mask) & qrem);
	debug("ridx(%lld), qidx(%lld)", ridx, qidx);

	/* calc reverse index */
	int64_t rridx = rlen - ridx;
	int64_t qridx = qlen - qidx;
	debug("rlen(%lld), qlen(%lld), rpos(%u), qpos(%u)", rlen, qlen, rpos.pos, qpos.pos);

	/* save lengths */
	int64_t len = sec->alen = sec->blen = MIN3(rridx, qridx, ctx->conf.params.k);
	sec->ppos = 0;
	// sec->plen = 2 * len;

	debug("rem(%lld), len(%lld), rridx(%lld), qridx(%lld)", rem, len, rridx, qridx);

	/* return if ptr reached the end of the path */
	if((rem -= len) <= 0) {
		return(1);
	}

	/* retrieve path info */
	uint32_t rlink = rrem>>GREF_K_MAX_BASE;
	uint32_t qlink = qrem>>GREF_K_MAX_BASE;
	struct gaba_path_section_s const *sec_base = sec;

	do {
		if((rridx -= len) <= 0) {
			rgid = gref_get_link(ctx->r, rgid).gid_arr[rlink & mask];
			rridx = gref_get_section(ctx->r, rgid)->len;
			rlink >>= GREF_K_MAX_BASE;
		}
		if((qridx -= len) <= 0) {
			qgid = gref_get_link(ctx->q, qgid).gid_arr[qlink & mask];
			qridx = gref_get_section(ctx->q, qgid)->len;
			qlink >>= GREF_K_MAX_BASE;
		}
		len = MIN2(rridx, qridx);
		*++sec = (struct gaba_path_section_s){
			.aid = rgid,
			.bid = qgid,
			.apos = 0,
			.bpos = 0,
			.alen = len,
			.blen = len,
			.ppos = 0
			/* .plen = 2 * len */
		};
	} while((rem -= len) > 0);
	return(sec - sec_base + 1);
}

/**
 * @fn dp_extend_seed
 */
static _force_inline
struct gaba_alignment_s const *dp_extend_seed(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* save stack */
	gaba_stack_t const *stack = gaba_dp_save_stack(ctx->dp);

	debug("extend seed rpos(%llx), qpos(%llx)", _cast_u(rpos), _cast_u(qpos));

	/* expand path on the seed */
	struct gaba_path_section_s sec[ctx->conf.params.k];
	int64_t len = dp_expand_pos(ctx, sec, rpos, qpos);

	/* extend */
	struct dp_fill_pair_s pair = dp_extend(ctx, sec, len);
	debug("fw_max(%lld), rv_max(%lld)", pair.fw->max, pair.rv->max);
	if(pair.fw->max + pair.rv->max <= ctx->conf.params.score_thresh) {
		debug("stack flushed, score(%lld, %lld)", pair.fw->max + pair.rv->max, ctx->conf.params.score_thresh);
		gaba_dp_flush_stack(ctx->dp, stack);
		return(NULL);
	}

	/* traceback */
	struct gaba_alignment_s const *aln = gaba_dp_trace(
		ctx->dp, pair.fw, pair.rv,
		GABA_TRACE_PARAMS(
			.lmm = ctx->res_lmm,
			.sec = sec,
			.slen = len,
			.k = ctx->conf.params.k
		));

	debug("trace finished, score(%lld), plen(%llu), slen(%u)", aln->score, aln->path->len, aln->slen);
	return(aln);
}


/* result containers and seed filtering */
/**
 * @fn load_u64
 */
static inline
uint64_t load_u64(
	uint64_t const *ptr,
	int64_t ridx)
{
	uint64_t const mask = 0xaaaaaaaaaaaaaaaa;
	int64_t idx = (-ridx)>>6;
	uint64_t farr = (idx <= 0) ? ptr[idx] : mask;
	uint64_t larr = (idx < 0) ? ptr[idx + 1] : mask;

	int64_t rem = (-ridx) & 63;
	return((farr>>rem) | ((larr<<(63 - rem))<<1));
}

/**
 * @fn rtree_append_result
 */
static _force_inline
struct rtree_node_s *rtree_append_result(
	struct ggsea_ctx_s *ctx,
	struct qtree_node_s *qhead,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	/* create rtree node */
	struct rtree_node_s *rn = (struct rtree_node_s *)
		rbtree_create_node(ctx->rtree);

	/* aln->path->array is always aligned on 8byte boundary */
	uint64_t const *path = (uint64_t const *)aln->path->array;

	/* root section */
	struct gaba_path_section_s const *rsec = &aln->sec[aln->rsidx];
	uint32_t plen = gaba_plen(rsec);
	*rn = (struct rtree_node_s){
		.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
			.gid = rsec->aid,
			.pos = aln->rapos + ctx->conf.overlap_width
		})),
		.prev_qpos = qpos.pos,
		.path_qpos = qpos.pos,
		.qlim = rsec->bpos + rsec->blen,

		.path_ridx = (plen - aln->rppos) & ~(64 - 1),
		.ptail = path + ((rsec->ppos + plen)>>6),

		.aln = aln,
		.sidx = aln->rsidx,

		.qhead = qhead
	};

	debug("append result, rn(%p), a(%u, %u), b(%u, %u)",
		rn, rsec->apos, rsec->alen, rsec->bpos, rsec->blen);
	rbtree_insert(ctx->rtree, (rbtree_node_t *)rn);
	return((struct rtree_node_s *)rbtree_right(ctx->rtree, (rbtree_node_t *)rn));
}

/**
 * @fn rtree_update
 * @brief update pos for the next r-iteration.
 */
static _force_inline
uint64_t rtree_update(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gref_gid_pos_s qpos)
{
	/* advance pos */
	int64_t inc = qpos.pos - rn->prev_qpos;

	debug("update rtree, rpos(%llu), inc(%lld), prev_qpos(%x), path_qpos(%x), qpos(%x)",
		0xffffffff & rn->h.key + inc, inc, rn->prev_qpos, rn->path_qpos, qpos.pos);

	rn->h.key += inc;
	rn->prev_qpos = qpos.pos;

	if(qpos.pos - rn->path_qpos >= 32) {
		/* adjust path */
		int64_t q = qpos.pos - rn->path_qpos;
		int64_t rpos = rn->h.key;
		int64_t ridx = rn->path_ridx;

		debug("q(%lld), ridx(%lld)", q, ridx);

		while(q > 0) {

			/* calculate advancing length on q side */
			int64_t qlen = MIN2(q, 32);

			/* load path array */
			uint64_t path_array = load_u64(rn->ptail, ridx);

			/* count vertical elements */
			int64_t dcnt = popcnt(path_array<<(64 - 2*qlen));

			debug("adjust path, qlen(%lld), dcnt(%lld), rpos(%llu), diff(%lld), qrem(%lld), ridx(%lld), path_array(%llx)",
				qlen, dcnt, 0xffffffff & rpos + 2*(qlen - dcnt), 2*(qlen - dcnt), q, ridx, path_array);

			rpos += 2 * (qlen - dcnt);
			q -= qlen;
			ridx -= 2*qlen;
		}

		/* write back indices */
		rn->h.key = rpos;
		rn->path_qpos = qpos.pos;
		rn->path_ridx = ridx;
	}
	return(rn->h.key);
}

/**
 * @fn rtree_advance
 */
static _force_inline
struct rtree_node_s *rtree_advance(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gref_gid_pos_s qpos)
{
	/* fetch next */
	struct rtree_node_s *next = (struct rtree_node_s *)rbtree_right(
		ctx->rtree, (rbtree_node_t *)rn);
	debug("fetched rnode, rn(%p, %lld), next(%p, %lld)",
		rn, rn->h.key,
		next, (next != NULL) ? next->h.key : -1);

	/* remove node if it reached the end */
	// if(rn->path_ridx == 0) {
	if(qpos.pos > rn->qlim) {
		debug("remove rnode, rn(%p), pridx(%lld), next(%p)", rn, rn->path_ridx, next);
		rbtree_remove(ctx->rtree, (rbtree_node_t *)rn);
	}
	return(next);
}

/**
 * @fn rtree_replace
 */
static _force_inline
void rtree_replace(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct qtree_node_s *qhead,
	// struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	/* replace path array */
	uint64_t const *path = (uint64_t const *)aln->path->array;

	/* replace root section */
	struct gaba_path_section_s const *rsec = &aln->sec[aln->rsidx];
	uint32_t plen = gaba_plen(rsec);

	debug("replace rn(%p), prev(%p, %llu), new(%p, %llu), diff(%lld)",
		rn, rn->aln, rn->h.key,
		aln, _cast_u(((struct gref_gid_pos_s){
			.gid = rsec->aid,
			.pos = aln->rapos + ctx->conf.overlap_width
		})),
		rn->h.key - _cast_u(((struct gref_gid_pos_s){
			.gid = rsec->aid,
			.pos = aln->rapos + ctx->conf.overlap_width
		})));

	rn->h.key = _cast_u(((struct gref_gid_pos_s){
		.gid = rsec->aid,
		.pos = aln->rapos + ctx->conf.overlap_width
	}));
	rn->path_ridx = (plen - aln->rppos) & ~(64 - 1),
	rn->ptail = path + ((rsec->ppos + plen)>>6);

	rn->aln = aln;
	rn->sidx = aln->rsidx;

	rn->qhead = qhead;
	return;
}

/**
 * @fn rtree_adjust_head
 */
static _force_inline
void rtree_adjust_head(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln)
{
	/* head of the root section is aligned */
	rn->h.key = _cast_u(((struct gref_gid_pos_s){
		.gid = _cast_p(rn->h.key).gid,
		.pos = aln->rapos + ctx->conf.overlap_width
	}));
	rn->path_ridx = gaba_plen(&rn->aln->sec[rn->sidx]) - aln->rppos;
	return;
}

/**
 * @fn rtree_adjust_tail
 */
static _force_inline
void rtree_adjust_tail(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln)
{
	/* tail of the root section is aligned */
	rn->h.key = _cast_u(((struct gref_gid_pos_s){
		.gid = _cast_p(rn->h.key).gid,
		.pos = aln->rapos + ctx->conf.overlap_width
	}));
	rn->path_ridx = gaba_plen(&aln->sec[aln->rsidx]) - aln->rppos;
	return;
}

/**
 * @fn qtree_advance
 * @brief search qtree with qpos, create rnode from matched qnode (if found).
 */
static _force_inline
struct qtree_node_s *qtree_advance(
	struct ggsea_ctx_s *ctx,
	struct qtree_node_s *qn,
	struct gref_gid_pos_s qpos)
{
	while(qn != NULL && (uint64_t)qn->h.key == _cast_u(qpos)) {
		/*
		 * current qpos hit at least one node in qtree
		 */
		struct rtree_node_s *rn = (struct rtree_node_s *)
			rbtree_create_node(ctx->rtree);

		/* build rnode from qnode */
		struct gaba_path_section_s const *sec = &qn->aln->sec[qn->sidx];
		uint32_t plen = gaba_plen(sec);

		*rn = (struct rtree_node_s){
			.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
				.gid = sec->aid,
				.pos = sec->apos + ctx->conf.overlap_width
			})),
			.prev_qpos = qpos.pos,
			.path_qpos = qpos.pos,
			.qlim = sec->bpos + sec->blen,

			.path_ridx = plen & ~(64 - 1),
			.ptail = (uint64_t const *)qn->aln->path + ((sec->ppos + plen)>>6),

			.aln = qn->aln,
			.sidx = qn->sidx
		};
		rbtree_insert(ctx->rtree, (rbtree_node_t *)rn);

		/* fetch the next node */
		qn = (struct qtree_node_s *)rbtree_right(ctx->qtree, (rbtree_node_t *)qn);
	}
	return(qn);
}

/**
 * @fn qtree_refresh_node
 */
static _force_inline
struct qtree_node_s *qtree_refresh_node(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s qpos)
{
	return((struct qtree_node_s *)rbtree_search_key_right(ctx->qtree, _cast_u(qpos)));
}

/**
 * @fn qtree_append_result
 * @brief append alignment to qtree, update qtree pointer for the next r-iteration.
 */
static _force_inline
struct qtree_node_s *qtree_append_result(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln,
	uint32_t res_id)
{
	struct qtree_node_s *prev = NULL, *head = NULL;

	/* append all sections to qtree except root */
	for(int64_t i = 0; i < aln->slen; i++) {

		/* create qnode */
		struct qtree_node_s *qn = (struct qtree_node_s *)
			rbtree_create_node(ctx->qtree);

		/* set index and result pointer */
		struct gaba_path_section_s const *sec = &aln->sec[i];
		*qn = (struct qtree_node_s){
			.h.key = _cast_u(((struct gref_gid_pos_s){ sec->aid, sec->apos })),
			.sidx = i,
			.aln = aln,
			.res_id = res_id
		};

		/* set link pointer */
		if(prev == NULL) {
			prev = head = qn;
		} else {
			prev->next = qn; prev = qn;
		}

		/* append it to tree */
		rbtree_insert(ctx->qtree, (rbtree_node_t *)qn);
	}

	/* set tail pointer */
	if(prev != NULL) { prev->next = NULL; }

	/* update qtree pointer for the next iteration */
	return(head);
}

/**
 * @fn qtree_replace
 */
static _force_inline
struct qtree_node_s *qtree_replace(
	struct ggsea_ctx_s *ctx,
	struct qtree_node_s *_qn,
	// struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln,
	uint32_t res_id,
	int64_t ofs)
{
	#define _set_qn(qn, i) { \
		struct gaba_path_section_s const *sec = &aln->sec[(i)]; \
		(qn)->h.key = _cast_u(((struct gref_gid_pos_s){ sec->aid, sec->apos })); \
		(qn)->sidx = (i); \
		(qn)->aln = aln; \
		(qn)->res_id = res_id; \
		if(prev == NULL) { prev = head = (qn); } else { prev->next = (qn); prev = (qn); } \
	}

	struct qtree_node_s *prev = NULL, *head = NULL;
	int64_t i;

	/* create qnodes for head sections */
	for(i = 0; i < ofs; i++) {
		struct qtree_node_s *qn = (struct qtree_node_s *)rbtree_create_node(ctx->qtree);
		_set_qn(qn, i);
		rbtree_insert(ctx->qtree, (rbtree_node_t *)qn);
	}

	/* replace aligned sections */
	for(struct qtree_node_s *qn = _qn;
		qn != NULL && i < aln->slen;
		qn = qn->next, i++) {
		_set_qn(qn, i);
	}

	/* create qnodes for tail sections */
	for(; i < aln->slen; i++) {
		struct qtree_node_s *qn = (struct qtree_node_s *)rbtree_create_node(ctx->qtree);
		_set_qn(qn, i);
		rbtree_insert(ctx->qtree, (rbtree_node_t *)qn);
	}

	/* set tail pointer */
	if(prev != NULL) { prev->next = NULL; }

	#undef _set_qn
	return(head);
}


/* adjacent filter */
/**
 * @fn adjacent_filter_skip_nodes
 */
static _force_inline
struct gref_gid_pos_s const *adjacent_filter_skip_nodes(
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s const *parr,
	struct gref_gid_pos_s const *ptail)
{
	while(parr < ptail && _cast_u(*parr) + 1 < _cast_u(rpos)) {
		parr++;
	}
	return(parr);
}

/**
 * @fn adjacent_filter_test
 */
static _force_inline
int64_t adjacent_filter_test(
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s const *parr,
	struct gref_gid_pos_s const *ptail)
{
	debug("adjacent filter test, ppos(%u), rpos(%u)", parr->pos, rpos.pos);
	return(parr < ptail && _cast_u(*parr) + 1 == _cast_u(rpos));
}


/* overlap filters */
/**
 * @struct rtree_node_pair_s
 */
struct rtree_node_pair_s {
	struct rtree_node_s *left;
	struct rtree_node_s *right;
};

/**
 * @fn overlap_filter_skip_nodes
 */
static _force_inline
struct rtree_node_pair_s overlap_filter_skip_nodes(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_pair_s r,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	while(r.right != NULL && rtree_update(ctx, r.right, qpos) < _cast_u(rpos)) {
		r.right = rtree_advance(ctx, (r.left = r.right), qpos);
	}
	return(r);
}

/**
 * @fn overlap_filter_test
 */
static _force_inline
int64_t overlap_filter_test(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_pair_s r,
	struct gref_gid_pos_s rpos)
{
	/* check if seed overlaps with previous results */
	debug("overlap filter test, rn(%p), rpos(%llu), diff(%lld)",
		r.right,
		(r.right != NULL) ? (0xffffffff & r.right->h.key) : (int64_t)-1,
		(r.right != NULL) ? (r.right->h.key - _cast_u(rpos)) : (int64_t)-1);

	/*
	if(r.right != NULL && (uint64_t)(r.right->h.key + 1000 - _cast_u(rpos)) < 2000) {
		debug("%lld", r.right->h.key - _cast_u(rpos));
	}
	*/

	uint64_t window = 2 * ctx->conf.overlap_width;
	return(r.right != NULL && (uint64_t)(r.right->h.key - _cast_u(rpos)) < window);
}


/* result vector handling */
/**
 * @struct resv_score_pos_s
 */
struct resv_score_pos_s {
	uint32_t idx;
	uint32_t pos;
	int64_t score;
};

/**
 * @fn resv_register
 */
static _force_inline
uint32_t resv_register(
	struct ggsea_ctx_s *ctx,
	struct gaba_alignment_s const *aln)
{
	lmm_kv_push(ctx->res_lmm, ctx->aln, aln);
	return(lmm_kv_size(ctx->aln) - 1);
}

/**
 * @fn resv_replace
 */
static _force_inline
void resv_replace(
	struct ggsea_ctx_s *ctx,
	uint32_t idx,
	struct gaba_alignment_s const *aln)
{
	lmm_kv_at(ctx->aln, idx) = aln;
	return;
}

/**
 * @fn resv_unregister
 */
static _force_inline
void resv_unregister(
	struct ggsea_ctx_s *ctx,
	uint32_t idx)
{
	lmm_kv_at(ctx->aln, idx) = NULL;
	return;
}

/**
 * @fn resv_cmp_result
 */
static _force_inline
int64_t resv_cmp_result(
	struct gaba_alignment_s const *const *aln,
	struct resv_score_pos_s const *karr,
	int64_t i,
	int64_t j)
{
	debug("compare i(%lld) and j(%lld), [i].score(%lld), [j].score(%lld), [i].aid(%u), [i].bid(%u), [j].aid(%u), [j].bid(%u)",
		i, j,
		aln[karr[i].idx]->score, aln[karr[i].idx]->score,
		aln[karr[i].idx]->sec[0].aid, aln[karr[i].idx]->sec[0].bid,
		aln[karr[j].idx]->sec[0].aid, aln[karr[j].idx]->sec[0].bid);
	return(aln[karr[i].idx]->score - aln[karr[j].idx]->score);
}

/**
 * @fn resv_dedup_result
 */
static _force_inline
int64_t resv_dedup_result(
	struct ggsea_ctx_s *ctx,
	struct gaba_alignment_s const **aln,
	int64_t const cnt)
{
	int64_t eff_cnt = 0;
	struct resv_score_pos_s karr[cnt];
	for(int64_t i = 0; i < cnt; i++) {
		if(aln[i] == NULL) { continue; }

		debug("score(%lld)", aln[i]->score);
		karr[eff_cnt++] = (struct resv_score_pos_s){
			.score = -aln[i]->score,		/* score in descending order */
			.pos = aln[i]->sec[0].aid + aln[i]->sec[0].bid,
			.idx = i
		};
		debug("pushed, i(%u), score(%lld), pos(%u)",
			karr[eff_cnt - 1].idx, karr[eff_cnt - 1].score, karr[eff_cnt - 1].pos);
	}
	psort_full(karr, eff_cnt, 16, 0);

	for(int64_t i = 0; i < eff_cnt; i++) {
		debug("sorted, i(%lld), idx(%u), score(%lld), pos(%u), score(%lld)",
			i, karr[i].idx, karr[i].score, karr[i].pos, aln[karr[i].idx]->score);
	}

	/* dedup */
	int64_t j = 0;
	for(int64_t i = 0; i < eff_cnt; i++) {
		if(resv_cmp_result(aln, karr, i, j) == 0) {
			continue;
		}

		debug("move to next, j(%lld)", j + 1);
		karr[++j] = karr[i];
	}
	int64_t dedup_cnt = j + 1;

	/* build shrinked result array */
	#if 0
	for(int64_t i = 0; i < dedup_cnt; i++) {
		/* swap elem */
		struct gaba_alignment_s const *tmp = aln[i];
		aln[i] = aln[karr[i].idx];
		aln[karr[i].idx] = tmp;
		debug("i(%lld) push(%u)", i, karr[i].idx);
	}
	debug("dedup finished, ptr(%p), cnt(%lld), dedup_cnt(%lld)", aln, cnt, dedup_cnt);
	#endif
	return(dedup_cnt);
}

/**
 * @fn resv_pack_result
 */
static _force_inline
struct ggsea_result_s *resv_pack_result(
	struct ggsea_ctx_s *ctx)
{
	/* build array and sort */
	struct gaba_alignment_s const **aln = lmm_kv_ptr(ctx->aln);
	int64_t const cnt = lmm_kv_size(ctx->aln);

	/* dedup result */
	int64_t dedup_cnt = (cnt != 0) ? resv_dedup_result(ctx, aln, cnt) : 0;

	/* pack pointer and length */
	struct ggsea_result_s *res = lmm_malloc(ctx->res_lmm, sizeof(struct ggsea_result_s));
	*res = (struct ggsea_result_s){
		.reserved1 = (void *)ctx->res_lmm,
		.ref = ctx->r,
		.query = ctx->q,
		.aln = aln,
		.cnt = dedup_cnt,
		.reserved2 = cnt
	};
	debug("result, ptr(%p), aln(%p)", res, res->aln);
	return(res);
}


/* postprocessing */
/**
 * @struct pp_match_s
 */
struct pp_match_s {
	int64_t cmp;
	uint32_t xidx, yidx;
};

/**
 * @struct pp_recomb_s
 */
struct pp_recomb_s {
	uint32_t hidx, tidx;
};

/**
 * @fn pp_clip_cmp
 */
static _force_inline
int64_t pp_clip_cmp(
	int64_t cmp)
{
	return((cmp > 0) ? 1 : ((cmp == 0) ? 0 : -1));
}

/**
 * @fn pp_match_section_forward
 */
static _force_inline
struct pp_match_s pp_match_section_forward(
	struct gaba_alignment_s const *x,
	uint32_t xsid,
	struct gaba_alignment_s const *y,
	uint32_t ysid)
{
	struct gaba_path_section_s const *xp = &x->sec[xsid];
	struct gaba_path_section_s const *yp = &y->sec[ysid];

	int64_t yofs = (int64_t)(x->slen - xsid) - (int64_t)(y->slen - ysid);
	struct gaba_path_section_s const *xlim = &x->sec[x->slen - MIN2(yofs, 0)];

	/* iterate until xp reaches either tail */
	while(xp < xlim) {

		debug("xp(%p), yp(%p), apos(%u, %u), alen(%u, %u), bpos(%u, %u), blen(%u, %u), atail(%u, %u), btail(%u, %u), diff(%lld, %lld)",
			xp, yp,
			xp->apos, yp->apos, xp->alen, yp->alen,
			xp->bpos, yp->bpos, xp->blen, yp->blen,
			xp->apos + xp->alen,
			yp->apos + yp->alen,
			xp->bpos + xp->blen,
			yp->bpos + yp->blen,
			(int64_t)xp->apos + xp->alen - yp->apos - yp->alen,
			(int64_t)xp->bpos + xp->blen - yp->bpos - yp->blen);

		/* test tail of the sections are aligned */
		if(((xp->aid - yp->aid) | (xp->bid - yp->bid)
		| ((xp->apos + xp->alen) - (yp->apos + yp->alen))
		| ((xp->bpos + xp->blen) - (yp->bpos + yp->blen))) != 0) {

			int64_t pofs = (int64_t)gaba_plen(xp) - (int64_t)gaba_plen(yp);
			debug("yofs(%lld), xplen(%u), yplen(%u), det(%lld)",
				yofs, gaba_plen(xp), gaba_plen(yp), (yofs<<1) + (pofs>>63));

			return((struct pp_match_s){
				.cmp = pp_clip_cmp((yofs == 0) ? pofs : yofs),
				.xidx = xp - x->sec,
				.yidx = yp - y->sec
			});
		}

		/* aligned, go next */
		xp++; yp++;
	}

	/* all the sections are aligned */
	debug("yofs(%lld)", yofs);
	return((struct pp_match_s){
		.cmp = pp_clip_cmp(yofs),
		.xidx = xp - x->sec,
		.yidx = yp - y->sec
	});
}

/**
 * @fn pp_match_section_reverse
 */
static _force_inline
struct pp_match_s pp_match_section_reverse(
	struct gaba_alignment_s const *x,
	uint32_t xsid,
	struct gaba_alignment_s const *y,
	uint32_t ysid)
{
	struct gaba_path_section_s const *xp = &x->sec[xsid];
	struct gaba_path_section_s const *yp = &y->sec[ysid];

	int64_t yofs = xsid - ysid;
	struct gaba_path_section_s const *xlim = &x->sec[MIN2(yofs, 0)];

	/* iterate until xp reaches either head */
	while(xp >= xlim) {

		debug("xp(%p), yp(%p), ahead(%u, %u), bhead(%u, %u), diff(%lld, %lld)",
			xp, yp,
			xp->apos,
			yp->apos,
			xp->bpos,
			yp->bpos,
			(int64_t)xp->apos - yp->apos,
			(int64_t)xp->bpos - yp->bpos);

		/* test head of the sections are aligned */
		if(((xp->aid - yp->aid) | (xp->bid - yp->bid)
		| (xp->apos - yp->apos) | (xp->bpos - yp->bpos)) != 0) {

			int64_t pofs = (int64_t)gaba_plen(xp) - (int64_t)gaba_plen(yp);
			debug("yofs(%lld), xplen(%u), yplen(%u), det(%lld)",
				yofs, gaba_plen(xp), gaba_plen(yp), (yofs<<1) + (pofs>>63));

			return((struct pp_match_s){
				.cmp = pp_clip_cmp((yofs == 0) ? pofs : yofs),
				.xidx = xp - x->sec + 1,
				.yidx = yp - y->sec + 1
			});
		}

		/* aligned, go next */
		xp--; yp--;
	}

	/* all the sections are aligned */
	debug("yofs(%lld)", yofs);
	return((struct pp_match_s){
		.cmp = pp_clip_cmp(yofs),
		.xidx = xp - x->sec + 1,
		.yidx = yp - y->sec + 1
	});
}

/**
 * @fn pp_is_head_aligned
 */
static _force_inline
int64_t pp_is_head_aligned(
	struct gaba_alignment_s const *x,
	struct pp_match_s h)
{
	return(h.xidx <= x->rsidx);
}

/**
 * @fn pp_is_tail_aligned
 */
static _force_inline
int64_t pp_is_tail_aligned(
	struct gaba_alignment_s const *x,
	struct pp_match_s t)
{
	return(x->rsidx < t.xidx);
}

/**
 * @fn pp_id
 */
static
struct gaba_alignment_s const *pp_id(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	/* do nothing */
	return(aln);
}

/**
 * @fn pp_adjust_head
 */
static
struct gaba_alignment_s const *pp_adjust_head(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	rtree_adjust_head(ctx, rn, aln);
	gaba_dp_res_free((struct gaba_alignment_s *)aln);
	return(NULL);
}

/**
 * @fn pp_adjust_tail
 */
static
struct gaba_alignment_s const *pp_adjust_tail(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	rtree_adjust_tail(ctx, rn, aln);
	gaba_dp_res_free((struct gaba_alignment_s *)aln);
	return(NULL);
}

/**
 * @fn pp_calc_recomb_pos
 */
static _force_inline
struct pp_recomb_s pp_calc_recomb_pos(
	struct gaba_alignment_s const *h,
	uint32_t hidx,
	struct gaba_alignment_s const *t,
	uint32_t tidx,
	uint32_t len)
{
	int64_t pacc = 0, pmin = 0;
	uint32_t hmin = hidx, tmin = tidx;

	struct gaba_path_section_s const *hp = &h->sec[hidx];
	struct gaba_path_section_s const *tp = &t->sec[tidx];

	for(uint32_t i = 0; i < len; i++) {
		pacc += gaba_plen(hp); hp++;
		pacc -= gaba_plen(tp); tp++;
		if(pacc < pmin) {
			pmin = pacc;
			hmin = hidx + i + 1;
			tmin = tidx + i + 1;
		}
	}
	return((struct pp_recomb_s){
		.hidx = hmin,
		.tidx = tmin
	});
}

/**
 * @fn pp_recomb_head
 */
static
struct gaba_alignment_s const *pp_recomb_head(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	struct gaba_alignment_s const *x = aln, *y = rn->aln;
	uint32_t xsidx = aln->rsidx, ysidx = rn->sidx;

	struct pp_recomb_s p = pp_calc_recomb_pos(x, xidx, y, yidx, len);
	x = gaba_dp_recombine(ctx->dp,
		(struct gaba_alignment_s *)x, p.hidx,
		(struct gaba_alignment_s *)y, p.tidx);

	/* replace nodes */
	struct qtree_node_s *qhead = rn->qhead;
	resv_replace(ctx, qhead->res_id, x);
	qhead = qtree_replace(ctx, qhead, x, qhead->res_id, (int64_t)xsidx - (int64_t)ysidx);
	rtree_replace(ctx, rn, qhead, x);
	return(NULL);
}

/**
 * @fn pp_recomb_tail
 */
static
struct gaba_alignment_s const *pp_recomb_tail(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	struct gaba_alignment_s const *x = aln, *y = rn->aln;
	uint32_t xsidx = aln->rsidx, ysidx = rn->sidx;

	struct pp_recomb_s p = pp_calc_recomb_pos(y, yidx, x, xidx, len);
	x = gaba_dp_recombine(ctx->dp,
		(struct gaba_alignment_s *)y, p.hidx,
		(struct gaba_alignment_s *)x, p.tidx);

	/* replace nodes */
	struct qtree_node_s *qhead = rn->qhead;
	resv_replace(ctx, qhead->res_id, x);
	qhead = qtree_replace(ctx, qhead, x, qhead->res_id, (int64_t)ysidx - (int64_t)xsidx);
	rtree_replace(ctx, rn, qhead, x);
	return(NULL);
}

/**
 * @fn pp_replace
 */
static
struct gaba_alignment_s const *pp_replace(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gaba_alignment_s const *aln,
	uint32_t xidx,
	uint32_t yidx,
	uint32_t len)
{
	uint32_t xsidx = aln->rsidx, ysidx = rn->sidx;

	/* remove previous result */
	struct qtree_node_s *qhead = rn->qhead;
	gaba_dp_res_free((struct gaba_alignment_s *)rn->aln);

	/* replace nodes */
	resv_replace(ctx, qhead->res_id, aln);
	qhead = qtree_replace(ctx, qhead, aln, qhead->res_id, (int64_t)xsidx - (int64_t)ysidx);
	rtree_replace(ctx, rn, qhead, aln);
	return(NULL);
}

/**
 * @fn pp_process_alignment
 */
static _force_inline
struct gaba_alignment_s const *pp_process_alignment(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *rn,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	struct gaba_alignment_s const *x = aln, *y = rn->aln;
	uint32_t xsidx = aln->rsidx, ysidx = rn->sidx;

	struct pp_match_s h = pp_match_section_reverse(x, xsidx, y, ysidx);
	struct pp_match_s t = pp_match_section_forward(x, xsidx, y, ysidx);

	if(!pp_is_head_aligned(x, h) && !pp_is_tail_aligned(x, t)) {
		/* both ends are not aligned; ignore */
		return(x);
	}

	struct gaba_alignment_s const *(*process[3][3])(
		struct ggsea_ctx_s *ctx,
		struct rtree_node_s *rn,
		struct gaba_alignment_s const *aln,
		uint32_t xidx,
		uint32_t yidx,
		uint32_t len) = {
		[0] = {
			[0] = pp_id,
			[1] = pp_adjust_tail,
			[2] = pp_recomb_head
		},
		[1] = {
			[0] = pp_adjust_head,
			[1] = pp_adjust_tail,
			[2] = pp_replace
		},
		[2] = {
			[0] = pp_recomb_tail,
			[1] = pp_replace,
			[2] = pp_replace
		},
	};
	return(process[h.cmp + 1][t.cmp + 1](ctx, rn, aln, h.xidx, h.yidx, t.xidx - h.xidx));
	// return(aln);
}

/**
 * @fn ggsea_append_result
 */
static _force_inline
struct rtree_node_s *ggsea_append_result(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	uint32_t res_id = resv_register(ctx, aln);
	struct qtree_node_s *qhead = qtree_append_result(ctx, qpos, aln, res_id);
	return(rtree_append_result(ctx, qhead, qpos, aln));
}

/**
 * @fn ggsea_evaluate_alignment
 */
static _force_inline
struct rtree_node_pair_s ggsea_evaluate_alignment(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_pair_s r,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	debug("right(%p), left(%p)", r.right, r.left);

	if(r.right != NULL) {
		if((aln = pp_process_alignment(ctx, r.right, rpos, qpos, aln)) == NULL) {
			debug("replaced with right");
			return(r);
		}		
	}

	if(r.left != NULL) {
		if((aln = pp_process_alignment(ctx, r.left, rpos, qpos, aln)) == NULL) {
			debug("replaced with left");
			return(r);
		}		
	}

	return((struct rtree_node_pair_s){
		.left = r.left,
		.right = ggsea_append_result(ctx, qpos, aln)
	});

	#if 0
	/* current rnode */
	if(rn == NULL) {
		return(ggsea_append_result(ctx, rn, qpos, aln));
	}
	if((aln = pp_process_alignment(ctx, rn, rpos, qpos, aln)) == NULL) {
		return(rn);
	}

	/* previous rnode */
	rn = (struct rtree_node_s *)rbtree_left(ctx->rtree, (rbtree_node_t *)rn);
	if(rn == NULL) {
		return(ggsea_append_result(ctx, rn, qpos, aln));
	}
	if((aln = pp_process_alignment(ctx, rn, rpos, qpos, aln)) == NULL) {
		return(rn);
	}

	return(ggsea_append_result(ctx, rn, qpos, aln));
	#endif
}


/**
 * @fn ggsea_evaluate_seeds
 */
static _force_inline
struct qtree_node_s *ggsea_evaluate_seeds(
	struct ggsea_ctx_s *ctx,
	struct qtree_node_s *qn,
	uint64_t kmer,
	struct gref_gid_pos_s const *rarr,		/* current seeds on reference side */
	int64_t rlen,
	struct gref_gid_pos_s const *parr,		/* previous seeds (for adjacent filter) */
	int64_t plen,
	struct gref_gid_pos_s qpos)
{
	debug("kmer(%llx), gid(%u), pos(%u), m.len(%lld)",
		kmer, qpos.gid, qpos.pos, rlen);

	/* init ptail for adjacent filter */
	struct gref_gid_pos_s const *ptail = parr + plen;
	debug("init adjacent filter, parr(%p), ptail(%p), plen(%lld)", parr, ptail, plen);

	/* iterate over rtree */
	struct rtree_node_pair_s r = {
		.left = NULL,
		.right = (struct rtree_node_s *)rbtree_search_key_right(
			ctx->rtree, INT64_MIN)
	};
	debug("init rnode, rn(%p, %lld)", r.right, (r.right != NULL) ? r.right->h.key : -1);
	for(int64_t i = 0; i < rlen; i++) {
		struct gref_gid_pos_s rpos = rarr[i];

		debug("seed: i(%lld), rpos(%llx)", i, _cast_u(rpos));

		/* first check the adjacency to the previous seeds */
		parr = adjacent_filter_skip_nodes(rpos, parr, ptail);
		if(adjacent_filter_test(rpos, parr, ptail)) { continue; }

		/* next overlap filter */
		r = overlap_filter_skip_nodes(ctx, r, rpos, qpos);
		if(overlap_filter_test(ctx, r, rpos)) { continue; }
		debug("filter passed, i(%lld), rpos(%u)", i, rpos.pos);

		/* extend */
		struct gaba_alignment_s const *aln = dp_extend_seed(ctx, rpos, qpos);
		if(aln == NULL) { continue; }

		/* postprocess */
		r = ggsea_evaluate_alignment(ctx, r, rpos, qpos, aln);
		qn = qtree_refresh_node(ctx, qpos);
		debug("extend finished, rn(%p), qn(%p), score(%lld)", r.right, qn, (aln != NULL) ? aln->score : 0);
	}
	return(qn);
}

/**
 * @fn ggsea_align
 */
ggsea_result_t *ggsea_align(
	ggsea_ctx_t *_ctx,
	gref_acv_t const *query,
	gref_iter_t *iter,
	void *_lmm)
{
	struct ggsea_ctx_s *ctx = (struct ggsea_ctx_s *)_ctx;
	lmm_t *lmm = (lmm_t *)_lmm;
	debug("align entry, check iter(%p)", iter);

	/* flush buffer */
	ggsea_flush(ctx, query, lmm);

	/* no results are stored in qtree on init */
	struct qtree_node_s *qn = NULL;

	/* initialize adjacent filter */
	struct gref_gid_pos_s dec = { .gid = -1, .pos = 0 };
	struct gref_match_res_s const init = { .gid_pos_arr = &dec, .len = 0 };
	struct gref_match_res_s p = init;

	/* iterate seeds over query sequence */
	struct gref_kmer_tuple_s t;
	while((t = gref_iter_next(iter)).gid_pos.gid != (uint32_t)-1) {

		/* fetch the next intersecting region */
		qn = qtree_advance(ctx, qn, t.gid_pos);

		/* match and evaluate */
		struct gref_match_res_s m = gref_match_2bitpacked(ctx->r, t.kmer);

		debug("fetched kmer iterator ptr(%p), len(%lld)", m.gid_pos_arr, m.len);

		/* skip if no seeds found */
		if(m.len == 0) {
			p = init; continue;
		}

		/* skip if too many seeds found (mark repetitive) */
		if(m.len > ctx->conf.params.kmer_cnt_thresh) {
			rep_save_pos(ctx, t.kmer, m.gid_pos_arr[0], t.gid_pos);
			p = init; continue;
		}

		/* evaluate */
		qn = ggsea_evaluate_seeds(ctx, qn, t.kmer,
			m.gid_pos_arr, m.len,
			p.gid_pos_arr, p.len,
			t.gid_pos);

		/* save previous seeds */
		p = m;
	}

	/* cleanup iterator */
	debug("done. %llu alignments generated", lmm_kv_size(ctx->aln));
	return((ggsea_result_t *)resv_pack_result(ctx));
}

/**
 * @fn ggsea_aln_free
 */
void ggsea_aln_free(
	ggsea_result_t *res)
{
	if(res == NULL) { return; }

	lmm_t *lmm = (lmm_t *)res->reserved1;
	int64_t cnt = res->reserved2;

	for(int64_t i = 0; i < cnt; i++) {
		gaba_dp_res_free((struct gaba_alignment_s *)res->aln[i]);
	}
	debug("free, ptr(%p), aln(%p)", res, res->aln);
	lmm_free(lmm, (void *)res->aln);
	lmm_free(lmm, (void *)res);
	return;
}


/* unittests */
/* global configuration */
unittest_config(
	.name = "ggsea",
	.depends_on = { "hmap", "psort", "tree", "gref", "gaba" }
);

#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

#define _pool(_k)		( GREF_PARAMS( .k = (_k), .seq_head_margin = 32, .seq_tail_margin = 32 ) )

/* create conf */
unittest()
{
	ggsea_conf_t *conf = ggsea_conf_init(GGSEA_PARAMS(
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1),
		.xdrop = 10));

	assert(conf != NULL, "%p", conf);

	ggsea_conf_clean(conf);
}

/* create context */
unittest()
{
	ggsea_conf_t *conf = ggsea_conf_init(GGSEA_PARAMS(
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1),
		.xdrop = 10));

	gref_pool_t *pool = gref_init_pool(_pool(4));
	gref_append_segment(pool, _str("seq1"), _seq("ACGTACGTACGTAACCACGTACGTACGT"));
	gref_acv_t *acv = gref_freeze_pool(pool);
	gref_idx_t *idx = gref_build_index(acv);
	assert(idx != NULL, "%p", idx);

	/* ctx_init fails without reference index */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, NULL);
	assert(sea == NULL, "%p", sea);

	/* with valid reference index */
	sea = ggsea_ctx_init(conf, idx);
	assert(sea != NULL, "%p", sea);

	/* cleanup */
	ggsea_ctx_clean(sea);
	ggsea_conf_clean(conf);
}

/* omajinais */
#define with_default_conf(_k) \
	.init = (void *(*)(void *))ggsea_conf_init, \
	.clean = (void (*)(void *))ggsea_conf_clean, \
	.params = (void *)GGSEA_PARAMS( \
		.score_matrix = GABA_SCORE_SIMPLE(2, 3, 5, 1), \
		.xdrop = 10, \
		.k = (_k))

#define omajinai() \
	ggsea_conf_t *conf = (ggsea_conf_t *)ctx;

/* omajinai test */
unittest(with_default_conf(4))
{
	omajinai();
	assert(conf != NULL, "%p", conf);
}

/* single linear sequence */
unittest(with_default_conf(4))
{
	omajinai();

	/* build sequence pool */
	gref_pool_t *pool = gref_init_pool(_pool(4));
	gref_append_segment(pool, _str("seq1"), _seq("ACGTACGTACGTAACCACGTACGTACGT"));
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));

	/* build ggsea context */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, idx);

	/* align */
	gref_iter_t *iter = gref_iter_init(idx, NULL);
	ggsea_result_t *r = ggsea_align(sea, (gref_acv_t *)idx, iter, NULL);
	assert(r->aln != NULL, "%p", r->aln);

	ggsea_aln_free(r);
	ggsea_ctx_clean(sea);
	gref_iter_clean(iter);
	gref_clean(idx);
}

/* longer sequence */
unittest(with_default_conf(14))
{
	omajinai();

	/* build reference index object */
	gref_pool_t *rpool = gref_init_pool(_pool(14));
	gref_append_segment(rpool, _str("ref1"),
		_seq("CTCACCTCGCTCAAAAGGGCTGCCTCCGAGCGTGTGGGCGAGGACAACCGCCCCACAGTCAAGCTCGAATGGGTGCTATTGCGTAGCTAGGACCGGCACT"));
	gref_idx_t *ref = gref_build_index(gref_freeze_pool(rpool));

	/* build query sequence iterator */
	gref_pool_t *qpool = gref_init_pool(_pool(14));
	gref_append_segment(qpool, _str("query1"),
		_seq("GGCTGCCTCCGAGCGTGTGGGCGAGGACAACCGCCCCACAGTCAAGCTCGAA"));
	gref_acv_t *query = gref_freeze_pool(qpool);


	/* build ggsea context */
	ggsea_ctx_t *sea = ggsea_ctx_init(conf, ref);

	/* align */
	gref_iter_t *iter = gref_iter_init(query, NULL);
	ggsea_result_t *r = ggsea_align(sea, query, iter, NULL);
	assert(r->aln != NULL, "%p", r->aln);

	/* cleanup */
	ggsea_aln_free(r);
	ggsea_ctx_clean(sea);

	gref_clean(ref);
	gref_iter_clean(iter);
	gref_clean(query);
}

/**
 * end of ggsea.c
 */

#if 0

/**
 * @fn ggsea_check_alignment
 */
static _force_inline
void ggsea_check_alignment(
	struct ggsea_ctx_s *ctx,
	struct gaba_alignment_s const *aln)
{
	int64_t acnt = 0, bcnt = 0;
	for(int64_t i = 0; i < aln->path->len; i++) {
		if(((aln->path->array[i / 32]>>(i & 31)) & 0x01) == 0) {
			acnt++;
		} else {
			bcnt++;
		}
	}

	debug("cnt(%lld, %lld), sec[0].len(%u, %u), diff(%lld, %lld), i(%d, %u, %lld), path(%x, %x, %x, %x, %x)",
		bcnt, acnt, aln->sec[0].blen, aln->sec[0].alen,
		bcnt - aln->sec[0].blen, acnt - aln->sec[0].alen,
		0, aln->rppos, aln->path->len,
		aln->path->array[0], aln->path->array[1],
		aln->path->array[aln->rppos / 32], aln->path->array[aln->rppos / 32 + 1],
		aln->path->array[aln->path->len / 32]);

	// ((struct gaba_path_section_s *)aln->sec)[0].alen = acnt;
	// ((struct gaba_path_section_s *)aln->sec)[0].blen = bcnt;
	return;
}

/**
 * @fn pp_match_section
 */
struct pp_match_s {
	struct gaba_alignment_s const *haln;
	struct gaba_alignment_s const *taln;
	struct gaba_path_section_s const *hsec;
	struct gaba_path_section_s const *tsec;
	struct gaba_path_section_s const *xp;
	struct gaba_path_section_s const *yp;
};
static _force_inline
struct pp_match_s pp_match_section(
	struct ggsea_ctx_s *ctx,
	struct gaba_alignment_s const *x,
	uint32_t xsid,
	struct gaba_alignment_s const *y,
	uint32_t ysid)
{
	struct pp_match_s h = pp_match_section_reverse(x, xsid, y, ysid);
	struct pp_match_s t = pp_match_section_forward(x, xsid, y, ysid);
	struct pp_match_s p = pp_detect_recomb_pos(x, xsid, y, ysid, h, t);
	return((struct pp_match_s){
		.haln = h.aln,
		.taln = t.aln,
		.hsec = h.sec,
		.tsec = t.sec,
		.xp = p.xp,
		.yp = p.yp,
	});
}

#endif

