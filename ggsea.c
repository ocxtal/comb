
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
	int64_t init_rep_hash_size;
	int64_t max_rep_vec_size;		/* max kmer vector size */
	int64_t overlap_width;			/* overlap filter width */
	int64_t res_lmm_size;			/* result memory manager size */
	struct ggsea_params_s params;
};

/**
 * @struct ggsea_rep_seed_s
 */
struct ggsea_rep_seed_s {
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
	rbtree_node_t h;		/* key: (id, pos) pair on reference side */
	uint64_t qpos;			/* (id, pos) pair on query side */
	uint32_t uridx;			/* remaining length until the next update */
	uint32_t pridx;			/* ramaining length until the end of the section */
	int64_t plim;			/* tail index on the path string */
	uint64_t const *path;	/* path string */
	struct gaba_alignment_s const *aln;
};
_static_assert(sizeof(struct rtree_node_s) == 80);
#else
struct rtree_node_s {
	rbtree_node_t h;		/* (40) */
	uint32_t prev_qpos;		/* q-coordinate at the previous r-index update */
	uint32_t path_qpos;		/* q-coordinate at the previous path adjustment */
	uint32_t path_ridx;		/* reverse index of the path string after the previous path adjustment */
	uint32_t path_rem;		/* rem length from the tail */
	uint64_t const *ptail;	/* path tail of the current section */
};
#endif

/**
 * @struct qtree_node_s
 */
struct qtree_node_s {
	rbtree_node_t h;		/* key: (id, pos) pair on query side */	
	struct gaba_alignment_s const *aln;
	uint32_t sidx;			/* section array index */
	uint32_t pad[3];
};
_static_assert(sizeof(struct qtree_node_s) == 64);

/**
 * @struct ggsea_fill_pair_s
 */
struct ggsea_fill_pair_s {
	gaba_fill_t const *fw;
	gaba_fill_t const *rv;
};

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
 * @struct ggsea_front_s
 * @brief segment info container (to push into heapqueue)
 */
struct ggsea_front_s {
	int64_t psum;
	gaba_fill_t const *fill;
	uint32_t rgid;
	uint32_t qgid;
};
_static_assert(sizeof(struct ggsea_front_s) == 24);

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
	kvec_t(struct ggsea_front_s) queue;	/* segment queue */
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
	conf->overlap_width = 32;
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
		int64_t hcnt = hmap_get_count(ctx->rep);
		for(int64_t i = 0; i < hcnt; i++) {
			struct ggsea_rep_seed_s *c = hmap_get_object(ctx->rep, i);
			debug("free rv(%p), qv(%p)", kv_ptr(c->rv), kv_ptr(c->qv));
			kv_destroy(c->rv);
			kv_destroy(c->qv);
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
		sizeof(struct ggsea_rep_seed_s),
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
		struct ggsea_rep_seed_s *c = hmap_get_object(ctx->rep, i);
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
	debug("flushed, aln(%p), lmm(%p), lim(%p)", lmm_kv_ptr(ctx->aln), ctx->res_lmm, ctx->res_lmm->lim);
	return;
}


/* repetitive kmer filters */
/**
 * @fn ggsea_dedup_seed_pos
 */
static _force_inline
int64_t ggsea_dedup_seed_pos(
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
	return(cnt);
}

/**
 * @fn ggsea_save_seed_pos
 */
static _force_inline
void ggsea_save_seed_pos(
	struct ggsea_ctx_s *ctx,
	uint64_t kmer,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* add kmer to the hashmap */
	uint32_t id = hmap_get_id(ctx->rep,
		(char const *)&kmer, sizeof(uint64_t));
	struct ggsea_rep_seed_s *c = hmap_get_object(ctx->rep, id);

	debug("id(%u), count(%u), diff(%d)", id, hmap_get_count(ctx->rep), id - hmap_get_count(ctx->rep) + 1);
	if(id == (hmap_get_count(ctx->rep) - 1)) {
		c->vec_size = ctx->conf.max_rep_vec_size;
		kv_init(c->rv);
		kv_init(c->qv);

		debug("malloc rv(%p), qv(%p)", kv_ptr(c->rv), kv_ptr(c->qv));
	}

	debug("save repetitive kmer(%llx), r(%u, %u), q(%u, %u), rv(%p, %llu), qv(%p, %llu)",
		kmer, rpos.gid, rpos.pos, qpos.gid, qpos.pos,
		kv_ptr(c->rv), kv_size(c->rv), kv_ptr(c->qv), kv_size(c->qv));

	/* add rpos */
	kv_push(c->rv, rpos);
	if(kv_size(c->rv) > c->vec_size) {
		kv_size(c->rv) = ggsea_dedup_seed_pos(ctx, kv_ptr(c->rv), kv_size(c->rv));
	}

	kv_push(c->qv, qpos);
	if(kv_size(c->qv) > c->vec_size) {
		kv_size(c->qv) = ggsea_dedup_seed_pos(ctx, kv_ptr(c->qv), kv_size(c->qv));
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
 * @fn ggsea_extend_leaf
 */
gaba_fill_t const *ggsea_extend_leaf(
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
 * @fn ggsea_extend_update_queue
 */
gaba_fill_t const *ggsea_extend_update_queue(
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
		return(ggsea_extend_leaf(ctx, fill, max, rsec, qsec, trigger_mask));
	}

	/* push section pairs */
	for(int64_t i = 0; i < rlink.len; i++) {
		for(int64_t j = 0; j < qlink.len; j++) {
			debug("push queue, fill(%p), psum(%lld), r(%u), q(%u)",
				fill, fill->psum, rlink.gid_arr[i], qlink.gid_arr[j]);

			kv_hq_push(ctx->queue, ((struct ggsea_front_s){
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
	max = ggsea_extend_update_queue(ctx, fill, max, rsec, qsec);

	/* loop */
	while(kv_hq_size(ctx->queue) > 0) {
		struct ggsea_front_s seg = kv_hq_pop(ctx->queue);
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
			max = ggsea_extend_update_queue(ctx, fill, max, rsec, qsec);
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
struct ggsea_fill_pair_s dp_extend(
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
	return((struct ggsea_fill_pair_s){
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
	sec->plen = 2 * len;

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
			.ppos = 0,
			.plen = 2 * len
		};
	} while((rem -= len) > 0);
	return(sec - sec_base + 1);
}

/**
 * @fn ggsea_extend_seed
 */
static _force_inline
struct gaba_alignment_s const *ggsea_extend_seed(
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
	struct ggsea_fill_pair_s pair = dp_extend(ctx, sec, len);
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
	lmm_kv_push(ctx->res_lmm, ctx->aln, aln);
	return(aln);
}


/* seed filtering */
/**
 * @fn load_u64
 */
static inline
uint64_t load_u64(
	uint64_t const *ptr,
	int64_t pos)
{
	int64_t rem = pos & 63;
	uint64_t a = (ptr[pos>>6]>>rem) | ((ptr[(pos>>6) + 1]<<(63 - rem))<<1);
	return(a);
}

#if 0
/**
 * @fn rtree_append_result
 */
static _force_inline
struct rtree_node_s *rtree_append_result(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *_rn,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	/* create rtree node */
	struct rtree_node_s *rn = (struct rtree_node_s *)
		rbtree_create_node(ctx->rtree);

	/* aln->path->array is always aligned on 8byte boundary */
	uint64_t const *path = (uint64_t const *)aln->path->array;

	/* root section */
	struct gaba_path_section_s const *rsec = &aln->sec[aln->rsid];
	*rn = (struct rtree_node_s){
		.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
			.gid = rsec->aid,
			.pos = aln->rapos + ctx->conf.overlap_width
		})),
		.qpos = _cast_u(qpos),
		.pridx = rsec->plen - aln->rppos,				/* rem. length until the end of the section */
		.plim = rsec->ppos + rsec->plen,				/* tail path index of the section */
		.path = path,
		.aln = aln
	};

	debug("append result, rn(%p), rpos(%lld), qpos(%lld), rsec(%x), h.key(%llx), pridx(%u), plim(%lld), plen(%d), rpos(%d)",
		rn, _cast_u(rpos), _cast_u(qpos), aln->rsid, rn->h.key, rn->pridx, rn->plim, rsec->plen, aln->rppos);
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
	int64_t inc = _cast_u(qpos) - rn->qpos;

	while()


	rn->h.key += (int64_t)(inc * 1.1);			/* increment index on reference side */
	rn->qpos = _cast_u(qpos);
	rn->pridx = MAX2(rn->pridx - 2 * inc, 0);
	debug("update rtree, rpos(%llx), a(%x, %x), b(%x, %x), r(%x, %x), diff(%llx, %llx), inc(%lld)",
		rn->h.key,
		rn->aln->sec[0].apos, rn->aln->sec[0].alen,
		rn->aln->sec[0].bpos, rn->aln->sec[0].blen,
		rn->aln->rapos, rn->aln->rbpos,
		rn->h.key - rn->aln->rapos,
		_cast_u(qpos) - rn->aln->rbpos,
		inc);

	#if 0
	/* update index if needed */
	if(rn->pridx > 0 && (_cast_u(qpos) & (32 - 1)) == 0) {
		int64_t len = MIN2(rn->pridx, 64);
		uint64_t path_array = load_u64(rn->path, rn->plim - rn->pridx);
		int64_t dcnt = popcnt((path_array)<<(64 - len));
		// rn->h.key += (len>>1) - dcnt;
		rn->h.key += len - 2 * dcnt;

		debug("update path, rn(%p), plim(%llx), pridx(%x), len(%lld), dcnt(%lld)", rn, rn->plim, rn->pridx, len, dcnt);
	}
	#endif
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

	if(rn != NULL && next != NULL && rn->h.key > next->h.key) {
		debug("invalid key order rn->key(%llx), next->key(%llx)", rn->h.key, next->h.key);
	}

	/* remove node if it reached the end */
	if(rn->pridx <= 0) {
		debug("remove rnode, rn(%p), pridx(%u)", rn, rn->pridx);
		rbtree_remove(ctx->rtree, (rbtree_node_t *)rn);
	}
	return(next);
}
#else
/**
 * @fn rtree_append_result
 */
static _force_inline
struct rtree_node_s *rtree_append_result(
	struct ggsea_ctx_s *ctx,
	struct rtree_node_s *_rn,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	/* create rtree node */
	struct rtree_node_s *rn = (struct rtree_node_s *)
		rbtree_create_node(ctx->rtree);

	/* aln->path->array is always aligned on 8byte boundary */
	uint64_t const *path = (uint64_t const *)aln->path->array;

	/* root section */
	struct gaba_path_section_s const *rsec = &aln->sec[aln->rsid];
	*rn = (struct rtree_node_s){
		.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
			.gid = rsec->aid,
			.pos = aln->rapos + ctx->conf.overlap_width
		})),
		.prev_qpos = qpos.pos,
		.path_qpos = qpos.pos,
		.path_ridx = rsec->plen - aln->rppos,
		.path_rem = (64 - 1) & rsec->plen,
		.ptail = path + (rsec->ppos + rsec->plen) / 64
	};

	debug("append result, rn(%p), rpos(%u), qpos(%u), rsec(%x), h.key(%llu), prev_qpos(%u), path_ridx(%u), path_rem(%u)",
		rn, rpos.pos, qpos.pos, aln->rsid, 0xffffffff & rn->h.key, rn->prev_qpos, rn->path_ridx, rn->path_rem);
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
		int64_t ridx = rn->path_ridx;
		int64_t rem = rn->path_rem;

		while(q > 0) {

			/* calculate advancing length on q side */
			int64_t qlen = MIN2(q, 32);

			/* load path array */
			uint64_t path_array = load_u64(rn->ptail, rem - ridx);

			/* append phantom diagonals */
			if(ridx < 64) {
				path_array = ((path_array<<(63 - ridx))<<1) | (0xaaaaaaaaaaaaaaaa>>ridx);
			}

			/* count vertical elements */
			int64_t dcnt = popcnt(path_array<<(64 - 2*qlen));

			debug("adjust path, qlen(%lld), dcnt(%lld), rpos(%llu), qrem(%lld), ridx(%lld), path_array(%llx)",
				qlen, dcnt, 0xffffffff & rn->h.key + 2*(qlen - dcnt), q, ridx, path_array);

			rn->h.key += 2 * (qlen - dcnt);
			q -= qlen;
			ridx -= MIN2(ridx, 64);

			#if 0
			int64_t len = MIN3(ridx, 2 * q, 64);
			uint64_t path_array = load_u64(rn->ptail, rem - ridx);
			int64_t dcnt = popcnt(path_array<<(64 - len));

			debug("adjust path, len(%lld), path_array(%llx), dcnt(%lld), rpos(%llu), qrem(%lld), ridx(%lld)",
				len, path_array, dcnt, 0xffffffff & rn->h.key + len - 2 * dcnt, q, ridx);

			rn->h.key += len - 2 * dcnt;
			q -= len>>1;
			ridx -= len;
			#endif
		}

		/* write back indices */
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

	/* remove node if it reached the end */
	if(rn->path_ridx <= 0) {
		debug("remove rnode, rn(%p), pridx(%u), next(%p)",
			rn, rn->path_ridx, next);
		rbtree_remove(ctx->rtree, (rbtree_node_t *)rn);
	}
	return(next);
}
#endif

#if 0
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
	while(qn != NULL && qn->h.key == _cast_u(qpos)) {
		/*
		 * current qpos hit at least one node in qtree
		 */
		struct rtree_node_s *rn = (struct rtree_node_s *)
			rbtree_create_node(ctx->rtree);

		/* build rnode from qnode */
		struct gaba_path_section_s const *sec = &qn->aln->sec[qn->sidx];
		*rn = (struct rtree_node_s){
			.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
				.gid = sec->aid,
				.pos = sec->apos + ctx->conf.overlap_width
			})),
			.qpos = _cast_u(qpos),
			.pridx = sec->blen,
			.plim = sec->ppos + sec->plen,
			.path = (uint64_t const *)qn->aln->path,
			.aln = qn->aln
		};
		rbtree_insert(ctx->rtree, (rbtree_node_t *)rn);

		/* fetch the next node */
		qn = (struct qtree_node_s *)rbtree_right(ctx->qtree, (rbtree_node_t *)qn);
	}
	return(qn);
}
#else
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
	while(qn != NULL && qn->h.key == _cast_u(qpos)) {
		/*
		 * current qpos hit at least one node in qtree
		 */
		struct rtree_node_s *rn = (struct rtree_node_s *)
			rbtree_create_node(ctx->rtree);

		/* build rnode from qnode */
		struct gaba_path_section_s const *sec = &qn->aln->sec[qn->sidx];
		*rn = (struct rtree_node_s){
			.h.key = _cast_u(((struct gref_gid_pos_s){		/* rpos */
				.gid = sec->aid,
				.pos = sec->apos + ctx->conf.overlap_width
			})),
			.prev_qpos = qpos.pos,
			.path_qpos = qpos.pos,
			.path_ridx = sec->plen,
			.path_rem = (64 - 1) & sec->plen,
			.ptail = (uint64_t const *)qn->aln->path + (sec->ppos + sec->plen) / 64
		};
		rbtree_insert(ctx->rtree, (rbtree_node_t *)rn);

		/* fetch the next node */
		qn = (struct qtree_node_s *)rbtree_right(ctx->qtree, (rbtree_node_t *)qn);
	}
	return(qn);
}
#endif

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
	struct qtree_node_s *_qn,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos,
	struct gaba_alignment_s const *aln)
{
	/* append all sections to qtree except root */
	for(int64_t i = 0; i < aln->slen; i++) {
		if(i == aln->rsid) { continue; }			/* skip root section */

		/* create qnode, set index and result pointer, append it to tree */
		struct gaba_path_section_s const *sec = &aln->sec[i];
		struct qtree_node_s *qn = (struct qtree_node_s *)
			rbtree_create_node(ctx->qtree);
		*qn = (struct qtree_node_s){
			.h.key = _cast_u(((struct gref_gid_pos_s){ sec->aid, sec->apos })),
			.sidx = i,
			.aln = aln
		};
		rbtree_insert(ctx->qtree, (rbtree_node_t *)qn);
	}

	/* update qtree pointer for the next iteration */
	return(qtree_refresh_node(ctx, qpos));
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
	struct rtree_node_s *rn = (struct rtree_node_s *)
		rbtree_search_key_right(ctx->rtree, INT64_MIN);
	for(int64_t i = 0; i < rlen; i++) {

		debug("seed: i(%lld), rpos(%llx)", i, _cast_u(rarr[i]));

		/* first check the adjacency to the previous seeds */
		while(parr < ptail && _cast_u(*parr) + 1 < _cast_u(rarr[i])) {
			parr++;
		}

		debug("adjacent filter test, ppos(%u), rpos(%u)", parr->pos, rarr[i].pos);
		if(parr < ptail && _cast_u(*parr) + 1 == _cast_u(rarr[i])) {
			debug("adjacent filter hit, skip seed");
			continue;
		}

		/* get the next rnode */
		while(rn != NULL && rtree_update(ctx, rn, qpos) < _cast_u(rarr[i])) {
			rn = rtree_advance(ctx, rn, qpos);
		}

		/* check if seed overlaps with previous results */
		debug("overlap filter test, rn(%p), i(%lld), rarr[i](%u), rpos(%llu), diff(%lld)",
			rn, i, rarr[i].pos,
			(rn != NULL) ? (0xffffffff & rn->h.key) : (int64_t)-1,
			(rn != NULL) ? (rn->h.key - _cast_u(rarr[i])) : (int64_t)-1);

		if(rn != NULL && (uint64_t)(rn->h.key + 1000 - _cast_u(rarr[i])) < 2000) {
			debug("%lld", rn->h.key - _cast_u(rarr[i]));
		}

		if(rn != NULL && (uint64_t)(rn->h.key - _cast_u(rarr[i])) < 2 * ctx->conf.overlap_width) {
			debug("overlap filter hit, skip seed");
			continue;		/* skip */
		}
		debug("filter passed, i(%lld), rarr[i](%u)", i, rarr[i].pos);

		/* extend */
		struct gaba_alignment_s const *aln = ggsea_extend_seed(ctx, rarr[i], qpos);
		if(aln != NULL) {
			rn = rtree_append_result(ctx, rn, rarr[i], qpos, aln);
			qn = qtree_append_result(ctx, qn, rarr[i], qpos, aln);
		}
		debug("extend finished, rn(%p), qn(%p), score(%lld)", rn, qn, (aln != NULL) ? aln->score : 0);
	}
	return(qn);
}

/**
 * @struct ggsea_score_pos_s
 */
struct ggsea_score_pos_s {
	uint32_t idx;
	uint32_t pos;
	int64_t score;
};

/**
 * @fn ggsea_cmp_result
 */
static _force_inline
int64_t ggsea_cmp_result(
	struct gaba_alignment_s const *const *aln,
	struct ggsea_score_pos_s const *karr,
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
 * @fn ggsea_dedup_result
 */
static _force_inline
int64_t ggsea_dedup_result(
	struct ggsea_ctx_s *ctx,
	struct gaba_alignment_s const **aln,
	int64_t const cnt)
{
	struct ggsea_score_pos_s karr[cnt];
	for(int64_t i = 0; i < cnt; i++) {
		debug("score(%lld)", aln[i]->score);
		karr[i] = (struct ggsea_score_pos_s){
			.score = -aln[i]->score,		/* score in descending order */
			.pos = aln[i]->sec[0].aid + aln[i]->sec[0].bid,
			.idx = i
		};
		debug("pushed, i(%u), score(%lld), pos(%u)",
			karr[i].idx, karr[i].score, karr[i].pos);
	}
	psort_full(karr, cnt, 16, 0);

	for(int64_t i = 0; i < cnt; i++) {
		debug("sorted, i(%lld), idx(%u), score(%lld), pos(%u), score(%lld)",
			i, karr[i].idx, karr[i].score, karr[i].pos, aln[karr[i].idx]->score);
	}

	/* dedup */
	int64_t j = 0;
	for(int64_t i = 0; i < cnt; i++) {
		if(ggsea_cmp_result(aln, karr, i, j) == 0) {
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
 * @fn ggsea_pack_result
 */
static _force_inline
struct ggsea_result_s *ggsea_pack_result(
	struct ggsea_ctx_s *ctx)
{
	/* build array and sort */
	struct gaba_alignment_s const **aln = lmm_kv_ptr(ctx->aln);
	int64_t const cnt = lmm_kv_size(ctx->aln);

	/* dedup result */
	int64_t dedup_cnt = (cnt != 0) ? ggsea_dedup_result(ctx, aln, cnt) : 0;

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
			ggsea_save_seed_pos(ctx, t.kmer, m.gid_pos_arr[0], t.gid_pos);
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
	return((ggsea_result_t *)ggsea_pack_result(ctx));
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

	gref_pool_t *pool = gref_init_pool(GREF_PARAMS( .k = 4 ));
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
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS( .k = 4 ));
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
	gref_pool_t *rpool = gref_init_pool(GREF_PARAMS( .k = 14 ));
	gref_append_segment(rpool, _str("ref1"),
		_seq("CTCACCTCGCTCAAAAGGGCTGCCTCCGAGCGTGTGGGCGAGGACAACCGCCCCACAGTCAAGCTCGAATGGGTGCTATTGCGTAGCTAGGACCGGCACT"));
	gref_idx_t *ref = gref_build_index(gref_freeze_pool(rpool));

	/* build query sequence iterator */
	gref_pool_t *qpool = gref_init_pool(GREF_PARAMS( .k = 14 ));
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

/* adjacent kmer filter */
/**
 * @fn ggsea_update_adjacent_filter
 */
static _force_inline
void ggsea_update_adjacent_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s *arr,
	int64_t len)
{
	debug("save adjacent filter arr(%p), len(%lld)", m.gid_pos_arr, m.len);
	ctx->prev_match_idx = 0;
	ctx->prev_match_len = len;
	ctx->prev_match_arr = arr;
	return;
}

/**
 * @fn ggsea_adjacent_filter
 * @brief filt out previously hit kmer
 */
static _force_inline
int64_t ggsea_adjacent_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s pos)
{
	union gid_pos_u {
		struct gref_gid_pos_s p;
		uint64_t u;
	};
	#define _cast_u(x)	( ((union gid_pos_u){ .p = (x) }).u )

	int64_t i = ctx->prev_match_idx;
	int64_t size = ctx->prev_match_len;
	struct gref_gid_pos_s *arr = ctx->prev_match_arr;
	uint64_t prev_u = 0xffffffff00000000;

	debug("i(%lld), size(%lld), arr(%p), pos(%llx)", i, size, arr, _cast_u(pos));
	while(i < size && (prev_u = _cast_u(arr[i])) < (_cast_u(pos) - 1)) {
		debug("i(%lld), size(%lld), arr[i](%llx), pos(%llx)", i, size, _cast_u(arr[i]), _cast_u(pos));
		i++;
	}
	ctx->prev_match_idx = i;
	debug("adjacent filter, prev_u(%llx), pos(%llx), ret(%d)", prev_u, _cast_u(pos), (prev_u == (_cast_u(pos) - 1)) ? 1 : 0);
	return((prev_u == (_cast_u(pos) - 1)) ? 1 : 0);

	#undef _cast_u
}


/**
 * @struct ggsea_region_s
 */
struct ggsea_region_s {
	ivtree_node_t h;				/* (56) header */
	uint32_t sp, ep;				/* (8) start and end p coordinate */
	int32_t grad, qbase;			/* (8) */
	int32_t depth;					/* (4) */
	int32_t pad[3];					/* (12) */
	int64_t score;					/* (8) */
};
_static_assert(sizeof(struct ggsea_region_s) == 96);


/**
 * @fn ggsea_calc_key
 * @brief calculate q coordinate for use in the overlap filter
 */
static _force_inline
int64_t ggsea_calc_key(
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	union ggsea_q_u {
		int64_t q;
		struct ggsea_q_elems_s {
			uint32_t pos;
			uint32_t gid;
		} e;
	};

	return(((union ggsea_q_u){
		.e = ((struct ggsea_q_elems_s){
			.pos = 0x80000000 ^ (rpos.pos - qpos.pos),
			.gid = rpos.gid ^ ((qpos.gid<<16) | (qpos.gid>>16))
		})
	}).q);
}

/**
 * @fn ggsea_overlap_filter
 */
static _force_inline
int64_t ggsea_overlap_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* calc p (pos) and q (key) */
	int64_t pos = rpos.pos + qpos.pos + ctx->conf.params.k;
	int64_t key = ggsea_calc_key(rpos, qpos);

	debug("pos(%lld), key(%lld)", pos, key);

	/* retrieve intersecting section */
	ivtree_iter_t *iter = ivtree_intersect(
		ctx->tree, key - ctx->envelope_width, key + ctx->envelope_width);

	int64_t depth = INT64_MAX;
	struct ggsea_region_s *n = NULL;
	while((n = (struct ggsea_region_s *)ivtree_next(iter)) != NULL) {
		debug("n(%p), n->h.lkey(%lld), n->h.rkey(%lld), n->sp(%u), n->ep(%u)",
			n, n->h.lkey, n->h.rkey, n->sp, n->ep);

		/* check if the seed is contained in the region */
		if((uint32_t)(pos - n->sp) >= (uint32_t)(n->ep - n->sp)) { continue; }
		if(((n->grad * (int32_t)(pos - n->sp))>>6) + n->qbase - key >= 2 * ctx->envelope_width) {
			continue;
		}

		/* update depth */
		depth = MIN2(depth, n->depth);
	}

	ivtree_iter_clean(iter);
	return((depth == INT64_MAX) ? 0 : depth);
}

/**
 * @fn ggsea_save_overlap_kmer
 */
static _force_inline
void ggsea_save_overlap_kmer(
	struct ggsea_ctx_s *ctx,
	int64_t score,
	struct gref_gid_pos_s rpos,
	struct gref_gid_pos_s qpos)
{
	/* discard (not implemented yet) */
	return;
}

/**
 * @fn ggsea_update_overlap_section
 */
static _force_inline
void ggsea_update_overlap_section(
	struct ggsea_ctx_s *ctx,
	struct ggsea_fill_pair_s pair,
	struct gaba_alignment_s const *r)
{
	debug("update overlap section r->slen(%u)", r->slen);

	for(int64_t i = 0; i < r->slen; i++) {
		debug("update overlap section i(%lld), a(%u, %u, %u), b(%u, %u, %u)",
			i,
			r->sec[i].aid, r->sec[i].apos, r->sec[i].alen,
			r->sec[i].bid, r->sec[i].bpos, r->sec[i].blen);

		/* calc p (pos) and q (key) */
		uint32_t sp = r->sec[i].apos + r->sec[i].bpos;
		int64_t skey = ggsea_calc_key(
			(struct gref_gid_pos_s){ r->sec[i].aid, r->sec[i].apos },
			(struct gref_gid_pos_s){ r->sec[i].bid, r->sec[i].bpos });
		int64_t ekey = ggsea_calc_key(
			(struct gref_gid_pos_s){ r->sec[i].aid, r->sec[i].apos + r->sec[i].alen },
			(struct gref_gid_pos_s){ r->sec[i].bid, r->sec[i].bpos + r->sec[i].blen });

		int64_t lkey = MIN2(skey, ekey);
		int64_t rkey = MAX2(skey, ekey) + 1;
		uint32_t qbase = skey + ctx->envelope_width;
		debug("skey(%lld), ekey(%lld), lkey(%lld), rkey(%lld), qbase(%u)",
			skey, ekey, lkey, rkey, qbase);

		/* search intersecting sections */
		ivtree_iter_t *iter = ivtree_intersect(ctx->tree, lkey, rkey);

		int64_t hit = 0;
		struct ggsea_region_s *n = NULL;
		while((n = (struct ggsea_region_s *)ivtree_next(iter)) != NULL) {
		debug("n(%p), n->h.lkey(%lld), n->h.rkey(%lld), n->sp(%u), n->ep(%u)",
			n, n->h.lkey, n->h.rkey, n->sp, n->ep);

			/* check if the seed is contained in the region */
			if((uint32_t)(sp - n->sp) >= (uint32_t)(n->ep - n->sp)) { continue; }
			if(n->qbase - skey >= 2 * ctx->envelope_width) { continue; }

			/* update node */
			hit++;
			n->depth++;
			n->score = MAX2(n->score, r->score);
		}
		ivtree_iter_clean(iter);

		/* add new region if no region found */
		if(hit == 0) {
			uint32_t len = r->sec[i].alen + r->sec[i].blen;
			struct ggsea_region_s *nn = (struct ggsea_region_s *)ivtree_create_node(ctx->tree);
			*nn = (struct ggsea_region_s){
				.h.lkey = lkey,
				.h.rkey = rkey,
				.sp = sp,
				.ep = sp + len,
				.grad = 64 * (ekey - skey) / len,
				.qbase = qbase,
				.depth = 1,
				.score = r->score
			};

			debug("no hit found, create new region n(%p), n->h.lkey(%lld), n->h.lkey(%lld), n->sp(%u), n->ep(%u)",
				nn, nn->h.lkey, nn->h.rkey, nn->sp, nn->ep);
			ivtree_insert(ctx->tree, (ivtree_node_t *)nn);
		}
	}
	return;
}

/* overlap filters (not implemented yet) */
/**
 * @fn ggsea_update_overlap_filter
 */
static _force_inline
void ggsea_update_overlap_filter(
	struct ggsea_ctx_s *ctx)
{
	return;
}

/**
 * @fn ggsea_overlap_filter
 */
static _force_inline
struct rtree_node_s *ggsea_overlap_filter(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s pos)
{
	struct rtree_node_s *node = ctx->rtree_node;
	int64_t prev_u = INT64_MIN;

	while(node != NULL && (prev_u = _cast_u(node->pos)) < _cast_u(pos)) {
		node = rbtree_right(ctx->rtree, node);
	}
	// ctx->rtree_node = (node == NULL) ? NULL : rbtree_right(ctx->rtree, node);
	ctx->rtree_node = node;

	return((prev_u < _cast_u(pos) + ctx->envelope_width) ? node : NULL);
}

/**
 * @fn ggsea_overlap_next
 */
static _force_inline
struct rtree_node_s *ggsea_overlap_next(
	struct ggsea_ctx_s *ctx,
	struct gref_gid_pos_s pos)
{
	struct rtree_node_s *node = rbtree_right(ctx->rtree, node);
	if(node == NULL) {
		return(NULL);
	}
	return((_cast_u(node->pos) < _cast_u(pos) + ctx->envelope_width) ? node : NULL);
}

#endif
