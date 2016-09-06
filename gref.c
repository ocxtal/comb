
/**
 * @file gref.c
 *
 * @brief reference sequence indexer and searcher.
 *
 * @author Hajime Suzuki
 * @date 2016/3/25
 * @license MIT
 */

#define UNITTEST_UNIQUE_ID			50
#define UNITTEST 					1

#include "unittest.h"

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "hmap.h"
#include "psort.h"
#include "zf.h"
#include "arch/arch.h"
// #include "kvec.h"
#include "lmm.h"
#include "log.h"
#include "sassert.h"
#include "gref.h"


/* constants */
#define K_MIN_BASE					( 2 )
#define K_MIN						( 1<<K_MIN_BASE )
_static_assert(K_MIN_BASE == GREF_K_MIN_BASE);
_static_assert(K_MIN == GREF_K_MIN);

#define K_MAX_BASE					( 5 )
#define K_MAX 						( 1<<K_MAX_BASE )
_static_assert(K_MAX_BASE == GREF_K_MAX_BASE);
_static_assert(K_MAX == GREF_K_MAX);


/* inline directive */
#define _force_inline				inline

/* roundup */
#define _roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )

/* max, min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MAX3(x, y, z)				( MAX2(MAX2(x, y), z) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )
#define MIN3(x, y, z)				( MIN2(MIN2(x, y), z) )

/* encode and decode id */
#define _rev(_d)					( 0x01 ^ (_d) )
#define _encode_id(_x, _d)			( ((_x)<<1) | (0x01 & (_d)) )
#define _decode_id(_x)				( (_x)>>1 )
#define _decode_dir(_x)				( (_x) & 0x01 )


/**
 * structs and typedefs
 */
/* the code supposes NULL == (void *)0 */
// _static_assert(NULL == (void *)0);

/**
 * @struct gref_gid_pair_s
 */
struct gref_gid_pair_s {
	int32_t from;
	int32_t to;
};

/**
 * @struct gref_seq_interval_s
 */
struct gref_seq_interval_s {
	uint64_t base;
	uint64_t tail;
};

/**
 * @struct gref_section_intl_s
 * @brief sizeof(gref_section_intl_s) == 64
 */
struct gref_section_intl_s {
	hmap_header_t header;

	/* forward link index */
	uint32_t fw_link_idx_base;

	/* splitted section link (used to get original name) */
	uint32_t base_gid;

	/* gaba_section_s compatible */
	struct gref_section_s fw_sec;

	/* padding */
	uint64_t reserved1;

	/* reverse link index */
	uint32_t rv_link_idx_base;
	uint32_t reserved2;

	/* gaba_section_s compatible */
	struct gref_section_s rv_sec;
};
_static_assert(sizeof(hmap_header_t) == 8);
_static_assert(sizeof(struct gref_section_intl_s) == 64);

/**
 * @struct gref_section_half_s
 * @brief former half of the section_intl_s
 */
struct gref_section_half_s {
	uint64_t reserved1;
	uint32_t link_idx_base;
	uint32_t reserved2;
	struct gref_section_s sec;
};
_static_assert(sizeof(struct gref_section_half_s) == 32);

/**
 * @enum gref_type
 * @breif gref->type
 */
enum gref_type {
	GREF_POOL 	= 1,
	GREF_ACV 	= 2,
	GREF_IDX 	= 3
};

/**
 * @struct gref_s
 * @brief aliased to gref_pool_t, gref_acv_t, and gref_idx_t
 */
struct gref_s {
	lmm_t *lmm;

	/* name -> section mapping */
	hmap_t *hmap;					/* name -> section_info hashmap */
	uint32_t sec_cnt;

	/* status */
	int8_t type;
	int8_t kmer_available;
	uint8_t reserved2[2];

	/* internal params */
	int64_t iter_init_stack_size;

	/* params */
	struct gref_params_s params;

	/* sequence container */
	lmm_kvec_t(uint8_t) seq;
	uint64_t seq_len;
	uint8_t const *seq_lim;

	/* link info container */
	lmm_kvec_t(struct gref_gid_pair_s) link;

	/* lmm_kv_ptr(link) and link_table shares pointer */
	uint64_t mask;
	int64_t link_table_size;
	uint32_t *link_table;

	/* kmer index container */
	int64_t *kmer_idx_table;
	int64_t kmer_table_size;
	struct gref_gid_pos_s *kmer_table;

	/* sequence encoder */
	struct gref_seq_interval_s (*append_seq)(
		struct gref_s *gref,
		uint8_t const *seq,
		int64_t len);
};
_static_assert(sizeof(struct gref_params_s) == 24);

/**
 * @fn gref_encode_2bit
 * @brief mapping IUPAC amb. to 2bit encoding
 */
static _force_inline
uint8_t gref_encode_2bit(
	int c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x00, C = 0x01, G = 0x02, T = 0x03
	};
	static uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('N')] = A,		/* treat 'N' as 'A' */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn gref_encode_4bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t gref_encode_4bit(
	uint8_t c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x01, C = 0x02, G = 0x04, T = 0x08
	};
	static uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('R')] = A | G,
		[_b('Y')] = C | T,
		[_b('S')] = G | C,
		[_b('W')] = A | T,
		[_b('K')] = G | T,
		[_b('M')] = A | C,
		[_b('B')] = C | G | T,
		[_b('D')] = A | G | T,
		[_b('H')] = A | C | T,
		[_b('V')] = A | C | G,
		[_b('N')] = 0,		/* treat 'N' as a gap */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b(c)]);

	#undef _b
}

/**
 * @fn gref_copy_seq_ascii
 */
static
struct gref_seq_interval_s gref_copy_seq_ascii(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	uint64_t base = lmm_kv_size(gref->seq);
	lmm_kv_reserve(gref->lmm, gref->seq, base + len);

	/* append */
	for(int64_t i = 0; i < len; i++) {
		lmm_kv_at(gref->seq, base + i) = gref_encode_4bit(seq[i]);
	}

	/* resize array */
	lmm_kv_size(gref->seq) += len;
	return((struct gref_seq_interval_s){
		.base = base,
		.tail = base + len
	});	
}

/**
 * @fn gref_copy_seq_4bit
 */
static
struct gref_seq_interval_s gref_copy_seq_4bit(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	uint64_t base = lmm_kv_size(gref->seq);
	lmm_kv_pushm(gref->lmm, gref->seq, seq, len);
	return((struct gref_seq_interval_s){
		.base = base,
		.tail = base + len
	});
}

/**
 * @fn gref_nocopy_seq_4bit
 */
static
struct gref_seq_interval_s gref_nocopy_seq_4bit(
	struct gref_s *gref,
	uint8_t const *seq,
	int64_t len)
{
	return((struct gref_seq_interval_s){
		.base = (uint64_t)seq,
		.tail = (uint64_t)seq + len
	});
}

/* init / destroy pool */
/**
 * @fn gref_init_pool
 */
gref_pool_t *gref_init_pool(
	gref_params_t const *params)
{
	struct gref_params_s const default_params = { 0 };
	struct gref_params_s p = (params == NULL) ? default_params : *params;

	/* restore defaults */
	#define restore(param, def)		{ (param) = ((uint64_t)(param) == 0) ? (def) : (param); }
	
	restore(p.k, 14);
	restore(p.seq_direction, GREF_FW_ONLY);
	restore(p.seq_format, GREF_ASCII);
	restore(p.copy_mode, GREF_COPY);
	restore(p.num_threads, 0);
	restore(p.hash_size, 1024);
	restore(p.seq_head_margin, 0);
	restore(p.seq_tail_margin, 0);
	restore(p.lmm, NULL);

	#undef restore

	/* check sanity */
	if(p.k < K_MIN || p.k > K_MAX) { return(NULL); }
	if((uint8_t)p.seq_format > GREF_4BIT) { return(NULL); }
	if((uint8_t)p.copy_mode > GREF_NOCOPY) { return(NULL); }
	p.seq_head_margin = _roundup(p.seq_head_margin, 16);
	p.seq_tail_margin = _roundup(p.seq_tail_margin, 16);

	debug("init called with seq_direction(%x), seq_format(%x), copy_mode(%x), "
		"seq_head_margin(%x), seq_tail_margin(%x), num_threads(%x), hash_size(%x)",
		p.seq_direction, p.seq_format, p.copy_mode,
		p.seq_head_margin, p.seq_tail_margin, p.num_threads, p.hash_size);

	/* malloc mem */
	lmm_t *lmm = (lmm_t *)p.lmm;
	struct gref_s *pool = (struct gref_s *)lmm_malloc(lmm, sizeof(struct gref_s));
	if(pool == NULL) {
		return(NULL);
	}
	memset(pool, 0, sizeof(struct gref_s));
	pool->lmm = lmm;

	/* init hmap */
	pool->hmap = hmap_init(
		sizeof(struct gref_section_intl_s),
		HMAP_PARAMS( .hmap_size = p.hash_size, .lmm = lmm ));
	if(pool->hmap == NULL) {
		goto _gref_init_pool_error_handler;
	}
	pool->sec_cnt = 0;
	pool->type = GREF_POOL;

	/* calc iterator buffer size */
	int64_t buf_size = 1;
	for(int64_t i = 0; i < (p.k + 1) / 2; i++) {
		buf_size *= 3;
	}
	pool->iter_init_stack_size = MAX2(1024, buf_size);

	/* init seq vector */
	if(p.copy_mode != GREF_NOCOPY) {
		lmm_kv_init(lmm, pool->seq);

		/* make margin at the head (leave uninitialized) */
		lmm_kv_reserve(lmm, pool->seq, p.seq_head_margin);
		lmm_kv_size(pool->seq) = p.seq_head_margin;
	} else {
		lmm_kv_ptr(pool->seq) = NULL;

		/* ignore margin option in nocopy mode */
		p.seq_head_margin = 0;
		p.seq_tail_margin = 0;
	}
	pool->seq_len = 0;
	pool->seq_lim = NULL;

	/* init link vector */
	lmm_kv_init(lmm, pool->link);

	/* init seq encoder */
	struct gref_seq_interval_s (*table[][3])(
		struct gref_s *gref,
		uint8_t const *seq,
		int64_t len) = {
		[GREF_ASCII] = {
			[GREF_COPY] = gref_copy_seq_ascii,
			[GREF_NOCOPY] = NULL		/* nocopy mode with ascii input is not supported */
		},
		[GREF_4BIT] = {
			[GREF_COPY] = gref_copy_seq_4bit,
			[GREF_NOCOPY] = gref_nocopy_seq_4bit
		}
	};
	if((pool->append_seq = table[p.seq_format][p.copy_mode]) == NULL) {
		goto _gref_init_pool_error_handler;
	}

	/* copy params */
	pool->params = p;
	return((gref_pool_t *)pool);

_gref_init_pool_error_handler:;
	if(pool != NULL) {
		hmap_clean(pool->hmap); pool->hmap = NULL;
		lmm_kv_destroy(lmm, pool->seq);
		lmm_kv_destroy(lmm, pool->link);
		lmm_free(lmm, pool);
	}
	return(NULL);
}

/**
 * @fn gref_clean
 */
void gref_clean(
	gref_t *_gref)
{
	struct gref_s *gref = (struct gref_s *)_gref;

	if(gref != NULL) {
		/* cleanup, cleanup... */
		hmap_clean(gref->hmap); gref->hmap = NULL;
		lmm_kv_destroy(gref->lmm, gref->seq);
		lmm_kv_destroy(gref->lmm, gref->link);
		// free(gref->link_table); gref->link_table = NULL;
		lmm_free(gref->lmm, gref->kmer_idx_table); gref->kmer_idx_table = NULL;
		lmm_free(gref->lmm, gref->kmer_table); gref->kmer_table = NULL;
		lmm_free(gref->lmm, gref);
	}
	return;
}


/* pool modify operation */
/**
 * @fn gref_append_segment
 */
int gref_append_segment(
	gref_pool_t *_pool,
	char const *name,
	int32_t name_len,
	uint8_t const *seq,
	int64_t seq_len)
{
	struct gref_s *pool = (struct gref_s *)_pool;
	// debug("append segment");

	/* gref object is mutable only when type == POOL */
	if(pool == NULL || pool->type != GREF_POOL) { return(-1); }

	/* add sequence at the tail of the seq buffer */
	struct gref_seq_interval_s iv = pool->append_seq(pool, seq, seq_len);

	/* update length */
	pool->seq_len += iv.tail - iv.base;

	/* append the first section */
	uint64_t const max_sec_len = 0x80000000;
	uint64_t len = MIN2(iv.tail - iv.base, max_sec_len);

	uint32_t id = hmap_get_id(pool->hmap, name, name_len);
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, id);

	/* update sec_cnt */
	pool->sec_cnt = MAX2(pool->sec_cnt, id + 1);

	/* store section info */
	sec->base_gid = _encode_id(id, 0);
	sec->fw_link_idx_base = 0;
	sec->rv_link_idx_base = 0;
	sec->fw_sec = (struct gref_section_s){
		.gid = _encode_id(id, 0),
		.len = len,
		.base = (void *)(iv.base - pool->params.seq_head_margin)
	};
	sec->rv_sec = (struct gref_section_s){
		.gid = _encode_id(id, 1),
		.len = len,
		.base = NULL
	};
	return(0);
}

/**
 * @fn gref_append_link
 */
int gref_append_link(
	gref_pool_t *_pool,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori)
{
	struct gref_s *pool = (struct gref_s *)_pool;
	debug("append link, src(%s), dst(%s), cnt(%u)", src, dst, hmap_get_count(pool->hmap));

	/* gref object is mutable only when type == POOL */
	if(pool == NULL || pool->type != GREF_POOL) { return(-1); }

	/* get ids */
	uint32_t src_id = hmap_get_id(pool->hmap, src, src_len);
	uint32_t dst_id = hmap_get_id(pool->hmap, dst, dst_len);

	/* add forward link */
	lmm_kv_push(pool->lmm, pool->link, ((struct gref_gid_pair_s){
		.from = _encode_id(src_id, src_ori),
		.to = _encode_id(dst_id, dst_ori)
	}));

	/* add reverse link */
	lmm_kv_push(pool->lmm, pool->link, ((struct gref_gid_pair_s){
		.from = _encode_id(dst_id, _rev(dst_ori)),
		.to = _encode_id(src_id, _rev(src_ori))
	}));

	/* update tail_id */
	pool->sec_cnt = MAX3(pool->sec_cnt, src_id + 1, dst_id + 1);
	debug("sec_cnt(%u), src_id(%u), dst_id(%u)", pool->sec_cnt, src_id, dst_id);
	return(0);
}

/**
 * @fn gref_append_snp
 */
int gref_append_snp(
	gref_pool_t *_pool,
	char const *name,
	int32_t name_len,
	int64_t pos,
	uint8_t snp)
{
	struct gref_s *pool = (struct gref_s *)_pool;

	/* gref object is mutable only when type == POOL */
	if(pool == NULL || pool->type != GREF_POOL) { return(-1); }

	/* fixme: not implemented yet */
	return(0);
}

/**
 * @fn gref_split_segment
 * @brief split base section and give new name (splitted) to the latter section.
 */
int gref_split_segment(
	gref_pool_t *_pool,
	char const *base,
	int32_t base_len,
	int64_t pos,
	char const *splitted,
	int32_t splitted_len)
{
	struct gref_s *pool = (struct gref_s *)_pool;

	/* gref object is mutable only when type == POOL */
	if(pool == NULL || pool->type != GREF_POOL) { return(-1); }

	/* fixme: not implemented yet */
	return(0);
}

/**
 * @fn gref_merge_pools
 * @brief merge two pools
 */
gref_pool_t *gref_merge_pools(
	gref_pool_t *_pool1,
	gref_pool_t *_pool2)
{
	/* not implemented yet */
	return(NULL);
}


/* build link table (pool -> acv conversion) */

/**
 * @fn gref_add_tail_section
 */
static _force_inline
void gref_add_tail_section(
	struct gref_s *pool)
{
	uint32_t tail_id = pool->sec_cnt;
	if(hmap_get_count(pool->hmap) >= tail_id) {
		/* sentinel already exists */
		return;
	}

	/* push sentinel with unique name */
	char const *template = "tail_sentinel_";
	int64_t len = strlen(template);

	char buf[256];
	strcpy(buf, template);
	uint32_t id = (uint32_t)-1;
	do {
		buf[len++] = '0';
		buf[len] = '\0';

		/* push sentinel to section array */
		id = hmap_get_id(pool->hmap, buf, len);
	} while(id != tail_id && len < 256);

	/* make margin at the tail (leave uninitialized) */
	lmm_kv_reserve(pool->lmm, pool->seq, lmm_kv_size(pool->seq) + pool->params.seq_tail_margin);
	lmm_kv_size(pool->seq) += pool->params.seq_tail_margin;

	/* set info */
	struct gref_section_intl_s *tail_sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, tail_id);
	tail_sec->base_gid = _encode_id(tail_id, 0);
	tail_sec->fw_sec = (struct gref_section_s){
		.gid = _encode_id(tail_id, 0),
		.len = 0,
		.base = 0
	};
	tail_sec->rv_sec = (struct gref_section_s){
		.gid = _encode_id(tail_id, 1),
		.len = 0,
		.base = 0
	};
	return;
}

/**
 * @fn gref_build_link_idx_table
 * @brief build link_idx table. gref_add_tail_section must be called beforehand.
 */
static _force_inline
int gref_build_link_idx_table(
	struct gref_s *pool)
{
	int64_t link_idx_table_size = 2 * pool->sec_cnt;
	int64_t link_table_size = lmm_kv_size(pool->link);

	/* store link table size */
	pool->link_table_size = link_table_size;

	/* sort by src, build src->dst mapping */
	debug("sort src->dst mapping, size(%llu)", lmm_kv_size(pool->link));
	if(psort_half(lmm_kv_ptr(pool->link), lmm_kv_size(pool->link),
		sizeof(struct gref_gid_pair_s), pool->params.num_threads) != 0) {

		/* sort failed */
		return(-1);
	}

	/* init index */
	uint32_t prev_gid = 0;
	struct gref_section_half_s *sec_half =
		(struct gref_section_half_s *)hmap_get_object(pool->hmap, 0);

	/* store link index */
	sec_half[prev_gid].link_idx_base = 0;		/* head is always zero */
	for(int64_t i = 0; i < link_table_size; i++) {
		uint32_t gid = lmm_kv_at(pool->link, i).from;
		debug("i(%lld), gid(%u), prev_gid(%u)", i, gid, prev_gid);

		if(prev_gid == gid) { continue; }

		debug("index table for gid(%u) ends at %lld, next_gid(%u)", prev_gid, i, gid);

		/* sequence update detected */
		for(int64_t j = prev_gid + 1; j < gid + 1; j++) {
			debug("fill gaps j(%lld)", j);
			sec_half[j].link_idx_base = i;
		}
		prev_gid = gid;
	}

	/* store tail */
	for(int64_t j = prev_gid + 1; j < link_idx_table_size + 1; j++) {
		// debug("fill gaps j(%lld)", j);
		sec_half[j].link_idx_base = link_table_size;
	}
	return(0);
}

/**
 * @fn gref_*_*_modify_seq
 * @brief calculate pos and lim of reverse section, append rv seq if needed
 */
static
int gref_fw_copy_modify_seq(
	struct gref_s *pool)
{
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, 0);

	debug("gref_fw_copy_modify_seq called");

	/* lim == 0x800000000000 in fw_copy mode */
	pool->seq_lim = GREF_SEQ_LIM;
	uint8_t const *rv_lim = GREF_SEQ_LIM + (uint64_t)GREF_SEQ_LIM;

	/* add offset to fw pos base, then calc rv pos base */
	uint64_t seq_base = (uint64_t)lmm_kv_ptr(pool->seq) + pool->params.seq_head_margin;
	for(int64_t i = 0; i < pool->sec_cnt; i++) {
		sec[i].fw_sec.base += seq_base;			/* convert pos to valid pointer */
		sec[i].rv_sec.base = rv_lim - (uint64_t)sec[i].fw_sec.base - sec[i].fw_sec.len;
	}
	return(0);
}
static
int gref_fw_nocopy_modify_seq(
	struct gref_s *pool)
{
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, 0);

	debug("gref_fw_nocopy_modify_seq called");

	/* lim == 0x800000000000 in fw_copy mode */
	pool->seq_lim = GREF_SEQ_LIM;
	uint8_t const *rv_lim = GREF_SEQ_LIM + (uint64_t)GREF_SEQ_LIM;

	/* calc rv pos base */
	for(uint64_t i = 0; i < pool->sec_cnt; i++) {
		sec[i].rv_sec.base = rv_lim - (uint64_t)sec[i].fw_sec.base - sec[i].fw_sec.len;
	}
	return(0);
}
static
int gref_fr_copy_modify_seq(
	struct gref_s *pool)
{
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, 0);

	debug("gref_fr_copy_modify_seq called");

	/* append reverse sequence */
	// uint64_t size = lmm_kv_size(pool->seq);
	uint64_t size = 2 * pool->seq_len
		+ pool->params.seq_head_margin
		+ pool->params.seq_tail_margin;
	lmm_kv_reserve(pool->lmm, pool->seq, size);
	if(lmm_kv_ptr(pool->seq) == NULL) { return(-1); }

	/* append revcomp seq */
	static uint8_t const comp[16] = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	uint64_t fw_tail_pos = pool->seq_len + pool->params.seq_head_margin;
	for(uint64_t i = 0; i < pool->seq_len; i++) {
		lmm_kv_at(pool->seq, fw_tail_pos + i) = comp[lmm_kv_at(pool->seq, fw_tail_pos - 1 - i)];
	}

	/* lim == 0x800000000000 in fw_copy mode */
	uint8_t const *seq_base = lmm_kv_ptr(pool->seq) + pool->params.seq_head_margin;
	uint8_t const *rv_lim = pool->seq_lim = seq_base + 2 * pool->seq_len;

	/* add offset to fw pos base, then calc rv pos base */
	for(uint64_t i = 0; i < pool->sec_cnt; i++) {
		sec[i].rv_sec.base = rv_lim - (uint64_t)sec[i].fw_sec.base - sec[i].fw_sec.len;
		sec[i].fw_sec.base += (uint64_t)seq_base;	/* convert pos to valid pointer */
		debug("%llu", (uint64_t)(sec[i].rv_sec.base - sec[i].fw_sec.base));
	}
	return(0);
}
static
int gref_fr_nocopy_modify_seq(
	struct gref_s *pool)
{
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(pool->hmap, 0);

	debug("gref_fr_nocopy_modify_seq called");

	/* lim == 0x800000000000 in fw_copy mode */
	pool->seq_lim = GREF_SEQ_LIM;

	/* calc rv pos base */
	for(int64_t i = 0; i < pool->sec_cnt; i++) {
		sec[i].rv_sec.base = sec[i].fw_sec.base + sec[i].fw_sec.len;
	}
	return(0);
}

/**
 * @fn gref_flush_modified_seq
 */
static _force_inline
int gref_flush_modified_seq(
	struct gref_s *acv)
{
	if(acv->params.copy_mode == GREF_NOCOPY) {
		lmm_kv_ptr(acv->seq) = NULL;
		return(0);
	}

	/* resize */
	lmm_kv_resize(acv->lmm, acv->seq, acv->seq_len);
	lmm_kv_size(acv->seq) = acv->seq_len;
	if(lmm_kv_ptr(acv->seq) == NULL) { return(-1); }

	/* subtract seq_base from fw_sec, clear rv_sec with NULL */
	uint64_t seq_base = (uint64_t)lmm_kv_ptr(acv->seq) + acv->params.seq_head_margin;
	struct gref_section_intl_s *sec =
		(struct gref_section_intl_s *)hmap_get_object(acv->hmap, 0);
	for(int64_t i = 0; i < 2 * (acv->sec_cnt + 1); i++) {
		sec[i].fw_sec.base -= seq_base;
		sec[i].rv_sec.base = NULL;
	}
	return(0);
}

/**
 * @fn gref_shrink_link_table
 * @brief shrink link table. gref_build_link_idx_table must be called beforehand.
 */
static _force_inline
int gref_shrink_link_table(
	struct gref_s *pool)
{
	int64_t link_table_size = pool->link_table_size;
	uint32_t *link_table = (uint32_t *)lmm_kv_ptr(pool->link);

	/* pack */
	for(int64_t i = 0; i < link_table_size; i++) {
		link_table[i] = lmm_kv_at(pool->link, i).to;
	}
	lmm_kv_resize(pool->lmm, pool->link, link_table_size / 2);
	lmm_kv_size(pool->link) = link_table_size / 2;
	debug("shrink ptr(%p), size(%llu)", lmm_kv_ptr(pool->link), link_table_size / 2);
	
	/* store info */
	pool->link_table = (uint32_t *)lmm_kv_ptr(pool->link);
	return((pool->link_table == NULL) ? -1 : 0);
}

/**
 * @fn gref_expand_link_table
 */
static _force_inline
int gref_expand_link_table(
	struct gref_s *acv)
{
	/* resize mem */
	int64_t link_table_size = acv->link_table_size;
	lmm_kv_resize(acv->lmm, acv->link, link_table_size);
	uint32_t *link = (uint32_t *)lmm_kv_ptr(acv->link);

	if(link == NULL) { return(-1); }

	/* load section ptr */
	int64_t sec_cnt = acv->sec_cnt;
	struct gref_section_half_s *sec_half =
		(struct gref_section_half_s *)hmap_get_object(acv->hmap, 0);

	/* expand, iterate from tail to head */
	for(int64_t i = 2 * sec_cnt - 1; i >= 0; i--) {
		for(int64_t j = sec_half[i + 1].link_idx_base - 1;
			j >= sec_half[i].link_idx_base;
			j--) {

			lmm_kv_at(acv->link, j) = (struct gref_gid_pair_s){
				.from = i,
				.to = link[j]
			};
		}
	}
	return(0);
}

/**
 * @fn gref_freeze_pool
 */
gref_acv_t *gref_freeze_pool(
	gref_pool_t *pool)
{
	struct gref_s *gref = (struct gref_s *)pool;

	debug("gref_freeze_pool called, gref(%p)", gref);

	if(gref == NULL || gref->type != GREF_POOL) {
		goto _gref_freeze_pool_error_handler;
	}

	/* push tail sentinel */
	gref_add_tail_section(gref);

	/* modify seq */
	int (*modify_seq[][3])(struct gref_s *pool) = {
		[GREF_FW_ONLY] = {
			[GREF_COPY] = gref_fw_copy_modify_seq,
			[GREF_NOCOPY] = gref_fw_nocopy_modify_seq
		},
		[GREF_FW_RV] = {
			[GREF_COPY] = gref_fr_copy_modify_seq,
			[GREF_NOCOPY] = gref_fr_nocopy_modify_seq
		}
	};
	if(modify_seq[gref->params.seq_direction][gref->params.copy_mode](gref) != 0) {
		goto _gref_freeze_pool_error_handler;
	}

	/* build link array */
	if(gref_build_link_idx_table(gref) != 0) {
		/* sort failed */
		goto _gref_freeze_pool_error_handler;
	}

	/* shrink table */
	if(gref_shrink_link_table(gref) != 0) {
		goto _gref_freeze_pool_error_handler;
	}

	/* change type */
	gref->type = GREF_ACV;
	return((gref_acv_t *)gref);

_gref_freeze_pool_error_handler:;
	gref_clean((void *)gref);
	return(NULL);
}

/**
 * @fn gref_melt_archive
 */
gref_pool_t *gref_melt_archive(
	gref_acv_t *acv)
{
	struct gref_s *gref = (struct gref_s *)acv;

	if(gref == NULL || gref->type != GREF_ACV) {
		goto _gref_melt_archive_error_handler;
	}

	/* flush modifications */
	if(gref_flush_modified_seq(gref) != 0) {
		goto _gref_melt_archive_error_handler;
	}

	/* remove kmer table */
	lmm_free(acv->lmm, gref->kmer_table); gref->kmer_table = NULL;
	gref->kmer_table_size = 0;

	/* expand table */
	if(gref_expand_link_table(gref) != 0) {
		goto _gref_melt_archive_error_handler;
	}

	/* change type */
	gref->type = GREF_POOL;
	return((gref_pool_t *)gref);

_gref_melt_archive_error_handler:;
	gref_clean((void *)gref);
	return(NULL);
}


/* kmer enumeration */
/**
 * @struct gref_iter_kmer_s
 * @brief kmer container
 */
struct gref_iter_kmer_s {
	uint8_t vac_len;
	uint8_t shift_len;
	uint8_t init_len;
	uint8_t pad;

	/* kmer table */
	uint16_t idx;
	uint16_t lim;
	uint64_t cnt;
	// uint64_t arr[];
};
_static_assert(sizeof(struct gref_iter_kmer_s) == 16);
#define _kmer_size(x)			( sizeof(struct gref_iter_kmer_s) + (x).lim * sizeof(uint64_t) )
// #define _kmer_tail(x)			( &(x).arr[(x).lim] )
// #define _kmer_tail(x)			( (struct gref_iter_kmer_s *)(x) + 1 )
#define _kmer_arr(x)			( (uint64_t *)((struct gref_iter_kmer_s *)(x) + 1) )
#define _kmer_tail(x)			( _kmer_arr(x) + (x)->lim )

/**
 * @fn gref_iter_kmer_flush
 */
static _force_inline
void gref_iter_kmer_flush(
	struct gref_iter_kmer_s *kmer)
{
	kmer->vac_len = kmer->init_len;
	kmer->idx = 1;
	kmer->lim = 1;
	kmer->cnt = 0;
	_kmer_arr(kmer)[0] = 0;

	debug("flush kmer, vac_len(%u), idx(%u), lim(%u)", kmer->vac_len, kmer->idx, kmer->lim);
	return;
}

/**
 * @fn gref_iter_kmer_init
 */
static _force_inline
void gref_iter_kmer_init(
	struct gref_iter_kmer_s *kmer,
	int64_t k)
{
	kmer->shift_len = 2 * (k - 1);
	kmer->init_len = k;
	gref_iter_kmer_flush(kmer);
	return;
}

/**
 * @fn gref_iter_kmer_append
 */
static _force_inline
void gref_iter_kmer_append(
	struct gref_iter_kmer_s *kmer,
	uint8_t conv,
	uint8_t c)
{
	if(c == 0) {
		gref_iter_kmer_flush(kmer);
		return;
	}

	/* conversion tables */
	static uint8_t const popcnt_table[] = {
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 0	/* ignore 0x0f */
	};
	enum shift_size {
		A = 0x00, C = 0x02, G = 0x04, T = 0x06
	};
	static uint8_t const shift_table[][3] = {
		{ 0 },
		{ A },
		{ C },
		{ A, C },
		{ G },
		{ A, G },
		{ C, G },
		{ A, C, G },
		{ T },
		{ A, T },
		{ C, T },
		{ A, C, T },
		{ G, T },
		{ A, G, T },
		{ C, G, T },
		{ 0 },
	};

	/* update count array */
	uint64_t pcnt = popcnt_table[c];
	kmer->cnt = (kmer->cnt>>2) | (pcnt<<(kmer->shift_len + 2));
	uint64_t lim = kmer->lim;

	/* branch */
	switch(3 - pcnt) {
		case 0: memcpy(&_kmer_arr(kmer)[2 * lim], _kmer_arr(kmer), sizeof(uint64_t) * lim);
		case 1: memcpy(&_kmer_arr(kmer)[lim], _kmer_arr(kmer), sizeof(uint64_t) * lim);
		/* fall through */
		default: break;		/* return(-1); */
	}

	/* append to vector */
	uint64_t *p = _kmer_arr(kmer);
	uint64_t mask = 0x03<<kmer->shift_len;
	for(uint64_t j = 0; j < pcnt; j++) {
		uint64_t b = mask & (conv<<(kmer->shift_len - shift_table[c][j]));
		for(uint64_t k = 0; k < lim; k++) {
			*p = (*p>>2) | b; p++;
			debug("%lld, %lld, %lld, %x, %x, %llx",
				j, k, j * lim + k, shift_table[c][j], 0x03 & (conv>>shift_table[c][j]), p[-1]);
		}
	}

	/* update lim */
	lim *= pcnt;

	/* merge (shrink buffer) */
	uint64_t shrink_skip = 0x03 & kmer->cnt;
	if(shrink_skip > 1) {
		// lim /= shrink_skip;
		lim = (lim * ((shrink_skip == 2) ? 0x10000 : 0xaaab))>>17;
		for(uint64_t j = 0; j < lim; j++) {
			_kmer_arr(kmer)[j] = _kmer_arr(kmer)[j * shrink_skip];
		}
	}

	/* write back lim */
	kmer->vac_len -= (kmer->vac_len != 0);
	kmer->idx = (kmer->vac_len != 0) ? lim : 0;
	kmer->lim = lim;

	debug("vac_len(%u), cnt(%llx), lim(%llu), shrink_skip(%llu)", kmer->vac_len, kmer->cnt, lim, shrink_skip);
	return;
}

/**
 * @fn gref_iter_kmer_ready
 */
static _force_inline
uint64_t gref_iter_kmer_ready(
	struct gref_iter_kmer_s *kmer)
{
	debug("idx(%u), lim(%u)", kmer->idx, kmer->lim);
	return(kmer->idx < kmer->lim);
}

/**
 * @fn gref_iter_kmer_next
 */
static _force_inline
uint64_t gref_iter_kmer_next(
	struct gref_iter_kmer_s *kmer)
{
	return(_kmer_arr(kmer)[kmer->idx++]);
}

/**
 * @struct gref_iter_stack_s
 * @brief aliased to gref_iter_t
 */
struct gref_iter_stack_s {
	/* previous stack */
	struct gref_iter_stack_s *prev_stack;

	/* current section id */
	uint32_t sec_gid;

	/* sequence */
	uint32_t len;
	uint8_t const *seq_ptr;
	uint32_t rem_len;
	int8_t incr;
	uint8_t conv_table;
	uint8_t global_rem_len;

	/* link */
	uint8_t link_depth;
	uint32_t link_ridx, link_len;
	uint32_t const *link_base;

	/* kmer table */
	struct gref_iter_kmer_s kmer;
};
_static_assert(sizeof(struct gref_iter_stack_s) == 64);

/**
 * @struct gref_iter_s
 */
#define GREF_ITER_INTL_MEM_ARR_LEN			( 5 )
struct gref_iter_s {
	lmm_t *lmm;

	/* global info */
	uint32_t base_gid;
	uint32_t tail_gid;
	uint32_t step_gid;
	uint8_t seed_len;
	uint8_t shift_len;
	uint16_t reserved2;

	uint8_t const *seq_lim;
	uint32_t const *link_table;
	struct gref_section_half_s const *hsec;

	/* stack mem array */
	struct gref_iter_stack_s *stack;
	void *mem_arr[GREF_ITER_INTL_MEM_ARR_LEN];
};
_static_assert(sizeof(struct gref_iter_s) == 96);

/**
 * @fn gref_iter_init_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_init_stack(
	struct gref_iter_s *iter)
{
	/* create stack */
	struct gref_iter_stack_s *stack = (struct gref_iter_stack_s *)(iter + 1);
	iter->stack = stack;

	uint32_t gid = iter->base_gid;
	uint32_t len = iter->hsec[gid].sec.len;
	debug("init_stack called, stack(%p), gid(%u), len(%u)", stack, gid, len);

	/* init stack */
	stack->prev_stack = NULL;

	/* current section info */
	stack->sec_gid = gid;

	/* seq info */
	uint8_t const *base = iter->hsec[gid].sec.base;

	stack->len = len - iter->seed_len;
	stack->seq_ptr = (base < iter->seq_lim)
		? base : (iter->seq_lim + (iter->seq_lim - base - 1));
	stack->rem_len = len;
	stack->incr = (base < iter->seq_lim) ? 1 : -1;
	stack->conv_table = (base < iter->seq_lim) ? 0xe4 : 0x1b;
	stack->global_rem_len = iter->seed_len - 1;

	/* link info */
	int64_t link_idx_base = iter->hsec[gid].link_idx_base;
	int64_t link_idx_tail = iter->hsec[gid + 1].link_idx_base;
	uint32_t link_len = (uint32_t)(link_idx_tail - link_idx_base);

	stack->link_depth = K_MAX_BASE;
	stack->link_ridx = link_len;
	stack->link_len = link_len;
	stack->link_base = iter->link_table + link_idx_base;

	debug("init_stack finished, stack(%p), rem_len(%u)", stack, stack->rem_len);
	return(stack);
}

/**
 * @fn gref_iter_push_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_push_stack(
	struct gref_iter_s *iter,
	struct gref_iter_stack_s *stack)
{
	/* fixme: check remaining memory size and malloc if there is no room */
	// struct gref_iter_stack_s *new_stack = (struct gref_iter_stack_s *)((uint64_t *)(stack + 1) + stack->lim);
	struct gref_iter_stack_s *new_stack = (struct gref_iter_stack_s *)_kmer_tail(&stack->kmer);

	debug("stack(%p), lim(%u), new_stack(%p)", stack, stack->kmer.lim, new_stack);

	/* load next gid, adjust len for path encoding */
	int64_t link_idx = stack->link_len - stack->link_ridx;
	uint32_t gid = stack->link_base[link_idx];
	stack->len += link_idx<<stack->link_depth;
	stack->link_ridx--;

	/* link to previous stack */
	new_stack->prev_stack = stack;

	/* section id */
	new_stack->sec_gid = gid;

	/* init seq info */
	uint8_t const *base = iter->hsec[gid].sec.base;
	uint32_t rem_len = MIN2(stack->global_rem_len, iter->hsec[gid].sec.len);
	uint32_t global_rem_len = stack->global_rem_len - rem_len;
	debug("rem_len(%u), stack->global_rem_len(%u), sec.len(%u), global_rem_len(%u)",
		rem_len, stack->global_rem_len, iter->hsec[gid].sec.len, global_rem_len);

	new_stack->len = stack->len + rem_len;
	new_stack->seq_ptr = (base < iter->seq_lim)
		? base : (iter->seq_lim + (iter->seq_lim - base - 1));
	new_stack->rem_len = rem_len;
	new_stack->incr = (base < iter->seq_lim) ? 1 : -1;
	new_stack->conv_table = (base < iter->seq_lim) ? 0xe4 : 0x1b;
	new_stack->global_rem_len = global_rem_len;

	/* link info */
	int64_t link_idx_base = iter->hsec[gid].link_idx_base;
	int64_t link_idx_tail = iter->hsec[gid + 1].link_idx_base;
	uint32_t link_len = (uint32_t)(link_idx_tail - link_idx_base);

	new_stack->link_depth = stack->link_depth + K_MAX_BASE;
	new_stack->link_ridx = link_len;
	new_stack->link_len = link_len;
	new_stack->link_base = iter->link_table + link_idx_base;

	/* copy kmer buffer */
	memcpy(&new_stack->kmer, &stack->kmer, _kmer_size(stack->kmer));
	return(new_stack);
}

/**
 * @fn gref_iter_pop_stack
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_pop_stack(
	struct gref_iter_stack_s *stack)
{
	debug("remove_stack called stack(%p)", stack);
	return(stack->prev_stack);
}

/**
 * @fn gref_iter_fetch_base
 */
static _force_inline
uint8_t gref_iter_fetch_base(
	struct gref_iter_stack_s *stack)
{
	uint8_t c = *stack->seq_ptr;
	stack->rem_len--;
	stack->seq_ptr += stack->incr;
	return(c);
}

/**
 * @fn gref_iter_fetch_link
 */
static _force_inline
struct gref_iter_stack_s *gref_iter_fetch_link(
	struct gref_iter_s *iter,
	struct gref_iter_stack_s *stack)
{
	/* return if no more seq remains */
	if(stack->global_rem_len == 0) {
		debug("reached tail");
		// stack = stack->prev_stack;
		stack = gref_iter_pop_stack(stack);
	}

	/* remove stack if no more link remains */
	while(stack->link_ridx == 0) {
		if((stack = gref_iter_pop_stack(stack)) == NULL) {
			/* root (reached the end of enumeration) */
			debug("reached NULL");
			return(NULL);
		}
	}

	/* add stack for the next section */
	stack = gref_iter_push_stack(iter, stack);
	debug("stack added, gid(%u), seq_ptr(%p), len(%u), rem_len(%u), global_rem_len(%u)",
		stack->sec_gid, stack->seq_ptr, stack->len, stack->rem_len, stack->global_rem_len);

	/* fetch seq */
	// return(gref_iter_fetch_base(iter, stack));
	return(stack);
}

/**
 * @fn gref_iter_fetch
 */
static _force_inline
uint8_t gref_iter_fetch(
	struct gref_iter_s *iter)
{
	debug("iter_fetch called, check rem_len(%u)", iter->stack->rem_len);
	while(iter->stack->rem_len == 0) {
		/* fetch link */
		if((iter->stack = gref_iter_fetch_link(iter, iter->stack)) == NULL) {
			return(0);
		}
	}

	return(gref_iter_fetch_base(iter->stack));
}

/**
 * @fn gref_iter_init
 */
gref_iter_t *gref_iter_init(
	gref_acv_t const *acv,
	gref_iter_params_t const *params)
{
	struct gref_s const *gref = (struct gref_s const *)acv;

	debug("init_stack_size(%lld)", gref->iter_init_stack_size);	
	if(gref == NULL || gref->type == GREF_POOL) { return(NULL); }

	/* restore params */
	static struct gref_iter_params_s const default_params = {
		.step_size = 1,
		.seq_direction = GREF_FW_ONLY
	};
	params = (params == NULL) ? &default_params : params;

	/* malloc mem */
	struct gref_iter_s *iter = (struct gref_iter_s *)lmm_malloc(acv->lmm,
		sizeof(struct gref_iter_s) + sizeof(uint64_t) * gref->iter_init_stack_size);
	if(iter == NULL) {
		return(NULL);
	}
	iter->lmm = acv->lmm;

	/* init param container */
	memset(iter->mem_arr, 0, sizeof(void *) * GREF_ITER_INTL_MEM_ARR_LEN);

	/* iterate from section 0 */
	iter->base_gid = 0;
	iter->tail_gid = _encode_id(acv->sec_cnt, 0);
	iter->step_gid = (params->seq_direction == GREF_FW_RV) ? 1 : 2;
	
	/* set params */
	iter->seed_len = acv->params.k;
	// iter->shift_len = 2 * (acv->params.k - 1);
	iter->seq_lim = gref->seq_lim;
	iter->link_table = acv->link_table;
	iter->hsec = (struct gref_section_half_s const *)hmap_get_object(acv->hmap, 0);

	/* create stack */
	// struct gref_iter_stack_s *stack = gref_iter_prepare_stack(iter);
	struct gref_iter_stack_s *stack = gref_iter_init_stack(iter);

	/* init kmer array */
	gref_iter_kmer_init(&stack->kmer, acv->params.k);
	return((gref_iter_t *)iter);
}

/**
 * @fn gref_iter_next
 */
struct gref_kmer_tuple_s gref_iter_next(
	gref_iter_t *_iter)
{
	struct gref_iter_s *iter = (struct gref_iter_s *)_iter;
	// struct gref_iter_stack_s *stack = iter->stack;
	// struct gref_iter_kmer_s *kmer = stack->kmer;

	while(!gref_iter_kmer_ready(&iter->stack->kmer)) {
		/* no kmer is available in the kmer array */
		uint8_t c = gref_iter_fetch(iter);
		debug("called iter_fetch(%u)", c);

		if(iter->stack != NULL) {
			gref_iter_kmer_append(&iter->stack->kmer, iter->stack->conv_table, c);
			continue;
		}

		debug("fetch section, gid(%u)", iter->base_gid + iter->step_gid);

		/* fetch the next section */
		if((iter->base_gid += iter->step_gid) < iter->tail_gid
		&& (iter->stack = gref_iter_init_stack(iter)) != NULL) {
			gref_iter_kmer_flush(&iter->stack->kmer);
			continue;
		}

		/* reached the end of the graph */
		debug("reached end");
		return((struct gref_kmer_tuple_s){
			.kmer = GREF_ITER_KMER_TERM,
			.gid_pos = (struct gref_gid_pos_s){
				.pos = 0,
				.gid = (uint32_t)-1
			}
		});
	}

	uint64_t k = gref_iter_kmer_next(&iter->stack->kmer);
	debug("return kmer(%llx), gid(%u), pos(%u)", k, iter->stack->sec_gid, iter->stack->len - iter->stack->rem_len);
	return((struct gref_kmer_tuple_s){
		.kmer = k,
		.gid_pos = (struct gref_gid_pos_s){
			.pos = iter->stack->len - iter->stack->rem_len,
			.gid = iter->base_gid
		}
	});
}

/**
 * @fn gref_iter_clean
 */
void gref_iter_clean(
	gref_iter_t *_iter)
{
	struct gref_iter_s *iter = (struct gref_iter_s *)_iter;
	debug("iter(%p)", iter);

	if(iter != NULL) {
		debug("stack(%p), iter(%p)", iter->stack, iter);
		for(int64_t i = 0; i < GREF_ITER_INTL_MEM_ARR_LEN; i++) {
			lmm_free(iter->lmm, iter->mem_arr[i]); iter->mem_arr[i] = NULL;
		}
		lmm_free(iter->lmm, (void *)iter);
	}
	return;
}


/* build kmer index (acv -> idx conversion) */
/**
 * @fn gref_build_kmer_idx_table
 */
static _force_inline
int64_t *gref_build_kmer_idx_table(
	struct gref_s *acv,
	struct gref_kmer_tuple_s *arr,
	int64_t size)
{
	lmm_kvec_t(int64_t) kmer_idx;
	lmm_kv_init(acv->lmm, kmer_idx);

	/* may fail when main memory is small */
	uint64_t kmer_idx_size = 0x01 << (2 * acv->params.k);
	lmm_kv_reserve(acv->lmm, kmer_idx, kmer_idx_size + 1);
	debug("ptr(%p), size(%llu)", lmm_kv_ptr(kmer_idx), kmer_idx_size);
	if(lmm_kv_ptr(kmer_idx) == NULL) { return(NULL); }

	uint64_t prev_kmer = 0;
	lmm_kv_push(acv->lmm, kmer_idx, prev_kmer);
	for(int64_t i = 0; i < size; i++) {
		uint64_t kmer = arr[i].kmer;
		debug("i(%lld), kmer(%llx), id(%u), pos(%u), prev_kmer(%llx)",
			i, kmer, arr[i].gid_pos.gid, arr[i].gid_pos.pos, prev_kmer);

		if(prev_kmer == kmer) { continue; }

		debug("index table for kmer(%llx) ends at %lld, next_kmer(%llx)", prev_kmer, i, kmer);

		/* sequence update detected */
		for(uint64_t j = prev_kmer + 1; j < kmer + 1; j++) {
			lmm_kv_push(acv->lmm, kmer_idx, i);
		}
		prev_kmer = kmer;
	}
	for(uint64_t j = prev_kmer; j < kmer_idx_size; j++) {
		lmm_kv_push(acv->lmm, kmer_idx, size);
	}
	return(lmm_kv_ptr(kmer_idx));
}

/**
 * @fn gref_shrink_kmer_table
 */
static _force_inline
struct gref_gid_pos_s *gref_shrink_kmer_table(
	struct gref_s *acv,
	struct gref_kmer_tuple_s *kmer_table,
	int64_t kmer_table_size)
{
	struct gref_gid_pos_s *packed_pos = (struct gref_gid_pos_s *)kmer_table;
	
	for(int64_t i = 0; i < kmer_table_size; i++) {
		packed_pos[i] = kmer_table[i].gid_pos;
	}

	return(lmm_realloc(acv->lmm, kmer_table, sizeof(struct gref_gid_pos_s) * kmer_table_size));
}

/**
 * @fn gref_build_index
 */
gref_idx_t *gref_build_index(
	gref_acv_t *acv)
{
	struct gref_s *gref = (struct gref_s *)acv;

	if(gref == NULL || gref->type != GREF_ACV) {
		goto _gref_build_index_error_handler;
	}

	/* enumerate kmers and pack into vector */
	lmm_kvec_t(struct gref_kmer_tuple_s) v;
	lmm_kv_init(acv->lmm, v);

	static struct gref_iter_params_s const iter_params = {
		.step_size = 1,
		.seq_direction = GREF_FW_RV
	};
	struct gref_iter_s *iter = gref_iter_init(acv, &iter_params);
	struct gref_kmer_tuple_s t;
	int64_t i = 0;
	while((t = gref_iter_next(iter)).gid_pos.gid != (uint32_t)-1) {
		i++;
		lmm_kv_push(acv->lmm, v, t);
	}
	gref_iter_clean(iter);

	/* sort kmers */
	if(psort_half(lmm_kv_ptr(v), lmm_kv_size(v),
		sizeof(struct gref_kmer_tuple_s), acv->params.num_threads) != 0) {
		lmm_kv_destroy(acv->lmm, v);
		debug("sort failed");
		goto _gref_build_index_error_handler;
	}

	/* build index of kmer table */
	gref->kmer_idx_table = gref_build_kmer_idx_table(acv, lmm_kv_ptr(v), lmm_kv_size(v));
	if(gref->kmer_idx_table == NULL) {
		debug("failed to build index table");
		goto _gref_build_index_error_handler;
	}

	/* shrink table */
	gref->kmer_table_size = lmm_kv_size(v);
	gref->kmer_table = gref_shrink_kmer_table(acv, lmm_kv_ptr(v), lmm_kv_size(v));
	if(gref->kmer_table == NULL) {
		debug("failed to shrink");
		goto _gref_build_index_error_handler;
	}

	/* store misc constants for kmer matching */
	gref->mask = (0x01<<(2 * gref->params.k)) - 1;

	/* change state */
	gref->type = GREF_IDX;
	return((gref_idx_t *)gref);

_gref_build_index_error_handler:;
	gref_clean((gref_t *)gref);
	return(NULL);
}

/**
 * @fn gref_disable_index
 */
gref_acv_t *gref_disable_index(
	gref_idx_t *idx)
{
	struct gref_s *gref = (struct gref_s *)idx;

	if(gref == NULL || gref->type != GREF_IDX) {
		gref_clean((gref_t *)gref);
		return(NULL);
	}

	/* cleanup kmer_idx_table */
	lmm_free(idx->lmm, idx->kmer_idx_table); idx->kmer_idx_table = NULL;

	/* change state */
	gref->type = GREF_ACV;
	return((gref_acv_t *)gref);
}

/**
 * @fn gref_match_2bitpacked
 */
struct gref_match_res_s gref_match_2bitpacked(
	gref_idx_t const *_gref,
	uint64_t seq)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	seq &= gref->mask;
	int64_t base = gref->kmer_idx_table[seq];
	int64_t tail = gref->kmer_idx_table[seq + 1];

	debug("seq(%llx), mask(%llx), base(%lld), tail(%lld)",
		seq, gref->mask, base, tail);
	return((struct gref_match_res_s){
		.gid_pos_arr = &gref->kmer_table[base],
		.len = tail - base
	});
}

/**
 * @fn gref_match
 * @brief seq length must be equal to k.
 */
struct gref_match_res_s gref_match(
	gref_idx_t const *_gref,
	uint8_t const *seq)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	int64_t const seed_len = gref->params.k;
	int64_t const shift_len = 2 * (seed_len - 1);

	uint64_t packed_seq = 0;
	for(int64_t i = 0; i < seed_len; i++) {
		packed_seq = (packed_seq>>2) | (gref_encode_2bit(seq[i])<<shift_len);
	}
	return(gref_match_2bitpacked((gref_t const *)gref, packed_seq));
}

/* misc */
#if 0
/**
 * @fn gref_dump_index
 */
int gref_dump_index(
	gref_t const *gref,
	zf_t *outfp)
{
	return(0);
}

/**
 * @fn gref_load_index
 */
gref_t *gref_load_index(
	zf_t *infp)
{
	return(NULL);
}
#endif

/**
 * @fn gref_get_section_count
 */
int64_t gref_get_section_count(
	gref_t const *_gref)
{
	struct gref_s *gref = (struct gref_s *)_gref;
	return((int64_t)gref->sec_cnt);
}

/**
 * @fn gref_get_section
 * @brief type must be ACV or IDX, otherwise return value is invalid
 */
struct gref_section_s const *gref_get_section(
	gref_acv_t const *_gref,
	uint32_t gid)
{
	struct gref_s *gref = (struct gref_s *)_gref;

	struct gref_section_half_s *base =
		(struct gref_section_half_s *)hmap_get_object(gref->hmap, 0);
	return((struct gref_section_s const *)&base[gid].sec);
}

/**
 * @fn gref_get_link
 * @brief type must be ACV or IDX, otherwise return value is invalid
 */
struct gref_link_s gref_get_link(
	gref_acv_t const *_gref,
	uint32_t gid)
{
	struct gref_s *gref = (struct gref_s *)_gref;

	struct gref_section_half_s *base =
		(struct gref_section_half_s *)hmap_get_object(gref->hmap, 0);
	return((struct gref_link_s){
		.gid_arr = &gref->link_table[base[gid].link_idx_base],
		.len = base[gid + 1].link_idx_base - base[gid].link_idx_base
	});
}

/**
 * @fn gref_get_name
 */
struct gref_str_s gref_get_name(
	gref_t const *_gref,
	uint32_t gid)
{
	struct gref_s *gref = (struct gref_s *)_gref;
	struct hmap_key_s key = hmap_get_key(gref->hmap, _decode_id(gid));
	return((struct gref_str_s){
		.ptr = key.ptr,
		.len = key.len
	});
}

#if 0
/* deprecated */
/**
 * @fn gref_get_ptr
 */
uint8_t const *gref_get_ptr(
	gref_t const *_gref)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	return((uint8_t const *)lmm_kv_ptr(gref->seq) + gref->params.seq_head_margin);
}
#endif

/**
 * @fn gref_get_total_len
 */
int64_t gref_get_total_len(
	gref_t const *_gref)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	return(gref->seq_len);
}

/**
 * @fn gref_get_lim
 * @brief type must be ACV or IDX, otherwise return value is invalid
 */
uint8_t const *gref_get_lim(
	gref_acv_t const *_gref)
{
	struct gref_s const *gref = (struct gref_s const *)_gref;
	return(gref->seq_lim);
}

/**
 * unittests
 */
unittest_config(
	.name = "gref",
	.depends_on = { "psort", "hmap", "zf" }
);

/**
 * @fn unittest_random_base
 */
static _force_inline
char unittest_random_base(void)
{
	char const table[4] = {0x01, 0x02, 0x04, 0x08};
	return(table[rand() % 4]);
}

/**
 * @fn unittest_generate_random_sequence
 */
#define UNITTEST_SEQ_MARGIN			( 8 )
static _force_inline
char *unittest_generate_random_sequence(
	int64_t len)
{
	char *seq;		/** a pointer to sequence */
	seq = (char *)malloc(sizeof(char) * (len + UNITTEST_SEQ_MARGIN));

	if(seq == NULL) { return NULL; }
	for(int64_t i = 0; i < len; i++) {
		seq[i] = unittest_random_base();
	}
	seq[len] = '\0';
	return seq;
}


#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

#define _pack(x) ({ \
	uint64_t _packed_seq = 0; \
	int64_t _len = strlen(x); \
	uint8_t _shift_len = 2 * (_len - 1); \
	for(int64_t i = 0; i < _len; i++) { \
		_packed_seq = (_packed_seq>>2) | (gref_encode_2bit((x)[i])<<_shift_len); \
	} \
	_packed_seq; \
})

/* make pool context */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4));

	assert(pool != NULL);

	gref_clean(pool);
}

/* add segment */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	int ret = gref_append_segment(pool, _str("sec0"), _seq("AARA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_segment(pool, _str("sec1"), _seq("MAAA"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_segment(pool, _str("sec2"), _seq("ACGT"));
	assert(ret == 0, "ret(%d)", ret);

	ret = gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	assert(ret == 0, "ret(%d)", ret);

	/* section count */
	assert(gref_get_section_count(pool) == 3, "%lld", gref_get_section_count(pool));

	/* total len */
	assert(gref_get_total_len(pool) == 12, "len(%lld)", gref_get_total_len(pool));

	gref_clean(pool);
}

/* archive (fw_copy) */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_direction = GREF_FW_ONLY,
		.copy_mode = GREF_COPY,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_count(acv) == 3, "%lld", gref_get_section_count(acv));

	/* total len */
	assert(gref_get_total_len(acv) == 13, "len(%lld)", gref_get_total_len(acv));

	/* get lim */
	assert(gref_get_lim(acv) == GREF_SEQ_LIM, "lim(%p)", gref_get_lim(acv));

	/* check bases at the head of the sections */
	assert(gref_get_section(acv, 0)->base[0] == 0x04, "%x", gref_get_section(acv, 0)->base[0]);
	assert(gref_get_section(acv, 2)->base[0] == 0x03, "%x", gref_get_section(acv, 2)->base[0]);
	assert(gref_get_section(acv, 4)->base[0] == 0x01, "%x", gref_get_section(acv, 4)->base[0]);

	gref_clean(acv);
}

/* archive (fw_nocopy) */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_direction = GREF_FW_ONLY,
		.seq_format = GREF_4BIT,
		.copy_mode = GREF_NOCOPY,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("\x04\x04\x05\x01"));
	gref_append_segment(pool, _str("sec1"), _seq("\x03"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("\x01\x03\x07\x07\x04\x08\x04\x08"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_count(acv) == 3, "%lld", gref_get_section_count(acv));

	/* total len */
	assert(gref_get_total_len(acv) == 13, "len(%lld)", gref_get_total_len(acv));

	/* get lim */
	assert(gref_get_lim(acv) == GREF_SEQ_LIM, "lim(%p)", gref_get_lim(acv));

	/* check bases at the head of the sections */
	assert(gref_get_section(acv, 0)->base[0] == 0x04, "%x", gref_get_section(acv, 0)->base[0]);
	assert(gref_get_section(acv, 2)->base[0] == 0x03, "%x", gref_get_section(acv, 2)->base[0]);
	assert(gref_get_section(acv, 4)->base[0] == 0x01, "%x", gref_get_section(acv, 4)->base[0]);

	gref_clean(acv);
}

/* archive (fr_copy) */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_direction = GREF_FW_RV,
		.copy_mode = GREF_COPY,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_count(acv) == 3, "%lld", gref_get_section_count(acv));

	/* total len */
	assert(gref_get_total_len(acv) == 13, "len(%lld)", gref_get_total_len(acv));

	/* get lim */
	assert(gref_get_lim(acv) == lmm_kv_ptr(acv->seq) + 2 * acv->seq_len + 32, "lim(%p)", gref_get_lim(acv));

	/* check bases at the head of the sections */
	assert(gref_get_section(acv, 0)->base[0] == 0x04, "%x", gref_get_section(acv, 0)->base[0]);
	assert(gref_get_section(acv, 1)->base[0] == 0x08, "%x", gref_get_section(acv, 1)->base[0]);
	assert(gref_get_section(acv, 2)->base[0] == 0x03, "%x", gref_get_section(acv, 2)->base[0]);
	assert(gref_get_section(acv, 3)->base[0] == 0x0c, "%x", gref_get_section(acv, 3)->base[0]);
	assert(gref_get_section(acv, 4)->base[0] == 0x01, "%x", gref_get_section(acv, 4)->base[0]);
	assert(gref_get_section(acv, 5)->base[0] == 0x01, "%x", gref_get_section(acv, 5)->base[0]);

	gref_clean(acv);
}

/* archive (fr_nocopy) */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_direction = GREF_FW_RV,
		.seq_format = GREF_4BIT,
		.copy_mode = GREF_NOCOPY,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* append */
	#define _seqrv(x)		(uint8_t const *)(x), strlen(x)/2

	gref_append_segment(pool, _str("sec0"), _seqrv("\x04\x04\x05\x01\x08\x0a\x02\x02"));
	gref_append_segment(pool, _str("sec1"), _seqrv("\x03\x0c"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seqrv("\x01\x03\x07\x07\x04\x08\x04\x08\x01\x02\x01\x02\x0e\x0e\x0c\x08"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	#undef _seqrv

	/* build index */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_count(acv) == 3, "%lld", gref_get_section_count(acv));

	/* total len */
	assert(gref_get_total_len(acv) == 13, "len(%lld)", gref_get_total_len(acv));

	/* get lim */
	assert(gref_get_lim(acv) == GREF_SEQ_LIM, "lim(%p)", gref_get_lim(acv));

	/* check bases at the head of the sections */
	assert(gref_get_section(acv, 0)->base[0] == 0x04, "%x", gref_get_section(acv, 0)->base[0]);
	assert(gref_get_section(acv, 1)->base[0] == 0x08, "%x", gref_get_section(acv, 1)->base[0]);
	assert(gref_get_section(acv, 2)->base[0] == 0x03, "%x", gref_get_section(acv, 2)->base[0]);
	assert(gref_get_section(acv, 3)->base[0] == 0x0c, "%x", gref_get_section(acv, 3)->base[0]);
	assert(gref_get_section(acv, 4)->base[0] == 0x01, "%x", gref_get_section(acv, 4)->base[0]);
	assert(gref_get_section(acv, 5)->base[0] == 0x01, "%x", gref_get_section(acv, 5)->base[0]);

	gref_clean(acv);
}

/* archive (fr_copy / large sequence) */
unittest()
{
	int64_t const len = 100000;
	int64_t const cnt = 1000;

	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_direction = GREF_FW_RV,
		.seq_format = GREF_4BIT,
		.copy_mode = GREF_COPY,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* dump */
	for(int64_t i = 0; i < cnt; i++) {

		/* name */
		char buf[1024];
		sprintf(buf, "seq%" PRId64 "", i);

		/* seq */
		char *seq = unittest_generate_random_sequence(len);

		/* append */
		gref_append_segment(pool, buf, strlen(buf), (uint8_t const *)seq, strlen(seq));
		free(seq);
	}

	/* freeze */
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	assert(acv->type == GREF_ACV, "%d", acv->type);
	assert(acv->link_table != NULL, "%p", acv->link_table);

	/* section count */
	assert(gref_get_section_count(acv) == cnt, "%lld", gref_get_section_count(acv));

	/* total len */
	assert(gref_get_total_len(acv) == cnt * len, "len(%lld)", gref_get_total_len(acv));

	gref_clean(acv);
}


/* seed iteration */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32,
		.seq_direction = GREF_FW_ONLY));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	
	/* build archive */
	gref_acv_t *acv = gref_freeze_pool(pool);

	/* create iterator */
	gref_iter_t *iter = gref_iter_init(acv, NULL);
	assert(iter != NULL, "%p", iter);

	/* enumerate */
	#define _f(_i)					gref_iter_next(_i)
	#define _check_kmer(_t, _k, _id, _pos) ( \
		   (_t).kmer == _pack(_k) \
		&& (_t).gid_pos.gid == _encode_id(_id, 0) \
		&& (_t).gid_pos.pos == (_pos) \
	)
	#define _print_kmer(_t) \
		"kmer(%llx), sec(%u), pos(%u)", \
		(_t).kmer, _decode_id((_t).gid_pos.gid), (_t).gid_pos.pos

	struct gref_kmer_tuple_s t;

	/* sec0 */
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGGA", 0, 0), _print_kmer(t));

	/* sec0-sec1 */
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAC", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAC", 0, 1), _print_kmer(t));

	/* sec0-sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AAAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GACA", 0, 2), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAAC", 0, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACAC", 0, 3), _print_kmer(t));

	/* sec0-sec2 */
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 1 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 1 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AAAC", 0, 2 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAC", 0, 2 + K_MAX), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AACA", 0, 3 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACC", 0, 3 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACG", 0, 3 + K_MAX), _print_kmer(t));

	/* sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AACA", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACA", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACC", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACC", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACG", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACG", 1, 0), _print_kmer(t));

	/* sec2 */
	t = _f(iter); assert(_check_kmer(t, "ACAA", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCA", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGA", 2, 0), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACAC", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCC", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGC", 2, 0), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACAG", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCG", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGG", 2, 0), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "CAAG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCAG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGAG", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CACG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCCG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGCG", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CAGG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCGG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGGG", 2, 1), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CAGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GCGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGTG", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGTG", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGTG", 2, 3), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "GTGT", 2, 4), _print_kmer(t));

	t = _f(iter); assert(t.kmer == GREF_ITER_KMER_TERM, "%llx, %llx", t.kmer, GREF_ITER_KMER_TERM);
	assert(t.gid_pos.gid == (uint32_t)-1, "gid(%u)", t.gid_pos.gid);

	gref_iter_clean(iter);
	gref_clean(acv);
}


/* N's in section */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32,
		.seq_direction = GREF_FW_ONLY));
	gref_append_segment(pool, _str("sec0"), _seq("GGRANNNNGTTCANNNNNAAAAT"));
	gref_append_segment(pool, _str("sec1"), _seq("TNNNCCCCC"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);

	/* build archive */
	gref_acv_t *acv = gref_freeze_pool(pool);

	/* create iterator */
	gref_iter_t *iter = gref_iter_init(acv, NULL);
	assert(iter != NULL, "%p", iter);
	struct gref_kmer_tuple_s t;

	/* sec0 */
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGGA", 0, 0), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "GTTC", 0, 8), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "TTCA", 0, 9), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAAA", 0, 18), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AAAT", 0, 19), _print_kmer(t));

	/* sec0-sec1 */
	t = _f(iter); assert(_check_kmer(t, "AATT", 0, 20), _print_kmer(t));

	/* sec1 */
	t = _f(iter); assert(_check_kmer(t, "CCCC", 1, 4), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCCC", 1, 5), _print_kmer(t));

	t = _f(iter); assert(t.kmer == GREF_ITER_KMER_TERM, "%llx, %llx", t.kmer, GREF_ITER_KMER_TERM);
	assert(t.gid_pos.gid == (uint32_t)-1, "gid(%u)", t.gid_pos.gid);

	#undef _f
	#undef _check_kmer
	#undef _print_kmer

	gref_iter_clean(iter);
	gref_clean(acv);
}

/* build index */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));

	/* append */
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);

	/* build index */	
	gref_acv_t *acv = gref_freeze_pool(pool);
	assert(acv != NULL, "acv(%p)", acv);

	gref_idx_t *idx = gref_build_index(pool);
	assert(idx != NULL, "idx(%p)", idx);

	/* section count */
	assert(gref_get_section_count(idx) == 3, "%lld", gref_get_section_count(idx));

	/* total len */
	assert(gref_get_total_len(idx) == 16, "len(%lld)", gref_get_total_len(idx));

	gref_clean(idx);
}

/* build iterator from gref_idx_t */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("M"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	
	/* build archive */
	gref_acv_t *idx = gref_build_index(gref_freeze_pool(pool));

	/* create iterator */
	gref_iter_t *iter = gref_iter_init(idx, NULL);
	assert(iter != NULL, "%p", iter);

	/* enumerate */
	#define _f(_i)					gref_iter_next(_i)
	#define _check_kmer(_t, _k, _id, _pos) ( \
		   (_t).kmer == _pack(_k) \
		&& (_t).gid_pos.gid == _encode_id(_id, 0) \
		&& (_t).gid_pos.pos == (_pos) \
	)
	#define _print_kmer(_t) \
		"kmer(%llx), sec(%u), pos(%u)", \
		(_t).kmer, _decode_id((_t).gid_pos.gid), (_t).gid_pos.pos

	struct gref_kmer_tuple_s t;

	/* sec0 */
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGGA", 0, 0), _print_kmer(t));

	/* sec0-sec1 */
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAC", 0, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAC", 0, 1), _print_kmer(t));

	/* sec0-sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AAAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACA", 0, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GACA", 0, 2), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAAC", 0, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACAC", 0, 3), _print_kmer(t));

	/* sec0-sec2 */
	t = _f(iter); assert(_check_kmer(t, "GAAA", 0, 1 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGAA", 0, 1 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AAAC", 0, 2 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAAC", 0, 2 + K_MAX), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AACA", 0, 3 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACC", 0, 3 + K_MAX), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACG", 0, 3 + K_MAX), _print_kmer(t));

	/* sec1-sec2 */
	t = _f(iter); assert(_check_kmer(t, "AACA", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACA", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACC", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACC", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "AACG", 1, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CACG", 1, 0), _print_kmer(t));

	/* sec2 */
	t = _f(iter); assert(_check_kmer(t, "ACAA", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCA", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGA", 2, 0), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACAC", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCC", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGC", 2, 0), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACAG", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACCG", 2, 0), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "ACGG", 2, 0), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "CAAG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCAG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGAG", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CACG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCCG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGCG", 2, 1), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "CAGG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCGG", 2, 1), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGGG", 2, 1), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "AAGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CAGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GAGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "ACGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CCGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GCGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGGT", 2, 2), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGGT", 2, 2), _print_kmer(t));
	
	t = _f(iter); assert(_check_kmer(t, "AGTG", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "CGTG", 2, 3), _print_kmer(t));
	t = _f(iter); assert(_check_kmer(t, "GGTG", 2, 3), _print_kmer(t));

	t = _f(iter); assert(_check_kmer(t, "GTGT", 2, 4), _print_kmer(t));

	t = _f(iter); assert(t.kmer == GREF_ITER_KMER_TERM, "%llx, %llx", t.kmer, GREF_ITER_KMER_TERM);
	assert(t.gid_pos.gid == (uint32_t)-1, "gid(%u)", t.gid_pos.gid);

	#undef _f
	#undef _check_kmer
	#undef _print_kmer

	gref_iter_clean(iter);
	gref_clean(idx);
}

/* get_section */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));

	#define _calc_base_sum(_id) ( \
		  (uint64_t)gref_get_section(idx, (_id))->base \
		+ (uint64_t)gref_get_section(idx, (_id) + 1)->base \
		+ (uint64_t)gref_get_section(idx, (_id) + 1)->len \
	)

	/* section id is given in ascending order from 0 */
	assert(gref_get_section(idx, 0) != NULL, "%p", gref_get_section(idx, 0));
	assert(gref_get_section(idx, 0)->gid == 0, "gid(%u)", gref_get_section(idx, 0)->gid);
	assert(gref_get_section(idx, 0)->len == 4, "len(%u)", gref_get_section(idx, 0)->len);

	/* section 0 in reverse */
	assert(gref_get_section(idx, 1) != NULL, "%p", gref_get_section(idx, 1));
	assert(gref_get_section(idx, 1)->gid == 1, "gid(%u)", gref_get_section(idx, 1)->gid);
	assert(gref_get_section(idx, 1)->len == 4, "len(%u)", gref_get_section(idx, 1)->len);

	/* check section 0/1 ptr */
	assert(_calc_base_sum(0) == 2 * (uint64_t)GREF_SEQ_LIM, "%llx", _calc_base_sum(0));

	/* section 1 */
	assert(gref_get_section(idx, 2) != NULL, "%p", gref_get_section(idx, 2));
	assert(gref_get_section(idx, 2)->gid == 2, "gid(%u)", gref_get_section(idx, 2)->gid);
	assert(gref_get_section(idx, 2)->len == 4, "len(%u)", gref_get_section(idx, 2)->len);

	/* section 1 in reverse */
	assert(gref_get_section(idx, 3) != NULL, "%p", gref_get_section(idx, 3));
	assert(gref_get_section(idx, 3)->gid == 3, "gid(%u)", gref_get_section(idx, 3)->gid);
	assert(gref_get_section(idx, 3)->len == 4, "len(%u)", gref_get_section(idx, 3)->len);

	/* check section 2/3 ptr */
	assert(_calc_base_sum(2) == 2 * (uint64_t)GREF_SEQ_LIM, "%llx", _calc_base_sum(2));

	/* section 2 */
	assert(gref_get_section(idx, 4) != NULL, "%p", gref_get_section(idx, 4));
	assert(gref_get_section(idx, 4)->gid == 4, "gid(%u)", gref_get_section(idx, 4)->gid);
	assert(gref_get_section(idx, 4)->len == 8, "len(%u)", gref_get_section(idx, 4)->len);

	/* section 2 in reverse */
	assert(gref_get_section(idx, 5) != NULL, "%p", gref_get_section(idx, 5));
	assert(gref_get_section(idx, 5)->gid == 5, "gid(%u)", gref_get_section(idx, 5)->gid);
	assert(gref_get_section(idx, 5)->len == 8, "len(%u)", gref_get_section(idx, 5)->len);

	/* check section 4/5 ptr */
	assert(_calc_base_sum(4) == 2 * (uint64_t)GREF_SEQ_LIM, "%llx", _calc_base_sum(4));

	#undef _calc_base_sum

	gref_clean(idx);
}

/* get_link */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));

	struct gref_link_s l = gref_get_link(idx, 0);
	assert(l.len == 2, "len(%lld)", l.len);
	assert(l.gid_arr[0] == 2, "%u", l.gid_arr[0]);
	assert(l.gid_arr[1] == 4, "%u", l.gid_arr[1]);

	l = gref_get_link(idx, 1);
	assert(l.len == 0, "len(%lld)", l.len);

	l = gref_get_link(idx, 2);
	assert(l.len == 1, "len(%lld)", l.len);
	assert(l.gid_arr[0] == 4, "%u", l.gid_arr[0]);

	l = gref_get_link(idx, 3);
	assert(l.len == 1, "len(%lld)", l.len);
	assert(l.gid_arr[0] == 1, "%u", l.gid_arr[0]);

	l = gref_get_link(idx, 4);
	assert(l.len == 0, "len(%lld)", l.len);

	l = gref_get_link(idx, 5);
	assert(l.len == 2, "len(%lld)", l.len);
	assert(l.gid_arr[0] == 3, "%u", l.gid_arr[0]);
	assert(l.gid_arr[1] == 1, "%u", l.gid_arr[1]);

	gref_clean(idx);
}

/* get_name */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));

	/* section id is given in ascending order from 0 */
	assert(gref_get_name(idx, 0).len == 4, "%d", gref_get_name(idx, 0).len);
	assert(strcmp(gref_get_name(idx, 0).ptr, "sec0") == 0, "%s", gref_get_name(idx, 0).ptr);

	assert(gref_get_name(idx, 1).len == 4, "%d", gref_get_name(idx, 1).len);
	assert(strcmp(gref_get_name(idx, 1).ptr, "sec0") == 0, "%s", gref_get_name(idx, 1).ptr);

	/* section 1 */
	assert(gref_get_name(idx, 2).len == 4, "%d", gref_get_name(idx, 2).len);
	assert(strcmp(gref_get_name(idx, 2).ptr, "sec1") == 0, "%s", gref_get_name(idx, 2).ptr);

	assert(gref_get_name(idx, 3).len == 4, "%d", gref_get_name(idx, 3).len);
	assert(strcmp(gref_get_name(idx, 3).ptr, "sec1") == 0, "%s", gref_get_name(idx, 3).ptr);

	/* section 2 */
	assert(gref_get_name(idx, 4).len == 4, "%d", gref_get_name(idx, 4).len);
	assert(strcmp(gref_get_name(idx, 4).ptr, "sec2") == 0, "%s", gref_get_name(idx, 4).ptr);

	assert(gref_get_name(idx, 5).len == 4, "%d", gref_get_name(idx, 5).len);
	assert(strcmp(gref_get_name(idx, 5).ptr, "sec2") == 0, "%s", gref_get_name(idx, 5).ptr);

	gref_clean(idx);
}

/* match */
unittest()
{
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = 4,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	gref_append_segment(pool, _str("sec0"), _seq("GGRA"));
	gref_append_segment(pool, _str("sec1"), _seq("MGGG"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec1"), 0);
	gref_append_link(pool, _str("sec1"), 0, _str("sec2"), 0);
	gref_append_segment(pool, _str("sec2"), _seq("ACVVGTGT"));
	gref_append_link(pool, _str("sec0"), 0, _str("sec2"), 0);
	gref_idx_t *idx = gref_build_index(gref_freeze_pool(pool));


	/* without ambiguous bases */
	struct gref_match_res_s r = gref_match(idx, (uint8_t const *)"GTGT");
	assert(r.gid_pos_arr != NULL, "%p", r.gid_pos_arr);
	assert(r.len == 2, "%lld", r.len);

	/* check pos */
	assert(r.gid_pos_arr[0].pos == 4, "%u", r.gid_pos_arr[0].pos);

	/* check section */
	struct gref_section_s const *sec = gref_get_section(idx, r.gid_pos_arr[0].gid);
	assert(sec->gid == 4, "gid(%u)", sec->gid);
	assert(sec->len == 8, "len(%u)", sec->len);

	/* check pos */
	assert(r.gid_pos_arr[1].pos == 4, "%u", r.gid_pos_arr[1].pos);

	/* check section */
	sec = gref_get_section(idx, r.gid_pos_arr[1].gid);
	assert(sec->gid == 5, "gid(%u)", sec->gid);
	assert(sec->len == 8, "len(%u)", sec->len);


	/* with ambiguous bases */
	r = gref_match(idx, (uint8_t const *)"CGGG");
	assert(r.gid_pos_arr != NULL, "%p", r.gid_pos_arr);
	assert(r.len == 3, "%lld", r.len);

	/* check pos */
	assert(r.gid_pos_arr[0].pos == 0, "%u", r.gid_pos_arr[0].pos);

	/* check section */
	sec = gref_get_section(idx, r.gid_pos_arr[0].gid);
	assert(sec->gid == 2, "gid(%u)", sec->gid);
	assert(sec->len == 4, "len(%u)", sec->len);

	/* check pos */
	assert(r.gid_pos_arr[1].pos == 1, "%u", r.gid_pos_arr[1].pos);

	/* check section */
	sec = gref_get_section(idx, r.gid_pos_arr[1].gid);
	assert(sec->gid == 4, "gid(%u)", sec->gid);
	assert(sec->len == 8, "len(%u)", sec->len);

	/* check pos */
	assert(r.gid_pos_arr[2].pos == 3, "%u", r.gid_pos_arr[2].pos);

	/* check section */
	sec = gref_get_section(idx, r.gid_pos_arr[2].gid);
	assert(sec->gid == 5, "gid(%u)", sec->gid);
	assert(sec->len == 8, "len(%u)", sec->len);

	gref_clean(idx);
}

/**
 * end of gref.c
 */
