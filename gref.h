
/**
 * @file gref.h
 *
 * @brief a header of gref.c
 *
 * @author Hajime Suzuki
 * @date 2015/6/23
 * @license Apache v2.
 *
 * @detail
 * Ref is a hash-based sequence indexer and exact matcher. 
 */

#ifndef _GREF_H_INCLUDED
#define _GREF_H_INCLUDED

#include <stdint.h>

/**
 * @enum gref_error
 * @brief error flags
 */
enum gref_error {
	/** error codes */
	GREF_SUCCESS 				= 0,
	GREF_ERROR 					= 1,
	GREF_ERROR_INVALID_CONTEXT	= 2,
	GREF_ERROR_INVALID_ARGS 	= 3,
	GREF_ERROR_OVERWRITE		= 4,
	GREF_ERROR_FILE_NOT_FOUND	= 5,
	GREF_ERROR_BROKEN_FILE		= 6,

	/** return values */
	GREF_INDEX_VALID			= 0,
	GREF_INDEX_INVALID			= -1
};

/**
 * @enum gref_revcomp
 */
enum gref_revcomp {
	GREF_FW_ONLY				= 1,
	GREF_FW_RV					= 2
};

/**
 * @enum gref_format
 */
enum gref_format_flags {
	GREF_ASCII					= 1,
	GREF_4BIT					= 2,
};

/**
 * @enum gref_copy_mode
 *
 * @brief sequences passed to the gref object must remain in the same
 * location until the object is destoryed if GREF_NOCOPY is specified.
 */
enum gref_copy_mode {
	GREF_COPY 					= 1,
	GREF_NOCOPY					= 2
};

/**
 * @type gref_t
 */
typedef struct gref_s gref_t;

/**
 * @type gref_pool_t
 * @brief mutable sequence pool
 */
typedef struct gref_s gref_pool_t;

/**
 * @type gref_acv_t
 * @brief immutable sequence pool, providing kmer iterator, converted from gref_pool_t
 */
typedef struct gref_s gref_acv_t;

/**
 * @type gref_idx_t
 * @brief immutable sequence pool with kmer index, converted from gref_acv_t
 */
typedef struct gref_s gref_idx_t;

/**
 * @type gref_iter_t
 */
typedef struct gref_iter_s gref_iter_t;

/**
 * @struct gref_params_s
 */
struct gref_params_s {
	uint8_t k;						/* kmer length */
	uint8_t seq_direction;
	uint8_t seq_format;
	uint8_t copy_mode;
	uint16_t num_threads;
	uint16_t reserved;
	uint32_t hash_size;
	uint16_t seq_head_margin;
	uint16_t seq_tail_margin;
	void *lmm;
};
typedef struct gref_params_s gref_params_t;
#define GREF_PARAMS(...)			( &((struct gref_params_s const) { __VA_ARGS__ }) )

/**
 * @struct gref_section_s
 * @brief has equivalent fields to struct sea_section_s
 */
struct gref_section_s {
	uint32_t gid;
	uint32_t len;
	uint8_t const *base;
};
typedef struct gref_section_s gref_section_t;

/**
 * @struct gref_link_s
 */
struct gref_link_s {
	uint32_t const *gid_arr;
	int64_t len;
};
typedef struct gref_link_s gref_link_t;

/**
 * @struct gref_str_s
 */
struct gref_str_s {
	char const *str;
	int32_t len;
};
typedef struct gref_str_s gref_str_t;

/**
 * @struct gref_gid_pos_s
 */
struct gref_gid_pos_s {
	uint32_t gid;
	uint32_t pos;
};
typedef struct gref_gid_pos_s gref_gid_pos_t;

/**
 * @struct gref_kmer_tuple_s
 */
struct gref_kmer_tuple_s {
	uint64_t kmer;
	struct gref_gid_pos_s gid_pos;
};
typedef struct gref_kmer_tuple_s gref_kmer_tuple_t;

/* encode and decode id */
#define GREF_FW 				( 0 )
#define GREF_RV 				( 0x01 )
#define gref_rev_gid(_gid)		( 0x01 ^ (_gid) )
#define gref_gid(_id, _d)		( ((_id)<<1) | (0x01 & (_d)) )
#define gref_id(_gid)			( (_gid)>>1 )
#define gref_dir(_gid)			( (_gid) & 0x01 )

/**
 * @struct gref_match_res_s
 */
struct gref_match_res_s {
	struct gref_gid_pos_s *gid_pos_arr;
	int64_t len;
};
typedef struct gref_match_res_s gref_match_res_t;

/**
 * @fn gref_init_pool
 * @brief initialize mutable reference object (reference index precursor)
 */
gref_pool_t *gref_init_pool(
	gref_params_t const *params);

/**
 * @fn gref_freeze_pool
 */
gref_acv_t *gref_freeze_pool(
	gref_pool_t *pool);

/**
 * @fn gref_melt_archive
 */
gref_pool_t *gref_melt_archive(
	gref_acv_t *acv);

/**
 * @fn gref_build_index
 * @brief build index
 */
gref_idx_t *gref_build_index(
	gref_acv_t *acv);

/**
 * @fn gref_disable_index
 */
gref_acv_t *gref_disable_index(
	gref_idx_t *idx);

/**
 * @fn gref_clean
 * @brief cleanup object. gref can be pool, acv, or idx.
 */
void gref_clean(
	gref_t *gref);

/**
 * @fn gref_append_segment
 * @brief append a sequence block to the context.
 */
int gref_append_segment(
	gref_pool_t *pool,
	char const *name,
	int32_t name_len,
	uint8_t const *seq,
	int64_t seq_len);

/**
 * @fn gref_append_link
 *
 * @brief append a edge on graph
 */
int gref_append_link(
	gref_pool_t *pool,
	char const *src,
	int32_t src_len,
	int32_t src_ori,
	char const *dst,
	int32_t dst_len,
	int32_t dst_ori);

/**
 * @fn gref_append_snp
 * @brief not implemented yet (;_;)
 */
int gref_append_snp(
	gref_pool_t *_pool,
	char const *name,
	int32_t name_len,
	int64_t pos,
	uint8_t snp);

/**
 * @fn gref_split_segment
 * @brief not implemented yet (;_;)
 */
int gref_split_segment(
	gref_pool_t *_pool,
	char const *base,
	int32_t base_len,
	int64_t pos,
	char const *splitted,
	int32_t splitted_len);

#if 0
/**
 * @fn gref_load_index
 */
gref_idx_t *gref_load_index(
	zf_t *fp);

/**
 * @fn gref_dump_index
 */
int gref_dump_index(
	gref_idx_t const *gref,
	zf_t *fp);
#endif

/**
 * @fn gref_iter_init, gref_iter_next, gref_iter_clean
 *
 * @brief kmer iterator
 */
gref_iter_t *gref_iter_init(
	gref_acv_t const *gref);

/**
 * @fn gref_iter_next
 */
#define GREF_ITER_KMER_TERM 		( (uint64_t)0xffffffffffffffff )
gref_kmer_tuple_t gref_iter_next(
	gref_iter_t *iter);

/**
 * @fn gref_iter_clean
 */
void gref_iter_clean(
	gref_iter_t *iter);

/**
 * @fn gref_match
 *
 * @brief index must be build with kmer-hash mode.
 */
struct gref_match_res_s gref_match(
	gref_idx_t const *gref,
	uint8_t const *seq);
struct gref_match_res_s gref_match_2bitpacked(
	gref_idx_t const *gref,
	uint64_t seq);

/**
 * @fn gref_get_section_count
 */
int64_t gref_get_section_count(
	gref_t const *gref);

/**
 * @fn gref_get_section
 */
struct gref_section_s const *gref_get_section(
	gref_acv_t const *gref,
	uint32_t gid);

/**
 * @fn gref_get_link
 */
struct gref_link_s gref_get_link(
	gref_t const *gref,
	uint32_t gid);

/**
 * @fn gref_get_name
 */
struct gref_str_s gref_get_name(
	gref_t const *gref,
	uint32_t gid);

#if 0
/* deprecated */
/**
 * @fn gref_get_ptr
 */
uint8_t const *gref_get_ptr(
	gref_t const *gref);
#endif

/**
 * @fn gref_get_total_len
 */
int64_t gref_get_total_len(
	gref_t const *gref);

/**
 * @fn gref_get_lim
 */
uint8_t const *gref_get_lim(
	gref_t const *gref);
#define gref_rev_ptr(ptr, lim)		( (uint8_t const *)(lim) + (uint64_t)(lim) - (uint64_t)(ptr) - 1 )

#if 0
/**
 * @fn gref_is_amb
 */
int64_t gref_is_amb(
	gref_t const *gref,
	int64_t lb, int64_t ub);
#endif

#endif /** #ifndef _GREF_H_INCLUDED */
/**
 * end of gref.h
 */
