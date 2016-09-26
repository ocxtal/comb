
/**
 * @file sr.h
 *
 * @brief sequence reader implementation
 */
#ifndef _SR_H_INCLUDED
#define _SR_H_INCLUDED

#include <stdint.h>
#include "gref.h"

/**
 * @type sr_t
 */
typedef struct sr_s sr_t;

/**
 * @enum sr_format
 * @brief format flag constant
 */
enum sr_format {
	SR_UNKNOWN		= 0,
	SR_FASTA		= 1,
	SR_FASTQ		= 2,
	SR_FAST5		= 3,
	SR_GFA			= 4
};

/**
 * @enum sr_revcomp
 */
enum sr_revcomp {
	SR_FW_ONLY		= 1,
	SR_FW_RV		= 2
};

/**
 * @struct sr_params_s
 */
struct sr_params_s {
	uint8_t k;					/* kmer length */
	uint8_t seq_direction;		/* FW_ONLY or FW_RV */
	uint8_t format;				/* equal to fna_params_t.file_format */
	uint8_t reserved1;
	uint16_t num_threads;
	uint16_t reserved2;
	uint32_t pool_size;
	uint32_t read_mem_size;
	void *lmm;					/* lmm memory manager */
};
typedef struct sr_params_s sr_params_t;

#define SR_PARAMS(...)		( &((struct sr_params_s const){ __VA_ARGS__ }) )

/**
 * @struct sr_gref_s
 * @brief gref and iter container
 */
struct sr_gref_s {
	void *lmm;
	char const *path;
	gref_t const *gref;
	gref_iter_t *iter;
	void *reserved1[2];
	uint32_t reserved2[2];
};

/**
 * @fn sr_init
 */
sr_t *sr_init(
	char const *path,
	sr_params_t const *params);

/**
 * @fn sr_clean
 */
void sr_clean(
	sr_t *sr);

/**
 * @fn sr_get_index
 */
struct sr_gref_s *sr_get_index(
	sr_t *sr);

/**
 * @fn sr_get_iter
 */
struct sr_gref_s *sr_get_iter(
	sr_t *sr);

/**
 * @fn sr_gref_free
 */
void sr_gref_free(
	struct sr_gref_s *gref);

#endif /* _SR_H_INCLUDED */
/**
 * end of sr.h
 */
