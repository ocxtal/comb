
/**
 * @file fna.h
 *
 * @brief FASTA / FASTQ reader implementation.
 *
 * @author Hajime Suzuki
 * @date 2015/05/23
 * @license Apache v2.
 *
 * @detail
 * Supported formats:
 *   FASTA (raw and gzipped)
 *   FASTQ (raw and gzipped)
 *   FAST5 (unsupported for now!!!)
 *
 * List of APIs:
 *   Basic readers:
 *     fna_t *fna_init(char const *path, int pack);
 *     fna_seq_t *fna_read(fna_t const *fna, fna_seq_t *seq);
 *     void fna_seq_free(fna_seq_t *seq);
 *     void fna_close(fna_t *fna);
 *
 *   Sequence duplicators:
 *     fna_seq_t *fna_duplicate(fna_seq_t const *seq);
 *     fna_seq_t *fna_revcomp(fna_seq_t const *seq);
 *
 *   Sequence modifiers:
 *     void fna_append(fna_seq_t *dst, fna_seq_t const *src);
 *     void fna_append_revcomp(fna_seq_t *dst, fna_seq_t const *src);
 *
 * Types and members:
 *   fna_t (alias to struct _fna): sequence reader instance container.
 *     path: path to the file.
 *   fna_seq_t (alias to struct _seq): sequence container.
 *     name: sequence name container (kvec_t(char) instance)
 *       name.a: pointer to the sequence name (null-terminated ASCII).
 *     seq: sequence container (kvec_t(uint8_t) or kpvec_t(uint8_t) instance)
 *       seq.a: pointer to the sequence (null-terminated when fna->pack == FNA_RAW)
 */
#ifndef _FNA_H_INCLUDED
#define _FNA_H_INCLUDED

#include <stdint.h>

/**
 * @enum fna_flag_encode
 */
enum fna_flag_encode {
	FNA_RAW 		= 0,
	FNA_ASCII 		= 0,
	FNA_2BIT		= 1,
	FNA_2BITPACKED 	= 2,
	FNA_4BIT		= 3,
	FNA_4BITPACKED 	= 4,
};

/**
 * @enum fna_format
 *
 * @brief format flag constant
 */
enum fna_format {
	FNA_UNKNOWN		= 0,
	FNA_FASTA		= 1,
	FNA_FASTQ		= 2,
	FNA_FAST5		= 3,
	FNA_GFA			= 4
};

/**
 * @enum fna_options
 */
enum fna_options {
	FNA_SKIP_QUAL 	= 1
};

/**
 * @enum fna_seq_type
 * @brief distinguish struct fna_seq_s with struct fna_link_s
 */
enum fna_seq_type {
	FNA_SEGMENT 	= 1,
	FNA_LINK		= 2
};

/**
 * @enum fna_status
 */
enum fna_status {
	FNA_SUCCESS					= 0,
	FNA_ERROR_FILE_OPEN			= 1,
	FNA_ERROR_UNKNOWN_FORMAT	= 2,
	FNA_ERROR_BROKEN_FORMAT		= 3,
	FNA_ERROR_OUT_OF_MEM		= 4,
	FNA_ERROR_UNSUPPORTED_VERSION = 5,
	FNA_EOF 					= -1
};

/**
 * @struct fna_params_s
 * @brief options
 */
struct fna_params_s {
	uint8_t file_format;		/** see enum fna_format */
	uint8_t seq_encode;			/** see enum fna_flag_encode */
	uint16_t options;			/** see enum fna_options */
	uint16_t head_margin;		/** margin at the head of fna_seq_t */
	uint16_t tail_margin;		/** margin at the tail of fna_seq_t	*/
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */
	uint16_t reserved[2];
	void *lmm;					/** lmm memory manager */
};
typedef struct fna_params_s fna_params_t;

#define FNA_PARAMS(...)			( &((struct fna_params_s const) { __VA_ARGS__ }) )

/**
 * @struct fna_s
 *
 * @brief a struct for fna context
 */
struct fna_s {
	void *reserved1;
	char *path;
	uint8_t file_format;		/** see enum fna_format */
	uint8_t seq_encode;			/** see enum fna_flag_encode */
	uint16_t options;
	int32_t status;				/** see enum fna_status */
	uint32_t reserved2[7];
};
typedef struct fna_s fna_t;

/**
 * @struct fna_str_s
 */
struct fna_str_s {
	char const *ptr;
	int32_t len;
};

/**
 * @struct fna_sarr_s
 */
struct fna_sarr_s {
	uint8_t const *ptr;
	int64_t len;
};

/**
 * @struct fna_cigar_s
 */
struct fna_cigar_s {
	char const *ptr;
	int64_t len;
};

/**
 * @struct fna_segment_s
 */
struct fna_segment_s {
	struct fna_str_s name;
	struct fna_str_s comment;
	struct fna_sarr_s seq;
	struct fna_sarr_s qual;
};

/**
 * @struct fna_link_s
 */
struct fna_link_s {
	struct fna_str_s src;
	struct fna_str_s dst;
	int32_t src_ori;			/** 0: forward, 1: reverse */
	int32_t dst_ori;			/** 0: forward, 1: reverse */
	struct fna_cigar_s cigar;
	int32_t _pad[2];
};

/**
 * @struct fna_seq_s
 *
 * @brief a struct to contain parsed sequence.
 */
struct fna_seq_s {
	void *reserved1;
	uint8_t type;
	uint8_t seq_encode;			/** one of fna_flag_encode */
	uint16_t options;
	union fna_seq_body_u {
		struct fna_segment_s segment;
		struct fna_link_s link;
	} s;
	uint16_t reserved3[4];
};
typedef struct fna_seq_s fna_seq_t;

/**
 * @fn fna_init
 *
 * @brief create a sequence reader context
 *
 * @param[in] path : a path to file to open.
 * @param[in] pack : see struct fna_params_s
 *
 * @return a pointer to the context, NULL if an error occurred (may be invalid path or invalid format)
 */
fna_t *fna_init(char const *path, fna_params_t const *params);

/**
 * @fn fna_close
 *
 * @brief clean up sequence reader context
 */
void fna_close(fna_t *fna);

/**
 * @fn fna_set_lmm
 *
 * @brief replace local memory manager context
 */
void *fna_set_lmm(fna_t *fna, void *lmm);

/**
 * @fn fna_read
 *
 * @brief read a sequence
 *
 * @param[in] fna : a pointer to the context
 *
 * @return a pointer to a sequence object, NULL if the file pointer reached the end.
 */
fna_seq_t *fna_read(fna_t *fna);

/**
 * @fn fna_append
 *
 * @brief concatenate src sesquence after dst sequence
 */
void fna_append(fna_seq_t *dst, fna_seq_t const *src);

/**
 * @fn fna_duplicate
 *
 * @brief duplicate sequence
 */
fna_seq_t *fna_duplicate(fna_seq_t const *seq);

/**
 * @fn fna_append_revcomp
 *
 * @brief append reverse complemented sequence after the given sequence
 */
void fna_append_revcomp(fna_seq_t *seq, fna_seq_t const *src);

/**
 * @fn fna_revcomp
 *
 * @brief make reverse complemented sequence
 */
fna_seq_t *fna_revcomp(fna_seq_t const *seq);

/**
 * @fn fna_seq_free
 *
 * @brief clean up sequence object
 */
void fna_seq_free(fna_seq_t *seq);

#endif /* #ifndef _FNA_H_INCLUDED */
/**
 * end of fna.h
 */
