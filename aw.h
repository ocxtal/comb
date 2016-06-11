
/**
 * @file aw.h
 *
 * @brief alignment writer
 */
#ifndef _AW_H_INCLUDED
#define _AW_H_INCLUDED

#include <stdint.h>


/**
 * @enum aw_file_format
 */
enum aw_file_format {
	AW_SAM = 1,
	AW_BAM = 2,
	AW_MAF = 3,
	AW_GPA = 4			/* graphical pairwise alignment format */
};

/**
 * @enum aw_sam_flags
 */
enum aw_sam_flags {
	/** flags written to the sam file */
	SAM_MULTIPLE_SEGMENTS	= 0x0001,
	SAM_PROPERLY_ALIGNED 	= 0x0002,
	SAM_UNMAPPED 			= 0x0004,
	SAM_NEXT_UNMAPPED 		= 0x0008,
	SAM_REVCOMP				= 0x0010,
	SAM_NEXT_REVCOMP		= 0x0020,
	SAM_FIRST_SEGMENT		= 0x0040,
	SAM_LAST_SEGMENT 		= 0x0080,
	SAM_SECONDARY			= 0x0100,
	SAM_SUPPLEMENTARY		= 0x0800
};

/**
 * @struct aw_params_s
 */
struct aw_params_s {
	uint8_t format;
	char clip;
	uint8_t pad[2];
	uint32_t program_id;
	char const *program_name;
	char const *command;
	char const *name_prefix;
};
typedef struct aw_params_s aw_params_t;
#define AW_PARAMS(...)		( &((aw_params_t const){ __VA_ARGS__ }) )

/**
 * @type aw_t
 */
typedef struct aw_s aw_t;

/**
 * @fn aw_init
 *
 * @brief create a alignment writer context
 */
aw_t *aw_init(
	char const *path,
	gref_idx_t const *idx,
	aw_params_t const *params);

/**
 * @fn aw_clean
 */
void aw_clean(aw_t *aw);

/**
 * @fn aw_append_alignment
 */
void aw_append_alignment(
	aw_t *aw,
	gref_idx_t const *ref,
	gref_acv_t const *query,
	struct gaba_result_s const *const *aln,
	int64_t cnt);

#endif /* _SAM_H_INCLUDED */
/**
 * end of sam.h
 */
