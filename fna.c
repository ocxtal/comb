
/**
 * @file fna.c
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
 *     fna_t *fna_init(char const *path, fna_params_t *params);
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
 *   fna_seq_t (alias to struct _fna_seq): sequence container.
 *     name: sequence name container (lmm_kvec_t(char) instance)
 *       name.a: pointer to the sequence name (null-terminated ASCII).
 *     seq: sequence container (lmm_kvec_t(uint8_t) or kpvec_t(uint8_t) instance)
 *       seq.a: pointer to the sequence (null-terminated when fna->seq_encode == FNA_RAW)
 */

#define UNITTEST_UNIQUE_ID			38
#include "unittest.h"

#include <stdint.h>
#include <string.h>
#include "zf.h"
#include "lmm.h"
#include "log.h"
#include "sassert.h"
#include "fna.h"


/* inline directive */
#define _force_inline				inline

/* roundup */
#define _roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )

/* type aliasing for returning values */
typedef lmm_kvec_t(uint8_t) lmm_kvec_uint8_t;

/**
 * @struct fna_read_ret_s
 */
struct fna_read_ret_s {
	int64_t len;
	char c;
};

/**
 * @struct fna_context_s
 *
 * @brief a struct for fna context
 */
struct fna_context_s {
	lmm_t *lmm;					/** lmm memory manager */
	char *path;
	uint8_t file_format;
	uint8_t seq_encode;
	uint16_t options;
	int32_t status;
	zf_t *fp;					/** zf context (file pointer) */
	uint16_t head_margin;		/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */

	/* file format specific parser */
	struct fna_seq_intl_s *(*read)(struct fna_context_s *fna);

	/* output sequence format specific parser */
	struct fna_read_ret_s (*read_seq)(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);
};
_static_assert_offset(struct fna_s, path, struct fna_context_s, path, 0);
_static_assert_offset(struct fna_s, file_format, struct fna_context_s, file_format, 0);
_static_assert_offset(struct fna_s, seq_encode, struct fna_context_s, seq_encode, 0);
_static_assert_offset(struct fna_s, status, struct fna_context_s, status, 0);

/**
 * @struct fna_seq_intl_s
 *
 * @brief a struct which contains individual sequence.
 */
struct fna_seq_intl_s {
	lmm_t *lmm;					/** lmm memory manager */
	uint8_t type;				/** type, seq or link */
	uint8_t seq_encode;			/** one of _fna_flag_encode */
	uint16_t options;
	union fna_seq_body_intl_u {
		struct {
			struct fna_str_s name;
			struct fna_str_s comment;
			struct fna_sarr_s seq;
			struct fna_sarr_s qual;
		} segment;
		struct {
			struct fna_str_s src;
			int32_t src_ori;			/** 0: forward, 1: reverse */
			struct fna_str_s dst;
			int32_t dst_ori;			/** 0: forward, 1: reverse */
			struct fna_cigar_s cigar;
		} link;

		#if 0
		struct {
			char *name;					/** sequence name */
			char *comment;				/** space separated tail of the name */
			int32_t name_len;
			int32_t comment_len;
			uint8_t *seq;				/** sequence */
			int64_t seq_len;			/** sequence length */
			uint8_t *qual;
			int64_t qual_len;
		} segment;
		struct {
			char *from;
			void *pad;
			int32_t from_len;
			int32_t from_ori;
			char *to;
			int32_t to_len;
			int32_t to_ori;
			char *cigar;				/** contains link cigar string */
			int64_t cigar_len;			/** contains link cigar string length (== strlen(seq)) */
		} link;
		#endif
	} s;
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */
};
_static_assert_offset(struct fna_seq_s, type, struct fna_seq_intl_s, type, 0);
_static_assert_offset(struct fna_seq_s, seq_encode, struct fna_seq_intl_s, seq_encode, 0);
_static_assert_offset(struct fna_seq_s, options, struct fna_seq_intl_s, options, 0);

/* segment members */
_static_assert_offset(struct fna_seq_s, s.segment.name.ptr, struct fna_seq_intl_s, s.segment.name.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.segment.name.len, struct fna_seq_intl_s, s.segment.name.len, 0);
_static_assert_offset(struct fna_seq_s, s.segment.comment.ptr, struct fna_seq_intl_s, s.segment.comment.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.segment.comment.len, struct fna_seq_intl_s, s.segment.comment.len, 0);
_static_assert_offset(struct fna_seq_s, s.segment.seq.ptr, struct fna_seq_intl_s, s.segment.seq.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.segment.seq.len, struct fna_seq_intl_s, s.segment.seq.len, 0);
_static_assert_offset(struct fna_seq_s, s.segment.qual.ptr, struct fna_seq_intl_s, s.segment.qual.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.segment.qual.len, struct fna_seq_intl_s, s.segment.qual.len, 0);

/* link members */
_static_assert_offset(struct fna_seq_s, s.link.src.ptr, struct fna_seq_intl_s, s.link.src.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.link.src.len, struct fna_seq_intl_s, s.link.src.len, 0);
_static_assert_offset(struct fna_seq_s, s.link.dst.ptr, struct fna_seq_intl_s, s.link.dst.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.link.dst.len, struct fna_seq_intl_s, s.link.dst.len, 0);
_static_assert_offset(struct fna_seq_s, s.link.cigar.ptr, struct fna_seq_intl_s, s.link.cigar.ptr, 0);
_static_assert_offset(struct fna_seq_s, s.link.cigar.len, struct fna_seq_intl_s, s.link.cigar.len, 0);


/* function delcarations */
static int fna_read_head_fasta(struct fna_context_s *fna);
static int fna_read_head_fastq(struct fna_context_s *fna);
static int fna_read_head_fast5(struct fna_context_s *fna);
static int fna_read_head_gfa(struct fna_context_s *fna);

static struct fna_seq_intl_s *fna_read_fasta(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fastq(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fast5(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_gfa(struct fna_context_s *fna);

static struct fna_read_ret_s fna_read_seq_ascii(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_2bit(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_2bitpacked(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_4bit(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_4bitpacked(struct fna_context_s *fna, lmm_kvec_uint8_t *v, uint8_t const *delim_table);

/**
 * @fn fna_init
 *
 * @brief create a sequence reader context
 *
 * @param[in] path : a path to a file to open.
 *
 * @return a pointer to the context
 */
fna_t *fna_init(
	char const *path,
	fna_params_t const *params)
{
	struct fna_context_s *fna = NULL;

	/* extension determination */
	struct _ext {
		char const *ext;
		int32_t file_format;
	};
	static struct _ext const ext[] = {
		{".fasta", FNA_FASTA},
		{".fas",   FNA_FASTA},
		{".seq",   FNA_FASTA},
		{".fna",   FNA_FASTA},
		{".ffn",   FNA_FASTA},
		{".fa",    FNA_FASTA},
		{".fastq", FNA_FASTQ},
		{".fq",    FNA_FASTQ},
		{".fast5", FNA_FAST5},
		{".f5",    FNA_FAST5},
		{".gfa",   FNA_GFA},
		{NULL,     0}
	};
	struct _ext const *ep;

	/* read functions */
	int (*read_head[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_head_fasta,
		[FNA_FASTQ] = fna_read_head_fastq,
		[FNA_FAST5] = fna_read_head_fast5,
		[FNA_GFA]	= fna_read_head_gfa
	};
	struct fna_seq_intl_s *(*read[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_fasta,
		[FNA_FASTQ] = fna_read_fastq,
		[FNA_FAST5] = fna_read_fast5,
		[FNA_GFA]	= fna_read_gfa
	};

	/* pack functions */
	struct fna_read_ret_s (*read_seq[])(
		struct fna_context_s *fna,
		lmm_kvec_uint8_t *v,
		uint8_t const *delim_table) = {
		[FNA_ASCII] = fna_read_seq_ascii,
		[FNA_2BIT] = fna_read_seq_2bit,
		[FNA_2BITPACKED] = fna_read_seq_2bitpacked,
		[FNA_4BIT] = fna_read_seq_4bit,
		[FNA_4BITPACKED] = fna_read_seq_4bitpacked,
	};

	/* default params */
	struct fna_params_s default_params = {
		.lmm = NULL,
		.seq_encode = FNA_ASCII,
		.file_format = FNA_UNKNOWN,
		.options = 0,
		.head_margin = 0,
		.tail_margin = 0,
		.seq_head_margin = 0,
		.seq_tail_margin = 0
	};

	if(path == NULL) { return NULL; }
	if(params == NULL) { params = &default_params; }

	/* global context is malloc'd with global malloc */
	if((fna = (struct fna_context_s *)malloc(sizeof(struct fna_context_s))) == NULL) {
		goto _fna_init_error_handler;
	}
	fna->lmm = params->lmm;
	fna->path = NULL;
	fna->fp = NULL;

	/* copy params */
	fna->seq_encode = params->seq_encode;	/** encode sequence to 2-bit if encode == FNA_2BITPACKED */
	fna->file_format = params->file_format;	/** format (see enum FNA_FORMAT) */
	fna->options = params->options;
	fna->head_margin = _roundup(params->head_margin, 16);
	fna->tail_margin = _roundup(params->tail_margin, 16);
	fna->seq_head_margin = _roundup(params->seq_head_margin, 16);
	fna->seq_tail_margin = _roundup(params->seq_tail_margin, 16);

	/* restore defaults */
	if(fna->seq_encode == 0) { fna->seq_encode = FNA_ASCII; }

	/* open file */
	fna->fp = zfopen(path, "r");
	if(fna->fp == NULL) { goto _fna_init_error_handler; }

	/**
	 * if fna->file_format is not specified...
	 * 1. determine file format from the path extension
	 */
	if(fna->file_format == 0) {
		int64_t path_len = strlen(fna->fp->path);
		char const *path_tail = fna->fp->path + path_len;
		for(ep = ext; ep->ext != NULL; ep++) {
			/* skip if path string is shorter than extension string */
			if(path_len < strlen(ep->ext)) { continue; }

			/* compare ext */
			if(strncmp(path_tail - strlen(ep->ext), ep->ext, strlen(ep->ext)) == 0) {
				fna->file_format = ep->file_format; break;
			}
		}
	}

	/**
	 * 2. determine file format from the content of the file
	 */
	if(fna->file_format == 0) {
		/* peek the head of the file */
		char buf[33] = { 0 };
		uint64_t len = zfpeek(fna->fp, buf, 32);
		for(int64_t i = 0; i < len; i++) {
			switch(buf[i]) {
				case '>': fna->file_format = FNA_FASTA; break;
				case '@': fna->file_format = FNA_FASTQ; break;
				case 'H':
				if(buf[i + 1] == '\t') {
					fna->file_format = FNA_GFA; break;
				}
			}
		}
	}
	if(fna->file_format == 0) {
		fna->status = FNA_ERROR_UNKNOWN_FORMAT;
		goto _fna_init_error_handler;
	}
	#ifndef HAVE_HDF5
		if(fna->file_format == FNA_FAST5) {
			// log_error("Fast5 file format is not supported in this build.\n");
			fna->status = FNA_ERROR_UNKNOWN_FORMAT;
			goto _fna_init_error_handler;
		}
	#endif
	fna->read = read[fna->file_format];
	fna->read_seq = read_seq[fna->seq_encode];
	fna->path = strdup(path);

	/* parse header */
	if(read_head[fna->file_format](fna) != FNA_SUCCESS) {
		/* something is wrong */
		goto _fna_init_error_handler;
	}
	return((struct fna_s *)fna);

_fna_init_error_handler:
	if(fna != NULL) {
		zfclose(fna->fp); fna->fp = NULL;
		free(fna->path); fna->path = NULL;
		free(fna);
	}
	return(NULL);
}

/**
 * @fn fna_close
 *
 * @brief clean up sequence reader context
 */
void fna_close(fna_t *ctx)
{
	struct fna_context_s *fna = (struct fna_context_s *)ctx;

	if(fna != NULL) {
		zfclose(fna->fp); fna->fp = NULL;
		free(fna->path); fna->path = NULL;
		free(fna); fna = NULL;
	}
	return;
}

/**
 * @fn fna_set_lmm
 */
void *fna_set_lmm(
	fna_t *ctx,
	void *new)
{
	struct fna_context_s *fna = (struct fna_context_s *)ctx;

	lmm_t *old = fna->lmm;
	fna->lmm = (lmm_t *)new;
	return((void *)old);
}

/**
 * miscellaneous tables and functions
 */

/**
 * @macro _fna_pack_segment
 */
#define _fna_pack_segment(_seg, _name, _name_len, _com, _com_len, _seq, _seq_len, _qual, _qual_len) ({ \
	(_seg)->s.segment.name = (_name); \
	(_seg)->s.segment.name_len = (_name_len); \
	(_seg)->s.segment.comment = (_com); \
	(_seg)->s.segment.comment_len = (_com_len); \
	(_seg)->s.segment.seq = (_seq); \
	(_seg)->s.segment.seq_len = (_seq_len); \
	(_seg)->s.segment.qual = (_qual); \
	(_seg)->s.segment.qual_len = (_qual_len); \
	(_seg); \
})

/**
 * @macro _fna_pack_link
 */
#define _fna_pack_link(_link, _from, _from_len, _from_ori, _to, _to_len, _to_ori, _cigar, _cigar_len) ({ \
	(_link)->s.link.from = (_from); \
	(_link)->s.link.from_len = (_from_len); \
	(_link)->s.link.from_ori = (_from_ori); \
	(_link)->s.link.to = (_to); \
	(_link)->s.link.to_len = (_to_len); \
	(_link)->s.link.to_ori = (_to_ori); \
	(_link)->s.link.cigar = (_cigar); \
	(_link)->s.link.cigar_len = (_cigar_len); \
	(_link); \
})

/**
 * @val ascii_table, alpha_table, delim_space
 */
#define DELIM_TERM			( 1 )

/* non-printable, first 32 elements of the ASCII table */
#define _non_printable(n) \
	(n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), \
	(n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n)

static
uint8_t const delim_space[256] = {
	['\0'] = DELIM_TERM,
	[' '] = DELIM_TERM,
	['\t'] = DELIM_TERM,
	['\v'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_line[256] = {
	['\r'] = DELIM_TERM,
	['\n'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_fasta_fastq_name[256] = {
	[' '] = DELIM_TERM,
	['\r'] = DELIM_TERM,
	['\n'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_fasta_seq[256] = {
	_non_printable(2),
	['>'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_fastq_seq[256] = {
	_non_printable(2),
	['@'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_fastq_qual[256] = {
	_non_printable(2),
	['+'] = DELIM_TERM,
	[0xff] = 0xff
};
static
uint8_t const delim_gfa_field[256] = {
	['\t'] = DELIM_TERM,
	['\r'] = DELIM_TERM,
	['\n'] = DELIM_TERM,
	[0xff] = 0xff
};

/**
 * @fn fna_seq_make_margin
 */
static _force_inline
void fna_seq_make_margin(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	int64_t len)
{
	uint64_t base = lmm_kv_size(*v);
	lmm_kv_reserve(fna->lmm, *v, base + len);
	for(int64_t i = 0; i < len; i++) {
		lmm_kv_at(*v, base + i) = 0;
	}
	lmm_kv_size(*v) = base + len;
	return;
}

/**
 * @fn fna_parse_version_string
 * @brief parse version string in ("%d.%d.%d", major, minor, patch) format,
 * return 0x10000 * major + 0x100 * minor + patch
 */
static _force_inline
uint64_t fna_parse_version_string(
	char const *str)
{
	char buf[256];
	uint64_t v[3] = { 0, 0, 0 };

	for(int64_t i = 0; i < 3; i++) {
		int64_t j = 0;
		while(*str != '\0') {
			if(*str == '.') { break; }
			buf[j++] = *str++;
		}
		buf[j] = '\0';
		v[i] = (uint64_t)strtoll(buf, NULL, 10);
		if(*str++ == '\0') { break; }
	}

	return(0x10000 * v[0] + 0x100 * v[1] + v[2]);
}

/**
 * @fn fna_read_ascii
 * @brief read ascii until delim, push '\0' terminator
 */
static _force_inline
struct fna_read_ret_s fna_read_ascii(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int64_t len = 0;
	int c;

	/* strip spaces at the head */
	while(delim_space[(uint8_t)(c = zfgetc(fna->fp))] == 1) {
		debug("%c, %d", c, c);
	}
	if(c == EOF) {
		debug("reached eof");
		lmm_kv_push(fna->lmm, *v, '\0');
		return((struct fna_read_ret_s){
			.len = len,
			.c = (char)EOF
		});
	}

	/* read line until delim */
	debug("%c, %d", c, c);
	lmm_kv_push(fna->lmm, *v, c); len++;
	while(delim_table[(uint8_t)(c = zfgetc(fna->fp))] == 0) {
		debug("%c, %d", c, c);
		lmm_kv_push(fna->lmm, *v, c); len++;
	}
	char cret = (char)c;

	/* strip spaces at the tail */
	while(len-- > 0 && delim_space[(uint8_t)(c = lmm_kv_pop(fna->lmm, *v))] == 1) {
		debug("%c, %d", c, c);
	}
	lmm_kv_push(fna->lmm, *v, c); len++;	/* push back last non-space char */

	lmm_kv_push(fna->lmm, *v, '\0');		/* push null terminator */
	debug("finished, len(%lld)", len);
	return((struct fna_read_ret_s){
		.len = len,
		.c = cret
	});
}

/**
 * @fn fna_read_skip
 */
static _force_inline
struct fna_read_ret_s fna_read_skip(
	struct fna_context_s *fna,
	uint8_t const *delim_table)
{
	int c;
	while((delim_table[(uint8_t)(c = zfgetc(fna->fp))] & DELIM_TERM) == 0) {
		debug("%c, %d, %u", c, c, delim_table[(uint8_t)c]);
	}
	return((struct fna_read_ret_s){
		.len = 0,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_ascii
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_ascii(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & DELIM_TERM) { break; }
		if(type != 0) { continue; }
		debug("%c, %d", c, c);
		lmm_kv_push(fna->lmm, *v, (uint8_t)c); len++;
	}
	lmm_kv_push(fna->lmm, *v, '\0');
	debug("finished, len(%lld)", len);

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_encode_2bit
 * @brief mapping IUPAC amb. to 2bit encoding
 */
static _force_inline
uint8_t fna_encode_2bit(
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
 * @fn fna_read_seq_2bit
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_2bit(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c = 0;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & DELIM_TERM) { break; }
		if(type != 0) { continue; }
		lmm_kv_push(fna->lmm, *v, fna_encode_2bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_2bitpacked
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_2bitpacked(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c = 0;
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 4x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			int _c; \
			uint8_t _type; \
			while((_type = delim_table[(uint8_t)(_c = zfgetc(_fna->fp))]) != 0) { \
				if(_type & DELIM_TERM) { goto _fna_read_seq_2bitpacked_finish; } \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6);
		lmm_kv_push(fna->lmm, *v, arr); len += 4;

		#undef _fetch
	}
_fna_read_seq_2bitpacked_finish:;
	lmm_kv_push(fna->lmm, *v, arr>>rem); len += (8 - rem) / 2;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_encode_4bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t fna_encode_4bit(
	int c)
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
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn fna_read_seq_4bit
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_4bit(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c = 0;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & DELIM_TERM) { break; }
		if(type != 0) { continue; }
		lmm_kv_push(fna->lmm, *v, fna_encode_4bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_4bitpacked
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_4bitpacked(
	struct fna_context_s *fna,
	lmm_kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c = 0;
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 2x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			int _c; \
			uint8_t _type; \
			while((_type = delim_table[(uint8_t)(_c = zfgetc(_fna->fp))]) != 0) { \
				if(_type & DELIM_TERM) { goto _fna_read_seq_4bitpacked_finish; } \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_4bit(c = _fetch(fna))<<4); rem -= 4;
		arr = (arr>>2) | (fna_encode_4bit(c = _fetch(fna))<<4);
		lmm_kv_push(fna->lmm, *v, arr); len += 2;

		#undef _fetch
	}
_fna_read_seq_4bitpacked_finish:;
	lmm_kv_push(fna->lmm, *v, arr>>rem); len += (8 - rem) / 4;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_head_fasta
 */
static
int fna_read_head_fasta(
	struct fna_context_s *fna)
{
	/* eat '>' at the head */
	fna_read_skip(fna, delim_fasta_seq);
	return(zfeof(fna->fp) ? FNA_EOF :  FNA_SUCCESS);
}

/**
 * @fn fna_read_fasta
 *
 * @brief (internal) fasta format parser
 */
static
struct fna_seq_intl_s *fna_read_fasta(
	struct fna_context_s *fna)
{
	lmm_kvec_uint8_t v;
	lmm_kv_init(fna->lmm, v);

	/* make margin at the head of seq */
	fna_seq_make_margin(fna, &v, fna->head_margin);

	lmm_kv_pusha(fna->lmm, struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.lmm = fna->lmm,
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse name */
	struct fna_read_ret_s n;
	int64_t name_len = (n = fna_read_ascii(fna, &v, delim_fasta_fastq_name)).len;

	/* parse comment after name */
	int64_t com_len = (n.c == ' ')
		? fna_read_ascii(fna, &v, delim_line).len
		: ({ lmm_kv_push(fna->lmm, v, '\0'); 0; });

	/* parse seq */
	fna_seq_make_margin(fna, &v, fna->seq_head_margin);
	int64_t seq_len = (fna->read_seq(fna, &v, delim_fasta_seq)).len;

	/* check termination */
	if(name_len == 0 && com_len == 0 && seq_len == 0) {
		lmm_kv_destroy(fna->lmm, v);
		return(NULL);
	}

	/* make margin at the tail */
	fna_seq_make_margin(fna, &v, fna->seq_tail_margin);
	lmm_kv_push(fna->lmm, v, '\0');			/* qual string terminator */
	fna_seq_make_margin(fna, &v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		lmm_kv_ptr(v) + fna->head_margin);

	#define _next(x)		( (x).ptr + (x).len + 1 )
	r->s.segment.name = (struct fna_str_s){
		.ptr = (char const *)(r + 1),
		.len = name_len
	};
	r->s.segment.comment = (struct fna_str_s){
		.ptr = (char const *)_next(r->s.segment.name),
		.len = com_len
	};
	r->s.segment.seq = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.comment) + r->seq_head_margin),
		.len = seq_len
	};
	r->s.segment.qual = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.seq) + r->seq_tail_margin),
		.len = 0
	};
	#undef _next
	return(r);

	#if 0
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + (seq_len + 1) + r->seq_tail_margin);
	int64_t qual_len = 0;

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
	#endif
}

/**
 * @fn fna_read_head_fastq
 */
static
int fna_read_head_fastq(
	struct fna_context_s *fna)
{
	/* eat '@' at the head */
	fna_read_skip(fna, delim_fastq_seq);
	return(zfeof(fna->fp) ? FNA_EOF :  FNA_SUCCESS);
}

/**
 * @fn fna_read_fastq
 *
 * @brief (internal) fastq format parser
 */
static
struct fna_seq_intl_s *fna_read_fastq(
	struct fna_context_s *fna)
{
	lmm_kvec_uint8_t v;
	lmm_kv_init(fna->lmm, v);

	/* make margin at the head of seq */
	fna_seq_make_margin(fna, &v, fna->head_margin);

	lmm_kv_pusha(fna->lmm, struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.lmm = fna->lmm,
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	#if 0
	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, delim_line).len;
	#endif

	/* parse name */
	struct fna_read_ret_s n;
	int64_t name_len = (n = fna_read_ascii(fna, &v, delim_fasta_fastq_name)).len;

	/* parse comment after name */
	int64_t com_len = (n.c == ' ')
		? fna_read_ascii(fna, &v, delim_line).len
		: ({ lmm_kv_push(fna->lmm, v, '\0'); 0; });

	/* parse seq */
	fna_seq_make_margin(fna, &v, fna->seq_head_margin);
	int64_t seq_len = (fna->read_seq(fna, &v, delim_fastq_qual)).len;
	fna_seq_make_margin(fna, &v, fna->seq_tail_margin);

	/* skip name */
	fna_read_skip(fna, delim_line);

	/* parse qual */
	int64_t qual_len = (((fna->options & FNA_SKIP_QUAL) == 0)
		? fna->read_seq(fna, &v, delim_fastq_seq)
		: fna_read_skip(fna, delim_fastq_seq)).len;
	lmm_kv_push(fna->lmm, v, '\0');			/* push null terminator */

	/* check termination */
	if(name_len == 0 && seq_len == 0) {
		lmm_kv_destroy(fna->lmm, v);
		return(NULL);
	}

	/* make margin at the tail */
	fna_seq_make_margin(fna, &v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		lmm_kv_ptr(v) + fna->head_margin);

	#define _next(x)		( (x).ptr + (x).len + 1 )
	r->s.segment.name = (struct fna_str_s){
		.ptr = (char const *)(r + 1),
		.len = name_len
	};
	r->s.segment.comment = (struct fna_str_s){
		.ptr = (char const *)_next(r->s.segment.name),
		.len = com_len
	};
	r->s.segment.seq = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.comment) + r->seq_head_margin),
		.len = seq_len
	};
	r->s.segment.qual = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.seq) + r->seq_tail_margin),
		.len = qual_len
	};
	#undef _next
	return(r);

	#if 0
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + (seq_len + 1) + r->seq_tail_margin);

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
	#endif
}

/**
 * @fn fna_read_head_fast5
 */
static
int fna_read_head_fast5(
	struct fna_context_s *fna)
{
#ifdef HAVE_HDF5
	/** not implemented yet */
#endif /* HAVE_HDF5 */
	return(FNA_EOF);
}

/**
 * @fn fna_read_fast5
 *
 * @brief (internal) fast5 format parser
 */
static
struct fna_seq_intl_s *fna_read_fast5(
	struct fna_context_s *fna)
{
#ifdef HAVE_HDF5
	/** not implemented yet */
#endif /* HAVE_HDF5 */
	return(NULL);
}

/**
 * @fn fna_read_head_gfa
 */
static
int fna_read_head_gfa(
	struct fna_context_s *fna)
{
	lmm_kvec_uint8_t buf;
	lmm_kv_init(fna->lmm, buf);

	debug("parse gfa header");

	/* read until '\n' */
	fna_read_ascii(fna, &buf, delim_line);

	/* check prefix */
	char const *prefix = "H\tVN:Z:";
	if(strncmp((char const *)lmm_kv_ptr(buf), prefix, strlen(prefix)) != 0) {
		debug("broken");
		return(fna->status = FNA_ERROR_BROKEN_FORMAT);
	}

	/* check version */
	uint64_t version = fna_parse_version_string(&((char const *)lmm_kv_ptr(buf))[strlen(prefix)]);

	/* cleanup */
	lmm_kv_destroy(fna->lmm, buf);

	debug("%llx", version);
	return((version >= 0x10000) ? FNA_SUCCESS : FNA_ERROR_UNSUPPORTED_VERSION);
}

/**
 * @fn fna_read_gfa_seq
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_seq(
	struct fna_context_s *fna)
{
	lmm_kvec_uint8_t v;
	lmm_kv_init(fna->lmm, v);

	/* make margin at the head of seq */
	fna_seq_make_margin(fna, &v, fna->head_margin);

	lmm_kv_pusha(fna->lmm, struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.lmm = fna->lmm,
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, delim_gfa_field).len;

	/* comment is always blank */
	lmm_kv_push(fna->lmm, v, '\0');

	/* parse seq */
	fna_seq_make_margin(fna, &v, fna->seq_head_margin);
	struct fna_read_ret_s ret = fna->read_seq(fna, &v, delim_gfa_field);
	int64_t seq_len = ret.len;

	/* check if optional field remains */
	if(ret.c == '\t') {
		/* skip optional fields */
		fna_read_skip(fna, delim_line);
	}

	/* check termination */
	if(name_len == 0 && seq_len == 0) {
		lmm_kv_destroy(fna->lmm, v);
		return(NULL);
	}

	/* make margin at the tail */
	fna_seq_make_margin(fna, &v, fna->seq_tail_margin);
	lmm_kv_push(fna->lmm, v, '\0');
	fna_seq_make_margin(fna, &v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		lmm_kv_ptr(v) + fna->head_margin);

	#define _next(x)		( (x).ptr + (x).len + 1 )
	r->s.segment.name = (struct fna_str_s){
		.ptr = (char const *)(r + 1),
		.len = name_len
	};
	r->s.segment.comment = (struct fna_str_s){
		.ptr = (char const *)_next(r->s.segment.name),
		.len = 0
	};
	r->s.segment.seq = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.comment) + r->seq_head_margin),
		.len = seq_len
	};
	r->s.segment.qual = (struct fna_sarr_s){
		.ptr = (uint8_t const *)(_next(r->s.segment.seq) + r->seq_tail_margin),
		.len = 0
	};
	#undef _next

	return(r);

	#if 0
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + (seq_len + 1) + r->seq_tail_margin);
	int64_t qual_len = 0;

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
	#endif
}

/**
 * @fn fna_read_gfa_link
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_link(
	struct fna_context_s *fna)
{
	lmm_kvec_uint8_t v;
	lmm_kv_init(fna->lmm, v);

	/* make margin at the head of seq */
	fna_seq_make_margin(fna, &v, fna->head_margin);

	lmm_kv_pusha(fna->lmm, struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.lmm = fna->lmm,
		.type = FNA_LINK,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse from field */
	struct fna_read_ret_s ret_src = fna_read_ascii(fna, &v, delim_gfa_field);
	if(ret_src.c != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* direction */
	int64_t src_ori = (zfgetc(fna->fp) == '+') ? 0 : 1;
	if(zfgetc(fna->fp) != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* parse to field */
	struct fna_read_ret_s ret_dst = fna_read_ascii(fna, &v, delim_gfa_field);
	if(ret_dst.c != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* direction */
	int64_t dst_ori = (zfgetc(fna->fp) == '+') ? 0 : 1;
	if(zfgetc(fna->fp) != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* parse cigar field */
	struct fna_read_ret_s ret_cig = fna_read_ascii(fna, &v, delim_gfa_field);

	/* check if optional field remains */
	if(ret_cig.c == '\t') {
		/* skip optional fields */
		fna_read_skip(fna, delim_line);
	}

	/* make margin at the tail */
	fna_seq_make_margin(fna, &v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		lmm_kv_ptr(v) + fna->head_margin);

	#define _next(x)		( (x).ptr + (x).len + 1 )
	r->s.link.src = (struct fna_str_s){
		.ptr = (char const *)(r + 1),
		.len = ret_src.len
	};
	r->s.link.src_ori = src_ori;
	r->s.link.dst = (struct fna_str_s){
		.ptr = (char const *)_next(r->s.link.src),
		.len = ret_dst.len
	};
	r->s.link.dst_ori = dst_ori;
	r->s.link.cigar = (struct fna_cigar_s){
		.ptr = (char const *)_next(r->s.link.dst),
		.len = ret_cig.len
	};
	#undef _next

	/* replace "*" with "" */
	if(r->s.link.cigar.ptr[0] == '*') {
		((char *)r->s.link.cigar.ptr)[0] = '\0';
		r->s.link.cigar.len = 0;
	}
	return(r);

	#if 0
	char *from = (char *)(r + 1);
	char *to = from + (ret_src.len + 1);
	char *cig = to + (ret_to.len + 1);

	/* replace "*" with "" */
	if(cig[0] == '*') {
		cig[0] = '\0';
		ret_cig.len = 0;
	}

	return(_fna_pack_link(r,
		from, ret_src.len, from_ori,
		to, ret_to.len, to_ori,
		cig, ret_cig.len));
	#endif
}

/**
 * @fn fna_read_gfa_cont
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_cont(
	struct fna_context_s *fna)
{
	/* fixme: ignoring containment information */
	return(NULL);
}

/**
 * @fn fna_read_gfa_path
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_path(
	struct fna_context_s *fna)
{
	/* fixme: ignoring path line */
	return(NULL);
}

/**
 * @fn fna_read_gfa
 *
 * @brief gfa parser driver
 */
static
struct fna_seq_intl_s *fna_read_gfa(
	struct fna_context_s *fna)
{
	int c;
	while((c = zfgetc(fna->fp)) != EOF) {

		/* eat tab after type character */
		if(zfgetc(fna->fp) != '\t') {
			fna->status = FNA_ERROR_BROKEN_FORMAT;
			return(NULL);
		}

		/* examine type */
		switch(c) {
			case 'S': return(fna_read_gfa_seq(fna));
			case 'L': return(fna_read_gfa_link(fna));

			case 'C':	/* fall throught to 'P' */
			case 'P':
			fna_read_skip(fna, delim_line);
			break;

			/*
			case 'C': return(fna_read_gfa_cont(fna));
			case 'P': return(fna_read_gfa_path(fna));
			*/
			
			default:
			/* broken broken broken */
			fna->status = FNA_ERROR_BROKEN_FORMAT;
			return(NULL);
		}

	}
	fna->status = FNA_EOF;
	return(NULL);
}

/**
 * @fn fna_read
 *
 * @brief read a sequence
 *
 * @param[in] fna : a pointer to the context
 *
 * @return a pointer to a sequence object
 */
fna_seq_t *fna_read(fna_t *ctx)
{
	struct fna_context_s *fna = (struct fna_context_s *)ctx;
	if(fna == NULL) { return NULL; }
	return((fna_seq_t *)fna->read(fna));
}

/**
 * @fn fna_seq_free
 *
 * @brief clean up sequence object
 */
void fna_seq_free(fna_seq_t *seq)
{
	struct fna_seq_intl_s *s = (struct fna_seq_intl_s *)seq;
	if(s != NULL) {

		#define _next(x)		( (x).ptr + (x).len + 1 )

		/* free if external mem is used */
		if(s->type == FNA_SEGMENT) {

			/* segment */
			char const *name_base = (char const *)(s + 1);
			char const *comment_base = (char const *)_next(s->s.segment.name);
			uint8_t const *seq_base = (uint8_t const *)_next(s->s.segment.comment) + s->seq_head_margin; 
			uint8_t const *qual_base = (uint8_t const *)_next(s->s.segment.seq) + s->seq_tail_margin; 

			/* segment */
			if(s->s.segment.name.ptr != name_base) {
				lmm_free(s->lmm, (void *)s->s.segment.name.ptr);
			}
			if(s->s.segment.comment.ptr != comment_base) {
				lmm_free(s->lmm, (void *)s->s.segment.comment.ptr);
			}
			if(s->s.segment.seq.ptr != seq_base) {
				lmm_free(s->lmm, (void *)s->s.segment.seq.ptr);
			}
			if(s->s.segment.qual.ptr != qual_base) {
				lmm_free(s->lmm, (void *)s->s.segment.qual.ptr);
			}

			s->s.segment.name.ptr = NULL;
			s->s.segment.comment.ptr = NULL;
			s->s.segment.seq.ptr = NULL;
			s->s.segment.qual.ptr = NULL;

		} else if(s->type == FNA_LINK) {

			/* link */
			char const *src_base = (char const *)(s + 1);
			char const *dst_base = (char const *)_next(s->s.link.src);
			char const *cigar_base = (char const *)_next(s->s.link.dst);

			if(s->s.link.src.ptr != src_base) {
				lmm_free(s->lmm, (void *)s->s.link.src.ptr);
			}
			if(s->s.link.dst.ptr != dst_base) {
				lmm_free(s->lmm, (void *)s->s.link.dst.ptr);
			}
			if(s->s.link.cigar.ptr != cigar_base) {
				lmm_free(s->lmm, (void *)s->s.link.cigar.ptr);
			}

			s->s.link.src.ptr = NULL;
			s->s.link.dst.ptr = NULL;
			s->s.link.cigar.ptr = NULL;
		}

		#undef _next

		/* free context */
		lmm_free(s->lmm, (void *)((uint8_t *)s - s->head_margin));
	}
	return;
}

#if 0
/**
 * @fn fna_base_comp
 * @brief (internal) make complemented base
 */
static _force_inline
uint8_t fna_base_comp(uint8_t c)
{
	switch(c) {
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default: return 'N';
	}
}

/**
 * @fn fna_append
 *
 * @brief concatenate src sesquence after dst sequence
 */
void fna_append(
	fna_seq_t *_dst,
	fna_seq_t const *_src)
{
	struct fna_seq_intl_s *dst = (struct fna_seq_intl_s *)_dst;
	struct fna_seq_intl_s *src = (struct fna_seq_intl_s *)_src;

	int64_t i;
	int64_t size = src->len;

	if(src->seq_encode != dst->seq_encode || src->seq_encode != src->seq_encode) {
		return;
	}
	if(dst->seq_encode == FNA_RAW) {
		(void)lmm_kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->seq_encode == FNA_RAW) {
		lmm_kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			lmm_kv_push(dst->seq, lmm_kv_at(src->seq, i));
		}
		lmm_kv_push(dst->seq, '\0');
	} else if(src->seq_encode == FNA_2BIT) {
		lmm_kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			lmm_kv_push(dst->seq, lmm_kv_at(src->seq, i));
		}
	} else if(src->seq_encode == FNA_2BITPACKED) {
		kpv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kpv_push(dst->seq, kpv_at(src->seq, i));
		}
	}
	return;
}

/**
 * @fn fna_duplicate
 *
 * @brief duplicate sequence
 */
fna_seq_t *fna_duplicate(
	fna_seq_t const *_seq)
{
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)_seq;

	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ seq->head_margin + seq->tail_margin);
	struct fna_seq_intl_s *dup = (struct fna_seq_intl_s *)(ptr + seq->head_margin);
	dup->head_margin = seq->head_margin;		/** set head margin size */
	dup->tail_margin = seq->tail_margin;

	lmm_kv_init(dup->name);
	lmm_kv_init(dup->seq);
	dup->len = 0;
	if(seq->s.segment.seq_encode == FNA_RAW) {
		lmm_kv_push(dup->seq, '\0');
	}

	fna_append((fna_seq_t *)dup, (fna_seq_t *)seq);

	dup->name = seq->s.segment.name;
	dup->name.a = strdup(seq->s.segment.name.a);
	dup->seq_encode = seq->s.segment.seq_encode;
	return((fna_seq_t *)dup);
}

/**
 * @fn fna_append_revcomp
 *
 * @brief append reverse complement of src after dst sequence
 */
void fna_append_revcomp(
	fna_seq_t *_dst,
	fna_seq_t const *_src)
{
	struct fna_seq_intl_s *dst = (struct fna_seq_intl_s *)_dst;
	struct fna_seq_intl_s *src = (struct fna_seq_intl_s *)_src;

	int64_t i;
	int64_t size = src->len;

	if(src->seq_encode != dst->seq_encode || src->seq_encode != src->seq_encode) {
		return;
	}

	if(dst->seq_encode == FNA_RAW) {
		(void)lmm_kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->seq_encode == FNA_RAW) {
		lmm_kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			lmm_kv_push(dst->seq, fna_base_comp(lmm_kv_at(src->seq, size - i - 1)));
		}
		lmm_kv_push(dst->seq, '\0');
	} else if(src->seq_encode == FNA_2BIT) {
		lmm_kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			lmm_kv_push(dst->seq, 0x03 - lmm_kv_at(src->seq, size - i - 1));
		}
	} else if(src->seq_encode == FNA_2BITPACKED) {
		lmm_kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kpv_push(dst->seq, 0x03 - kpv_at(src->seq, size - i - 1));
		}
	}
	return;
}

/**
 * @fn fna_revcomp
 * @brief make reverse complemented sequence
 */
fna_seq_t *fna_revcomp(
	fna_seq_t const *_seq)
{
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)_seq;

	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ seq->head_margin + seq->tail_margin);
	struct fna_seq_intl_s *rev = (struct fna_seq_intl_s *)(ptr + seq->head_margin);
	rev->head_margin = seq->head_margin;		/** set head margin size */
	rev->tail_margin = seq->tail_margin;

	lmm_kv_init(rev->name);
	lmm_kv_init(rev->seq);
	rev->len = 0;
	if(seq->s.segment.seq_encode == FNA_RAW) {
		lmm_kv_push(rev->seq, '\0');
	}

	fna_append_revcomp((fna_seq_t *)rev, (fna_seq_t *)seq);

	rev->name = seq->s.segment.name;
	rev->name.a = strdup(seq->s.segment.name.a);
	rev->seq_encode = seq->s.segment.seq_encode;
	return((fna_seq_t *)rev);
}
#endif

/**
 * unittests
 */
#include <sys/stat.h>

unittest_config(
	.name = "fna",
	.depends_on = {"zf"}
);

/**
 * @fn fdump
 * @brief dump string to file, returns 1 if succeeded
 */
static _force_inline
int fdump(
	char const *filename,
	char const *content)
{
	FILE *fp = fopen(filename, "w");
	int64_t l = fprintf(fp, "%s", content);
	fclose(fp);
	return(l == strlen(content));
}

/**
 * @fn fcmp
 * @brief compare file, returns zero if the content is equivalent to arr
 */
static _force_inline
int fcmp(char const *filename, int64_t size, uint8_t const *arr)
{
	int64_t res;
	struct stat st;
	FILE *fp = NULL;
	uint8_t *buf = NULL;

	if((fp = fopen(filename, "rb")) == NULL) { return(1); }
	fstat(fileno(fp), &st);
	buf = malloc(sizeof(uint8_t) * st.st_size);

	if(fread(buf, sizeof(uint8_t), st.st_size, fp) != st.st_size) {
		return(0);
	}
	res = memcmp(buf, arr, size);
	free(buf);
	fclose(fp);
	return(res == 0);
}

/**
 * @fn unittest_random_base
 */
static _force_inline
char unittest_random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

/**
 * @fn unittest_generate_random_sequence
 */
static _force_inline
char *unittest_generate_random_sequence(
	int64_t len)
{
	char *seq;		/** a pointer to sequence */
	seq = (char *)malloc(sizeof(char) * (len + 1));

	if(seq == NULL) { return NULL; }
	for(int64_t i = 0; i < len; i++) {
		seq[i] = unittest_random_base();
	}
	seq[len] = '\0';
	return seq;
}


/* unittest for parse_version_string */
unittest()
{
	#define assert_parse_version_string(_str, _num) \
		assert(fna_parse_version_string(_str) == (_num), \
			"%u", \
			fna_parse_version_string(_str));

	assert_parse_version_string("0.0.0", 0x000000);
	assert_parse_version_string("0.0.1", 0x000001);
	assert_parse_version_string("0.1.0", 0x000100);
	assert_parse_version_string("1.2.3", 0x010203);
	assert_parse_version_string("100.200.50", 0x64c832);
	assert_parse_version_string("0.0.01", 0x000001);
	assert_parse_version_string("0.0.10", 0x00000a);
	assert_parse_version_string("0.0.15", 0x00000f);

	#undef assert_parse_version_string
}

/**
 * fails when failed to open file
 */
unittest()
{
	char const *fasta_filename = "test.fa";

	remove(fasta_filename);
	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna == NULL);
}

/**
 * basic FASTA parsing
 */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTA format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fasta_filename = "test_fna.fa";
	char const *fasta_content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n"
		">  test2\n\nAAAA\n"
		">\ttest3 \tcomment   \nACGT";
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));

	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTA, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test0") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test1") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 8, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test2") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test3") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "comment") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ACGT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fasta_filename);
}

/**
 * basic FASTQ parsing
 */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTQ format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fastq_filename = "test_fna.fq";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n"
		"@\ttest3\t  comment   \nACGT\n\n+\ttest3\nNNNN";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test0") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "NNNN") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 4, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test1") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 8, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "NNNNNNNN") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 8, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test2") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "NNNN") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 4, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test3") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp(seq->s.segment.comment.ptr, "comment") == 0, "comment(%s)", seq->s.segment.comment.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ACGT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "NNNN") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 4, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fastq_filename);
	return;
}

/* fastq with qual skipping */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTQ format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fastq_filename = "test_fna.fq";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n"
		"@\ttest3\nACGT\n\n+\ttest3\nNNNN";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, FNA_PARAMS(.options = FNA_SKIP_QUAL));
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test0") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 0, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test1") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 8, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 0, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test2") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "AAAA") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 0, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->s.segment.name.ptr, "test3") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ACGT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 4, "len(%lld)", seq->s.segment.seq.len);
	assert(strcmp((char const *)seq->s.segment.qual.ptr, "") == 0, "seq(%s)", (char const *)seq->s.segment.qual.ptr);
	assert(seq->s.segment.qual.len == 0, "len(%lld)", seq->s.segment.qual.len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fastq_filename);
	return;
}

/**
 * gfa parsing
 */
unittest()
{
	char const *gfa_filename = "test_gfa.gfa";
	char const *gfa_content =
		"H	VN:Z:1.0\n"
		"S	11	ACCTT\n"
		"S	12	TCAAGG\n"
		"S	13	CTTGATT\n"
		"L	11	+	12	-	4M\n"
		"L	12	-	13	+	5M\n"
		"L	11	+	13	+	3M\n"
		"P	14	11+,12-,13+	4M,5M\n"
		"S	15	CTTGATT\n";

	assert(fdump(gfa_filename, gfa_content));
	assert(fcmp(gfa_filename, strlen(gfa_content), (uint8_t const *)gfa_content));

	debug("start gfa");
	fna_t *fna = fna_init(gfa_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_GFA, "fna->file_format(%d)", fna->file_format);

	/* segment 11 */
	fna_seq_t *seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->s.segment.name.ptr, "11") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "ACCTT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 5, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* segment 12 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->s.segment.name.ptr, "12") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "TCAAGG") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 6, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* segment 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->s.segment.name.ptr, "13") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "CTTGATT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 7, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* link 11 to 12 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->s.link.src.ptr, "11") == 0, "from(%s)", seq->s.link.src.ptr);
	assert(seq->s.link.src_ori == 0, "src_ori(%d)", seq->s.link.src_ori);
	assert(strcmp(seq->s.link.dst.ptr, "12") == 0, "to(%s)", seq->s.link.dst.ptr);
	assert(seq->s.link.dst_ori == 1, "dst_ori(%d)", seq->s.link.dst_ori);
	assert(strcmp((char const *)seq->s.link.cigar.ptr, "4M") == 0, "cigar(%s)", (char const *)seq->s.link.cigar.ptr);
	fna_seq_free(seq);

	/* link 12 to 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->s.link.src.ptr, "12") == 0, "from(%s)", seq->s.link.src.ptr);
	assert(seq->s.link.src_ori == 1, "src_ori(%d)", seq->s.link.src_ori);
	assert(strcmp(seq->s.link.dst.ptr, "13") == 0, "to(%s)", seq->s.link.dst.ptr);
	assert(seq->s.link.dst_ori == 0, "dst_ori(%d)", seq->s.link.dst_ori);
	assert(strcmp((char const *)seq->s.link.cigar.ptr, "5M") == 0, "cigar(%s)", (char const *)seq->s.link.cigar.ptr);
	fna_seq_free(seq);

	/* link 11 to 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->s.link.src.ptr, "11") == 0, "from(%s)", seq->s.link.src.ptr);
	assert(seq->s.link.src_ori == 0, "src_ori(%d)", seq->s.link.src_ori);
	assert(strcmp(seq->s.link.dst.ptr, "13") == 0, "to(%s)", seq->s.link.dst.ptr);
	assert(seq->s.link.dst_ori == 0, "dst_ori(%d)", seq->s.link.dst_ori);
	assert(strcmp((char const *)seq->s.link.cigar.ptr, "3M") == 0, "cigar(%s)", (char const *)seq->s.link.cigar.ptr);
	fna_seq_free(seq);

	/* skip path line, not implemented yet */

	/* segment 15 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->s.segment.name.ptr, "15") == 0, "name(%s)", seq->s.segment.name.ptr);
	assert(strcmp((char const *)seq->s.segment.seq.ptr, "CTTGATT") == 0, "seq(%s)", (char const *)seq->s.segment.seq.ptr);
	assert(seq->s.segment.seq.len == 7, "len(%lld)", seq->s.segment.seq.len);
	fna_seq_free(seq);

	/* term */
	seq = fna_read(fna);
	assert(seq == NULL, "%p", seq);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(gfa_filename);
	return;
}

/**
 * format detection (fail)
 */
unittest()
{
	char const *fail_filename = "test_fail.txt";
	char const *fail_content = "A quick brown fox jumps over the lazy dog.\n";
	assert(fdump(fail_filename, fail_content));
	assert(fcmp(fail_filename, strlen(fail_content), (uint8_t const *)fail_content));

	/* returns NULL when failed to detect format */
	fna_t *fna = fna_init(fail_filename, NULL);
	assert(fna == NULL, "fna(%p)", fna);

	/** cleanup file */
	remove(fail_filename);
	return;
}

/**
 * format detection from the content of the file (FASTA)
 */
unittest()
{
	char const *fasta_filename = "test_fna.txt";
	char const *fasta_content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n"
		">  test2\n\nAAAA\n";
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));
	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTA, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(fasta_filename);
	return;
}

/**
 * format detection (FASTQ)
 */
unittest()
{
	char const *fastq_filename = "test_fna.txt";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(fastq_filename);
	return;
}

/* format detection */
unittest()
{
	char const *gfa_filename = "test_gfa.txt";
	char const *gfa_content =
		"H	VN:Z:1.0\n"
		"S	11	ACCTT\n"
		"S	12	TCAAGG\n"
		"S	13	CTTGATT\n"
		"L	11	+	12	-	4M\n"
		"L	12	-	13	+	5M\n"
		"L	11	+	13	+	3M\n"
		"P	14	11+,12-,13+	4M,5M\n"
		"S	15	CTTGATT\n";

	assert(fdump(gfa_filename, gfa_content));
	assert(fcmp(gfa_filename, strlen(gfa_content), (uint8_t const *)gfa_content));

	fna_t *fna = fna_init(gfa_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_GFA, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(gfa_filename);
	return;
}

/* large sequence */
unittest()
{
	char const *filename = "test.fa";

	/* dump */
	zf_t *fp = zfopen(filename, "w");
	for(int64_t i = 0; i < 100; i++) {
		char *seq = unittest_generate_random_sequence(100000);
		char *base = seq;

		zfprintf(fp, "> seq%lld\n", i);
		while(*seq != '\0') {
			for(int64_t j = 0; j < 80; j++) {
				zfputc(fp, *seq++);
				if(*seq == '\0') { break; }
			}
			zfputc(fp, '\n');
		}

		free(base);
	}
	zfclose(fp);

	/* read */
	fna_t *fna = fna_init(filename, NULL);
	fna_seq_t *seq = NULL;

	int64_t i = 0;
	while((seq = fna_read(fna)) != NULL) {
		char buf[1024];
		sprintf(buf, "seq%" PRId64 "", i++);

		assert(strcmp(seq->s.segment.name.ptr, buf) == 0, "name(%s, %s)", seq->s.segment.name.ptr, buf);
		assert(seq->s.segment.seq.len == 100000, "len(%lld)", seq->s.segment.seq.len);

		fna_seq_free(seq);
	}

	fna_close(fna);
	remove(filename);
}

#if 0
/**
 * sequence handling
 */
unittest()
{
	char const *fasta_filename = "test_fna_40.fa";
	char const *fasta_content = ">test0\nAACA\n";
	int32_t const margin = 32;
	char const *magic[3] = {
		"The quick brown fox jumps over the lazy dog.",
		"Lorem ipsum dolor sit amet, consectetur adipisicing elit,",
		"ETAOIN SHRDLU CMFWYP VBGKQJ XZ  "
	};
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));

	fna_t *fna = fna_init(fasta_filename,
		FNA_PARAMS(
			.head_margin = margin,
			.tail_margin = margin
		));
	assert(fna != NULL, "fna(%p)", fna);

	fna_seq_t *seq = fna_read(fna);
	/* fill margin with magic */
	memcpy((void *)seq - margin, magic[0], margin);
	memcpy((void *)(seq + 1), magic[0], margin);

	/* duplicate */
	fna_seq_t *dup = fna_duplicate(seq);
	memcpy((void *)dup - margin, magic[1], margin);
	memcpy((void *)(dup + 1), magic[1], margin);
	assert(strcmp(dup->name, "test0") == 0, "name(%s)", dup->name);
	assert(strcmp((char const *)dup->seq, "AACA") == 0, "dup(%s)", (char const *)dup->seq);
	assert(dup->len == 4, "len(%lld)", dup->len);

	/* generate reverse complement */
	fna_seq_t *rev = fna_revcomp(dup);
	memcpy((void *)rev - margin, magic[2], margin);
	memcpy((void *)(rev + 1), magic[2], margin);
	assert(strcmp(rev->name, "test0") == 0, "name(%s)", rev->name);
	assert(strcmp((char const *)rev->seq, "TGTT") == 0, "rev(%s)", (char const *)rev->seq);
	assert(rev->len == 4, "len(%lld)", rev->len);

	/* append */
	fna_append(dup, dup);
	assert(strcmp(dup->name, "test0") == 0, "name(%s)", dup->name);
	assert(strcmp((char const *)dup->seq, "AACAAACA") == 0, "dup(%s)", (char const *)dup->seq);
	assert(dup->len == 8, "len(%lld)", dup->len);

	/* append reverse complement */
	fna_append_revcomp(rev, rev);
	assert(strcmp(rev->name, "test0") == 0, "name(%s)", rev->name);
	assert(strcmp((char const *)rev->seq, "TGTTAACA") == 0, "rev(%s)", (char const *)rev->seq);
	assert(rev->len == 8, "len(%lld)", rev->len);

	/* check margin */
	assert(((struct fna_seq_intl_s *)seq)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);
	assert(((struct fna_seq_intl_s *)dup)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);
	assert(((struct fna_seq_intl_s *)rev)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);

	assert(strncmp((char const *)seq - margin, magic[0], margin) == 0, "");
	assert(strncmp((char const *)dup - margin, magic[1], margin) == 0, "");
	assert(strncmp((char const *)rev - margin, magic[2], margin) == 0, "");

	assert(((struct fna_seq_intl_s *)seq)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);
	assert(((struct fna_seq_intl_s *)dup)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);
	assert(((struct fna_seq_intl_s *)rev)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);

	assert(strncmp((char const *)(seq + 1), magic[0], margin) == 0, "");
	assert(strncmp((char const *)(dup + 1), magic[1], margin) == 0, "");
	assert(strncmp((char const *)(rev + 1), magic[2], margin) == 0, "");

	/* cleanup */
	fna_seq_free(seq);
	fna_seq_free(dup);
	fna_seq_free(rev);

	fna_close(fna);
	remove(fasta_filename);
}
#endif

/**
 * end of fna.c
 */
