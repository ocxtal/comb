
/**
 * @file aw.c
 *
 * @brief alignment writer implementation
 */

#define UNITTEST_UNIQUE_ID			21
#include "unittest.h"

#include <stdint.h>
#include <string.h>
#include "gref.h"
#include "gaba.h"
#include "zf.h"
#include "log.h"
#include "aw.h"


/* inline directive */
#define _force_inline				inline

/* max / min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )

/* constants */
#define SAM_VERSION_STRING			"1.0"
#define SAM_DEFAULT_READGROUP		( 1 )

#define GPA_VERSION_STRING			"0.1"

/**
 * @struct aw_conf_s
 */
struct aw_conf_s {
	char const *ext;
	char const *mode;
	void (*header)(
		aw_t *aw,
		gref_idx_t const *r,
		gref_acv_t const *q);
	void (*body)(
		aw_t *aw,
		gref_idx_t const *r,
		gref_acv_t const *q,
		gaba_alignment_t const *aln);
	void (*footer)(
		aw_t *aw,
		gref_idx_t const *r,
		gref_acv_t const *q);
};

/**
 * @struct aw_s
 */
struct aw_s {
	zf_t *fp;
	struct aw_conf_s conf;		/* function pointers and misc */
	uint8_t format;				/* AW_SAM, AW_GPA, ... */
	char clip;					/* cigar operation to represent clip operation ('S' or 'H') */
	// uint8_t include_qual;		/* include quality string in output if not zero */
	uint8_t pad1[2];
	uint32_t program_id;
	char *program_name;
	char *command;				/* '\t' must be substituted to ' ' */

	/* alignment name in gam format */
	char *aln_name_prefix;
	uint32_t aln_name_len;
	uint32_t pad2;
	int64_t aln_cnt;
};


/* formatter utils */
/**
 * @fn aw_print_str
 */
static _force_inline
void aw_print_str(
	zf_t *fp,
	char const *str,
	uint32_t len)
{
	debug("print_str %p, %s", str, str);

	char const *p = str;
	char const *lim = str + len;
	while(p < lim) {
		zfputc(fp, *p++);
	}
	return;
}

/**
 * @fn aw_print_num
 */
static _force_inline
void aw_print_num(
	zf_t *fp,
	int64_t n)
{
	debug("print_num %lld", n);

	zfprintf(fp, "%lld", n);
	return;
}

/**
 * @union aw_parse_cigar_table_u
 */
union aw_parse_cigar_table_u {
	struct aw_parse_cigar_table_s {
		char str[2];
		uint8_t len;
		uint8_t adv;
	} table;
	uint32_t all;
};

/**
 * @fn parse_get_cigar_elem
 * @brief get cigar element of length i
 */
static _force_inline
union aw_parse_cigar_table_u aw_parse_get_cigar_elem(
	uint8_t i)
{
	static struct aw_parse_cigar_table_s const conv_table[64] = {
		{{0 }, 0, 0}, {{'1'}, 1, 2}, {{'2'}, 1, 2}, {{'3'}, 1, 2},
		{{'4'}, 1, 2}, {{'5'}, 1, 2}, {{'6'}, 1, 2}, {{'7'}, 1, 2},
		{{'8'}, 1, 2}, {{'9'}, 1, 2}, {{'1', '0'}, 2, 3}, {{'1', '1'}, 2, 3},
		{{'1', '2'}, 2, 3}, {{'1', '3'}, 2, 3}, {{'1', '4'}, 2, 3}, {{'1', '5'}, 2, 3},
		{{'1', '6'}, 2, 3}, {{'1', '7'}, 2, 3}, {{'1', '8'}, 2, 3}, {{'1', '9'}, 2, 3},
		{{'2', '0'}, 2, 3}, {{'2', '1'}, 2, 3}, {{'2', '2'}, 2, 3}, {{'2', '3'}, 2, 3},
		{{'2', '4'}, 2, 3}, {{'2', '5'}, 2, 3}, {{'2', '6'}, 2, 3}, {{'2', '7'}, 2, 3},
		{{'2', '8'}, 2, 3}, {{'2', '9'}, 2, 3}, {{'3', '0'}, 2, 3}, {{'3', '1'}, 2, 3},
		{{'3', '2'}, 2, 3}, {{'3', '3'}, 2, 3}, {{'3', '4'}, 2, 3}, {{'3', '5'}, 2, 3},
		{{'3', '6'}, 2, 3}, {{'3', '7'}, 2, 3}, {{'3', '8'}, 2, 3}, {{'3', '9'}, 2, 3},
		{{'4', '0'}, 2, 3}, {{'4', '1'}, 2, 3}, {{'4', '2'}, 2, 3}, {{'4', '3'}, 2, 3},
		{{'4', '4'}, 2, 3}, {{'4', '5'}, 2, 3}, {{'4', '6'}, 2, 3}, {{'4', '7'}, 2, 3},
		{{'4', '8'}, 2, 3}, {{'4', '9'}, 2, 3}, {{'5', '0'}, 2, 3}, {{'5', '1'}, 2, 3},
		{{'5', '2'}, 2, 3}, {{'5', '3'}, 2, 3}, {{'5', '4'}, 2, 3}, {{'5', '5'}, 2, 3},
		{{'5', '6'}, 2, 3}, {{'5', '7'}, 2, 3}, {{'5', '8'}, 2, 3}, {{'5', '9'}, 2, 3},
		{{'6', '0'}, 2, 3}, {{'6', '1'}, 2, 3}, {{'6', '2'}, 2, 3}, {{'6', '3'}, 2, 3}
	};
	return((union aw_parse_cigar_table_u){ .table = conv_table[i] });
}

/**
 * @fn aw_nprintf
 */
static
int aw_cigar_printf(
	void *ctx,
	char const *fmt,
	...)
{
	va_list l;
	va_start(l, fmt);

	zf_t *fp = (zf_t *)ctx;
	uint64_t const fmt_base_len = strlen("%" PRId64 "");

	int64_t len = va_arg(l, int64_t);
	char op = (fmt[fmt_base_len] == '%') ? va_arg(l, int) : fmt[fmt_base_len];

	int adv = 0;
	if(len < 64) {
		union aw_parse_cigar_table_u c = aw_parse_get_cigar_elem(len);
		for(int64_t i = 0; i < c.table.len; i++) {
			zfputc(fp, c.table.str[i]);
		}
		zfputc(fp, op);
		adv = c.table.adv;
	} else {
		adv = zfprintf(fp, "%" PRId64 "%c", len, op);
	}

	va_end(l);
	return(adv);
}


/* sam format writers */

/**
 * @fn sam_write_header
 */
static
void sam_write_header(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q)
{
	/* write header */
	zfprintf(aw->fp, "@HD\tVN:%s\tSO:unsorted\n", SAM_VERSION_STRING);

	/* write reference sequence names */
	int64_t ref_cnt = gref_get_section_count(r);
	for(int64_t i = 0; i < ref_cnt; i++) {
		zfprintf(aw->fp, "@SQ\tSN:%s\tLN:%u\n",
			gref_get_name(r, gref_gid(i, 0)).ptr,
			gref_get_section(r, gref_gid(i, 0))->len);

		debug("i(%lld), gid(%u), name(%s), len(%u)", i,
			gref_get_section(r, gref_gid(i, 0))->gid,
			gref_get_name(r, gref_gid(i, 0)).ptr,
			gref_get_section(r, gref_gid(i, 0))->len);
	}

	/* write readgroup info */
	zfprintf(aw->fp, "@RG\tID:%d\n", SAM_DEFAULT_READGROUP);

	/* program info */
	if(aw->program_name != NULL || aw->command != NULL) {
		zfprintf(aw->fp, "@PG");

		if(aw->program_name != NULL) {
			zfprintf(aw->fp, "\tID:%d\tPN:%s",
				aw->program_id, aw->program_name);
		}

		if(aw->command != NULL) {
			zfprintf(aw->fp, "\tCL:%s", aw->command);
		}
		zfputc(aw->fp, '\n');
	}
	return;
}

/**
 * @fn sam_calc_flags
 */
static _force_inline
int64_t sam_calc_flags(
	gref_idx_t const *r,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr,
	struct gaba_path_section_s const *next)
{	
	return(0);
}

/**
 * @fn sam_print_option_tags
 */
static _force_inline
void sam_print_option_tags(
	aw_t *aw,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr,
	struct gaba_path_s const *path)
{
	/* print alignment score */
	zfprintf(aw->fp, "RG:Z:%d", SAM_DEFAULT_READGROUP);
	return;
}

/**
 * @fn sam_print_seq_forward
 */
static _force_inline
void sam_print_seq_forward(
	aw_t *aw,
	uint8_t const *seq,
	int64_t len)
{
	char const *decode = "NACMGRSVTWYHKDBN";

	uint8_t const *p = seq;
	uint8_t const *lim = seq + len;
	while(p < lim) {
		zfputc(aw->fp, decode[*p++]);
	}
	return;
}

/**
 * @fn sam_print_seq_reverse
 */
static _force_inline
void sam_print_seq_reverse(
	aw_t *aw,
	uint8_t const *seq,
	int64_t len)
{
	char const *decode = "NTGKCYSBAWRDMHVN";

	uint8_t const *p = seq + len;
	uint8_t const *lim = seq;
	while(p > lim) {
		zfputc(aw->fp, decode[*--p]);		
	}
	return;
}

/**
 * @fn sam_print_cigar_forward
 */
static _force_inline
void sam_print_cigar_forward(
	aw_t *aw,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr,
	struct gaba_path_s const *path)
{
	gref_section_t const *bsec = gref_get_section(q, curr->bid);
	int64_t hlen = curr->bpos;
	int64_t tlen = bsec->len - (curr->bpos + curr->blen);

	debug("blen(%u), hlen(%lld), len(%u), tlen(%lld)", curr->blen, hlen, bsec->len, tlen);

	/* print clip at the head */
	if(hlen > 0) {
		zfprintf(aw->fp, "%lld%c", hlen, aw->clip);
	}

	/* print cigar */
	gaba_dp_print_cigar_forward(
		(gaba_dp_fprintf_t)aw_cigar_printf,
		(void *)aw->fp,
		path->array,
		curr->ppos,
		gaba_plen(curr));

	/* print clip at the tail */
	if(tlen > 0) {
		zfprintf(aw->fp, "%lld%c", tlen, aw->clip);
	}
	zfputc(aw->fp, '\t');
	return;
}

/**
 * @fn sam_print_cigar_reverse
 */
static _force_inline
void sam_print_cigar_reverse(
	aw_t *aw,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr,
	struct gaba_path_s const *path)
{
	gref_section_t const *bsec = gref_get_section(q, curr->bid);
	int64_t hlen = curr->bpos;
	int64_t tlen = bsec->len - (curr->bpos + curr->blen);

	debug("blen(%u), hlen(%lld), len(%u), tlen(%lld)", curr->blen, hlen, bsec->len, tlen);

	/* print clip at the head */
	if(tlen > 0) {
		zfprintf(aw->fp, "%lld%c", tlen, aw->clip);
	}

	/* print cigar */
	gaba_dp_print_cigar_reverse(
		(gaba_dp_fprintf_t)aw_cigar_printf,
		(void *)aw->fp,
		path->array,
		curr->ppos,
		gaba_plen(curr));

	/* print clip at the tail */
	if(hlen > 0) {
		zfprintf(aw->fp, "%lld%c", hlen, aw->clip);
	}
	zfputc(aw->fp, '\t');
	return;
}

/**
 * @fn sam_print_seq_qual_forward
 * @brief seqa and seqb are aligned in parallel
 */
static _force_inline
void sam_print_seq_qual_forward(
	aw_t *aw,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr)
{
	/* determine direction and fix seq pointer */
	gref_section_t const *bsec = gref_get_section(q, gref_fw(curr->bid));

	/* include clipped sequence if clipping is 'S' */
	sam_print_seq_forward(aw,
		(aw->clip == 'S') ? bsec->base : &bsec->base[curr->bpos],
		(aw->clip == 'S') ? bsec->len : curr->blen);

	/* print quality string */
	zfprintf(aw->fp, "\t*\t");
	return;
}

/**
 * @fn sam_write_segment_forward
 * @brief seqa and seqb are aligned in parallel
 */
static _force_inline
void sam_write_segment_forward(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	struct gaba_path_s const *path,
	struct gaba_path_section_s const *curr,
	struct gaba_path_section_s const *next)
{
	/* query name */
	aw_print_str(aw->fp,
		gref_get_name(q, curr->bid).ptr,
		gref_get_name(q, curr->bid).len);
	zfputc(aw->fp, '\t');

	/* flags (revcomp indicator) */
	aw_print_num(aw->fp, sam_calc_flags(r, q, curr, next));
	zfputc(aw->fp, '\t');

	/* reference name and pos (name is skipped by default) */
	aw_print_str(aw->fp, 
		gref_get_name(r, curr->aid).ptr,
		gref_get_name(r, curr->aid).len);
	zfputc(aw->fp, '\t');
	aw_print_num(aw->fp, curr->apos + 1);
	zfputc(aw->fp, '\t');

	/* mapping quality */
	aw_print_num(aw->fp, 255);
	zfputc(aw->fp, '\t');

	/* cigar */
	sam_print_cigar_forward(aw, q, curr, path);

	/* ref name and pos of the next section */
	if(next != NULL) {
		aw_print_str(aw->fp,
			gref_get_name(r, next->aid).ptr,
			gref_get_name(r, next->aid).len);
		zfputc(aw->fp, '\t');
		aw_print_num(aw->fp, next->apos);
		zfputc(aw->fp, '\t');
	} else {
		/* tail */
		zfprintf(aw->fp, "*\t0\t");
	}

	/* template length */
	zfprintf(aw->fp, "0\t");

	/* seq and qual */
	sam_print_seq_qual_forward(aw, q, curr);

	/* print option tags */
	sam_print_option_tags(aw, q, curr, path);
	zfputc(aw->fp, '\n');
	return;
}

/**
 * @fn sam_print_seq_qual_reverse
 * @brief seqa and seqb are aligned antiparallel
 */
static _force_inline
void sam_print_seq_qual_reverse(
	aw_t *aw,
	gref_acv_t const *q,
	struct gaba_path_section_s const *curr)
{
	/* determine direction and fix seq pointer */
	gref_section_t const *bsec = gref_get_section(q, gref_rv(curr->bid));
	uint8_t const *lim = gref_get_lim(q);

	if(bsec->base < lim) {
		/* if reverse-complemented sequence is available */
		sam_print_seq_forward(aw,
			(aw->clip == 'S') ? bsec->base : &bsec->base[curr->bpos],
			(aw->clip == 'S') ? bsec->len : curr->blen);
	} else {
		/* if reverse-complemented sequence is not available */
		gref_section_t const *r = gref_get_section(q, gref_fw(curr->bid));
		sam_print_seq_reverse(aw,
			(aw->clip == 'S') ? r->base : &r->base[r->len - (curr->bpos + curr->blen)],
			(aw->clip == 'S') ? r->len : curr->blen);
	}

	/* print quality string */
	zfprintf(aw->fp, "\t*\t");
	return;
}

/**
 * @fn sam_write_segment_reverse
 * @brief seqa and seqb are aligned antiparallel
 */
static _force_inline
void sam_write_segment_reverse(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	struct gaba_path_s const *path,
	struct gaba_path_section_s const *curr,
	struct gaba_path_section_s const *next)
{
	/* query name */
	aw_print_str(aw->fp,
		gref_get_name(q, curr->bid).ptr,
		gref_get_name(q, curr->bid).len);
	zfputc(aw->fp, '\t');

	/* flags (revcomp indicator) */
	aw_print_num(aw->fp, 0x10 | sam_calc_flags(r, q, curr, next));
	zfputc(aw->fp, '\t');

	/* reference name and pos (name is skipped by default) */
	aw_print_str(aw->fp,
		gref_get_name(r, curr->aid).ptr,
		gref_get_name(r, curr->aid).len);
	zfputc(aw->fp, '\t');
	aw_print_num(aw->fp,
		gref_get_section(r, curr->aid)->len - curr->apos - curr->alen + 1);
	zfputc(aw->fp, '\t');

	/* mapping quality */
	aw_print_num(aw->fp, 255);
	zfputc(aw->fp, '\t');

	/* cigar */
	sam_print_cigar_reverse(aw, q, curr, path);

	/* ref name and pos of the next section */
	if(next != NULL) {
		aw_print_str(aw->fp,
			gref_get_name(r, next->aid).ptr,
			gref_get_name(r, next->aid).len);
		zfputc(aw->fp, '\t');
		aw_print_num(aw->fp, next->apos);
		zfputc(aw->fp, '\t');
	} else {
		/* tail */
		zfprintf(aw->fp, "*\t0\t");
	}

	/* template length */
	zfprintf(aw->fp, "0\t");

	/* seq and qual */
	sam_print_seq_qual_reverse(aw, q, curr);

	/* print option tags */
	sam_print_option_tags(aw, q, curr, path);
	zfputc(aw->fp, '\n');
	return;
}

/**
 * @fn sam_write_segment
 */
static _force_inline
void sam_write_segment(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	struct gaba_path_s const *path,
	struct gaba_path_section_s const *curr,
	struct gaba_path_section_s const *next)
{
	if(gref_dir(curr->aid) ^ gref_dir(curr->bid)) {
		sam_write_segment_reverse(aw, r, q, path, curr, next);
	} else {
		sam_write_segment_forward(aw, r, q, path, curr, next);
	}
	return;
}

/**
 * @fn sam_write_alignment
 */
static
void sam_write_alignment(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	gaba_alignment_t const *aln)
{
	debug("slen(%u)", aln->slen);
	for(int64_t i = 0; i < aln->slen - 1; i++) {
		debug("i(%lld), path(%p), &sec[i](%p), &sec[i+1](%p)",
			i, aln->path, &aln->sec[i], &aln->sec[i + 1]);
		sam_write_segment(aw, r, q, aln->path, &aln->sec[i], &aln->sec[i + 1]);
	}

	debug("i(%u), path(%p), &sec[i](%p), &sec[i+1](%p)",
		aln->slen - 1, aln->path, &aln->sec[aln->slen - 1], NULL);
	sam_write_segment(aw, r, q, aln->path, &aln->sec[aln->slen - 1], NULL);
	return;
}


/* gpa format writers */

/**
 * @fn gpa_write_header
 */
static
void gpa_write_header(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q)
{
	/**
	 * GPA header should contain paths to input GFA files,
	 * current signature can't take the path information
	 * (since the gref objects are ignorant of the source
	 * sequence files).
	 */
	zfprintf(aw->fp, "H\tVN:Z:%s\n", GPA_VERSION_STRING);
	return;
}

/**
 * @fn gpa_write_segment
 */
static _force_inline
void gpa_write_segment(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	struct gaba_path_s const *path,
	struct gaba_path_section_s const *sec,
	int head,
	int tail)
{
	/* write tag ('A': alignment) */
	zfprintf(aw->fp, "A\t");

	/* alignment name */
	aw_print_str(aw->fp, aw->aln_name_prefix, aw->aln_name_len);
	aw_print_num(aw->fp, aw->aln_cnt);
	zfputc(aw->fp, '\t');

	/* ref name */
	aw_print_str(aw->fp,
		gref_get_name(r, sec->aid).ptr,
		gref_get_name(r, sec->aid).len);
	zfputc(aw->fp, '\t');

	/* ref pos */
	aw_print_num(aw->fp,
		(gref_dir(sec->aid) == GREF_FW)
			? sec->apos
			: gref_get_section(r, sec->aid)->len - sec->apos);
	zfputc(aw->fp, '\t');

	/* ref len */
	aw_print_num(aw->fp, sec->alen);
	zfputc(aw->fp, '\t');

	/* ref direction */
	zfputc(aw->fp, (gref_dir(sec->aid) == GREF_FW) ? '+' : '-');
	zfputc(aw->fp, '\t');

	/* query name */
	aw_print_str(aw->fp,
		gref_get_name(q, sec->bid).ptr,
		gref_get_name(q, sec->bid).len);
	zfputc(aw->fp, '\t');

	/* query pos */
	aw_print_num(aw->fp,
		(gref_dir(sec->bid) == GREF_FW)
			? sec->bpos
			: gref_get_section(r, sec->bid)->len - sec->bpos);
	zfputc(aw->fp, '\t');

	/* query len */
	aw_print_num(aw->fp, sec->blen);
	zfputc(aw->fp, '\t');

	/* query direction */
	zfputc(aw->fp, (gref_dir(sec->bid) == GREF_FW) ? '+' : '-');
	zfputc(aw->fp, '\t');

	/* cigar string */
	gaba_dp_print_cigar_forward(
		// (gaba_dp_fprintf_t)zfprintf,
		(gaba_dp_fprintf_t)aw_cigar_printf,
		(void *)aw->fp,
		path->array,
		sec->ppos,
		gaba_plen(sec));
	zfputc(aw->fp, '\t');

	/* prev */
	if(head == 0){
		aw_print_str(aw->fp, aw->aln_name_prefix, aw->aln_name_len);
		aw_print_num(aw->fp, aw->aln_cnt - 1);
	} else {
		zfputc(aw->fp, '*');
	}
	zfputc(aw->fp, '\t');

	/* next */
	if(tail == 0){
		aw_print_str(aw->fp, aw->aln_name_prefix, aw->aln_name_len);
		aw_print_num(aw->fp, aw->aln_cnt + 1);
	} else {
		zfputc(aw->fp, '*');
	}
	zfputc(aw->fp, '\t');

	/* optional fields */
	/* mapping quality */
	zfprintf(aw->fp, "MQ:i:%d\n", 255);
	return;
}

/**
 * @fn gpa_write_alignment
 */
static
void gpa_write_alignment(
	aw_t *aw,
	gref_idx_t const *r,
	gref_acv_t const *q,
	gaba_alignment_t const *aln)
{
	debug("slen(%u)", aln->slen);

	for(int64_t i = 0; i < aln->slen; i++) {
		debug("i(%lld), path(%p), &sec[i](%p)", i, aln->path, &aln->sec[i]);
		gpa_write_segment(aw, r, q, aln->path, &aln->sec[i], i == 0, i == (aln->slen - 1));
		aw->aln_cnt++;
	}
	return;
}


/**
 * @fn aw_append_alignment
 */
void aw_append_alignment(
	aw_t *aw,
	gref_idx_t const *ref,
	gref_acv_t const *query,
	struct gaba_alignment_s const *const *aln,
	int64_t cnt)
{
	for(int64_t i = 0; i < cnt; i++) {
		debug("append i(%lld), ref(%p), query(%p), aln[i](%p)", i, ref, query, aln[i]);
		aw->conf.body(aw, ref, query, aln[i]);
	}
	return;
}

/**
 * @fn strdup_rm_tab
 */
char *strdup_rm_tab(
	char const *str)
{
	char *copy = strdup(str);

	for(int64_t i = 0; i < strlen(str); i++) {
		if(copy[i] == '\t') {
			copy[i] = ' ';
		}
	}
	return(copy);
}

/**
 * @fn aw_init
 *
 * @brief initialize alignment writer context
 */
aw_t *aw_init(
	char const *path,
	gref_idx_t const *idx,
	aw_params_t const *params)
{
	/* replace params if null */
	struct aw_params_s default_params = { 0 };
	params = (params != NULL) ? params : &default_params;

	/* malloc context */
	struct aw_s *aw = (struct aw_s *)malloc(sizeof(struct aw_s));
	if(aw == NULL) {
		goto _aw_init_error_handler;
	}
	memset(aw, 0, sizeof(struct aw_s));

	struct aw_conf_s conf[] = {
		[0] = {
			.ext = "-",
			.mode = "w",
			.header = gpa_write_header,
			.body = gpa_write_alignment,
			.footer = NULL
		},
		[AW_SAM] = {
			.ext = ".sam",
			.mode = "w",
			.header = sam_write_header,
			.body = sam_write_alignment,
			.footer = NULL
		},
		[AW_GPA] = {
			.ext = ".gpa",
			.mode = "w",
			.header = gpa_write_header,
			.body = gpa_write_alignment,
			.footer = NULL
		}
	};

	/* detect format */
	if(params->format != 0) {
		aw->conf = conf[params->format];
	} else {
		for(int64_t i = 0; i < sizeof(conf) / sizeof(struct aw_conf_s); i++) {
			/* skip if extension string is not provided */
			if(conf[i].ext == NULL) { continue; }

			/* if path string is shorter than extension string */
			if(strlen(path) < strlen(conf[i].ext)) { continue; }

			if(strncmp(path + strlen(path) - strlen(conf[i].ext), conf[i].ext, strlen(conf[i].ext)) == 0) {
				debug("format detected %s", conf[i].ext);

				aw->conf = conf[i];
			}
		}
	}
	if(aw->conf.ext == NULL) {
		goto _aw_init_error_handler;
	}

	/* copy params */
	aw->format = params->format;
	aw->clip = (params->clip == 'S' || params->clip == 'H') ? params->clip : 'S';
	// aw->include_qual = (params->include_qual == 0) ? 0 : 1;
	aw->program_id = params->program_id;
	aw->program_name = (params->program_name != NULL)
		? strdup_rm_tab(params->program_name) : NULL;
	aw->command = (params->command != NULL)
		? strdup_rm_tab(params->command) : NULL;
	aw->aln_name_prefix = (params->name_prefix != NULL)
		? strdup_rm_tab(params->name_prefix) : NULL;
	aw->aln_name_len = (params->name_prefix != NULL)
		? strlen(aw->aln_name_prefix) : 0;

	/* init name id counter */
	aw->aln_cnt = 0;

	/* open file */
	aw->fp = zfopen(path, aw->conf.mode);
	if(aw->fp == NULL) {
		goto _aw_init_error_handler;
	}

	if(aw->conf.header != NULL) {
		aw->conf.header(aw, idx, NULL);
	}
	return((aw_t *)aw);

_aw_init_error_handler:;
	if(aw != NULL) {
		zfclose(aw->fp); aw->fp = NULL;
		free(aw->program_name); aw->program_name = NULL;
		free(aw->command); aw->command = NULL;
		free(aw->aln_name_prefix); aw->aln_name_prefix = NULL;
	}
	free(aw);
	return(NULL);
}

/**
 * @fn aw_clean
 *
 * @brief flush the content of buffer and close file
 */
void aw_clean(
	aw_t *aw)
{
	if(aw != NULL) {
		if(aw->conf.footer != NULL) {
			aw->conf.footer(aw, NULL, NULL);
		}

		zfclose(aw->fp); aw->fp = NULL;
		free(aw->program_name); aw->program_name = NULL;
		free(aw->command); aw->command = NULL;
		free(aw->aln_name_prefix); aw->aln_name_prefix = NULL;
	}
	free(aw);
	return;
}


/* unittest */
#include <unistd.h>
#include <fcntl.h>

#define _str(x)		x, strlen(x)
#define _seq(x)		(uint8_t const *)(x), strlen(x)

struct aw_unittest_ctx_s {
	gref_idx_t *idx;
	gaba_alignment_t **res;
	int64_t cnt;
};

void *aw_unittest_init(
	void *params)
{
	struct aw_unittest_ctx_s *ctx = (struct aw_unittest_ctx_s *)malloc(
		sizeof(struct aw_unittest_ctx_s));

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

	ctx->idx = idx;


	struct gaba_alignment_s **res = (struct gaba_alignment_s **)malloc(
		3 * sizeof(struct gaba_alignment_s *));
	/* aln 0 */ {
		res[0] = (struct gaba_alignment_s *)malloc(
			  sizeof(struct gaba_alignment_s)
			+ 3 * sizeof(struct gaba_path_section_s)
			+ sizeof(struct gaba_path_s)
			+ 3 * sizeof(uint32_t));

		struct gaba_path_section_s *s = (struct gaba_path_section_s *)(res[0] + 1);
		s[0] = (struct gaba_path_section_s){ 0, 0, 0, 0, 4, 4, 0 /*, 8*/ };
		s[1] = (struct gaba_path_section_s){ 2, 2, 0, 0, 4, 4, 8 /*, 8*/ };
		s[2] = (struct gaba_path_section_s){ 4, 4, 0, 0, 8, 8, 16 /*, 16*/ };

		struct gaba_path_s *p = (struct gaba_path_s *)(s + 3);
		p->len = 32;
		p->array[0] = 0x55555555;
		p->array[1] = 0x01;
		p->array[2] = 0;

		res[0]->sec = s;
		res[0]->path = p;
		res[0]->score = 10;
		res[0]->slen = 3;
		// res[0]->qual = 100;
	}

	/* aln 1 */ {
		res[1] = (struct gaba_alignment_s *)malloc(
			  sizeof(struct gaba_alignment_s)
			+ 3 * sizeof(struct gaba_path_section_s)
			+ sizeof(struct gaba_path_s)
			+ 3 * sizeof(uint32_t));

		struct gaba_path_section_s *s = (struct gaba_path_section_s *)(res[1] + 1);
		s[0] = (struct gaba_path_section_s){ 0, 5, 0, 4, 4, 4, 0 /*, 8*/ };
		s[1] = (struct gaba_path_section_s){ 2, 3, 0, 0, 4, 4, 8 /*, 8*/ };
		s[2] = (struct gaba_path_section_s){ 4, 1, 0, 0, 2, 2, 16 /*, 4*/ };

		struct gaba_path_s *p = (struct gaba_path_s *)(s + 3);
		p->len = 24;
		p->array[0] = 0x01555555;
		p->array[1] = 0;
		p->array[2] = 0;

		res[1]->sec = s;
		res[1]->path = p;
		res[1]->score = 8;
		res[1]->slen = 3;
		// res[1]->qual = 110;
	}

	/* aln 2 */ {
		res[2] = (struct gaba_alignment_s *)malloc(
			  sizeof(struct gaba_alignment_s)
			+ 3 * sizeof(struct gaba_path_section_s)
			+ sizeof(struct gaba_path_s)
			+ 3 * sizeof(uint32_t));

		struct gaba_path_section_s *s = (struct gaba_path_section_s *)(res[2] + 1);
		s[0] = (struct gaba_path_section_s){ 0, 0, 0, 0, 4, 4, 0 /*, 8*/ };
		s[1] = (struct gaba_path_section_s){ 4, 4, 0, 0, 8, 8, 8 /*, 16*/ };

		struct gaba_path_s *p = (struct gaba_path_s *)(s + 3);
		p->len = 24;
		p->array[0] = 0x01555555;
		p->array[1] = 0;
		p->array[2] = 0;

		res[2]->sec = s;
		res[2]->path = p;
		res[2]->score = 6;
		res[2]->slen = 2;
		// res[2]->qual = 90;
	}
	ctx->res = res;
	ctx->cnt = 3;
	return((void *)ctx);
}

void aw_unittest_clean(
	void *_ctx)
{
	struct aw_unittest_ctx_s *ctx = (struct aw_unittest_ctx_s *)_ctx;
	gref_clean(ctx->idx); ctx->idx = NULL;
	free(ctx->res[0]); ctx->res[0] = NULL;
	free(ctx->res[1]); ctx->res[1] = NULL;
	free(ctx->res[2]); ctx->res[2] = NULL;
	free(ctx->res); ctx->res = NULL;
	free(ctx);
	return;
}

#define omajinai() \
	struct aw_unittest_ctx_s *c = (struct aw_unittest_ctx_s *)gctx;

unittest_config(
	.name = "aw",
	.depends_on = { "gaba", "gref", "zf" },
	.init = aw_unittest_init,
	.clean = aw_unittest_clean
);

unittest()
{
	assert(SAM_DEFAULT_READGROUP == 1);
}

/* redirect to stdout */
unittest()
{
	omajinai();

	/* redirect stdout to /dev/null */
	fflush(stdout);
	int b = dup(1), n = open("/dev/null", O_WRONLY);
	dup2(n, 1); close(n);

	aw_t *aw = aw_init("-", c->idx, NULL);
	assert(aw != NULL, "%p", aw);
	aw_clean(aw);

	/* cleanup */
	fflush(stdout);
	dup2(b, 1); close(b);
}

/* sam format writer */
unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, NULL);

	assert(aw != NULL, "%p", aw);

	aw_clean(aw);
	remove(path);
}

unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, NULL);
	aw_clean(aw);

	char const *sam =
		"@HD\tVN:1.0\tSO:unsorted\n"
		"@SQ\tSN:sec0\tLN:4\n"
		"@SQ\tSN:sec1\tLN:4\n"
		"@SQ\tSN:sec2\tLN:8\n"
		"@RG\tID:1\n";
	char *rbuf = (char *)malloc(strlen(sam) + 1);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, strlen(sam) + 1);

	assert(size == strlen(sam), "size(%lld, %lld)", size, strlen(sam));
	assert(memcmp(rbuf, sam, MIN2(size, strlen(sam))) == 0, "%s%s", dump(rbuf, size), dump(sam, strlen(sam)));

	zfclose(fp);
	free(rbuf);
	remove(path);
}

unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, AW_PARAMS(
		.program_name = "hoge"));
	aw_clean(aw);

	char const *sam =
		"@HD\tVN:1.0\tSO:unsorted\n"
		"@SQ\tSN:sec0\tLN:4\n"
		"@SQ\tSN:sec1\tLN:4\n"
		"@SQ\tSN:sec2\tLN:8\n"
		"@RG\tID:1\n"
		"@PG\tID:0\tPN:hoge\n";
	char *rbuf = (char *)malloc(strlen(sam) + 1);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, strlen(sam) + 1);

	assert(size == strlen(sam), "size(%lld, %lld)", size, strlen(sam));
	assert(memcmp(rbuf, sam, MIN2(size, strlen(sam))) == 0, "%s%s", dump(rbuf, size), dump(sam, strlen(sam)));

	zfclose(fp);
	free(rbuf);
	remove(path);
}

unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, AW_PARAMS(
		.command = "--hoge=aaa --fuga=bbb\t--piyo=ccc"));
	aw_clean(aw);

	char const *sam =
		"@HD\tVN:1.0\tSO:unsorted\n"
		"@SQ\tSN:sec0\tLN:4\n"
		"@SQ\tSN:sec1\tLN:4\n"
		"@SQ\tSN:sec2\tLN:8\n"
		"@RG\tID:1\n"
		"@PG\tCL:--hoge=aaa --fuga=bbb --piyo=ccc\n";
	char *rbuf = (char *)malloc(strlen(sam) + 1);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, strlen(sam) + 1);

	assert(size == strlen(sam), "size(%lld, %lld)", size, strlen(sam));
	assert(memcmp(rbuf, sam, MIN2(size, strlen(sam))) == 0, "%s%s", dump(rbuf, size), dump(sam, strlen(sam)));

	zfclose(fp);
	free(rbuf);
	remove(path);
}

/* append alignment */
unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, NULL);
	aw_append_alignment(aw, c->idx, c->idx, (gaba_alignment_t const *const *)c->res, c->cnt);
	aw_clean(aw);

	char const *sam =
		"@HD\tVN:1.0\tSO:unsorted\n"
		"@SQ\tSN:sec0\tLN:4\n"
		"@SQ\tSN:sec1\tLN:4\n"
		"@SQ\tSN:sec2\tLN:8\n"
		"@RG\tID:1\n"
		"sec0\t0\tsec0\t1\t255\t4M\tsec1\t0\t0\tGGRA\t*\tRG:Z:1\n"
		"sec1\t0\tsec1\t1\t255\t4M\tsec2\t0\t0\tMGGG\t*\tRG:Z:1\n"
		"sec2\t0\tsec2\t1\t255\t8M\t*\t0\t0\tACVVGTGT\t*\tRG:Z:1\n"
		"sec2\t16\tsec0\t1\t255\t4M4S\tsec1\t0\t0\tACACBBGT\t*\tRG:Z:1\n"
		"sec1\t16\tsec1\t1\t255\t4M\tsec2\t0\t0\tCCCK\t*\tRG:Z:1\n"
		"sec0\t16\tsec2\t7\t255\t2S2M\t*\t0\t0\tTYCC\t*\tRG:Z:1\n"
		"sec0\t0\tsec0\t1\t255\t4M\tsec2\t0\t0\tGGRA\t*\tRG:Z:1\n"
		"sec2\t0\tsec2\t1\t255\t8M\t*\t0\t0\tACVVGTGT\t*\tRG:Z:1\n";
	char *rbuf = (char *)malloc(1024);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, 1024);

	assert(size == strlen(sam), "size(%lld, %lld)", size, strlen(sam));
	assert(memcmp(rbuf, sam, MIN2(size, strlen(sam))) == 0, "%s%s, %s, %s", dump(rbuf, size), dump(sam, strlen(sam)), rbuf, sam);

	zfclose(fp);
	free(rbuf);
	remove(path);
}

/* append alignment (hard clip) */
unittest()
{
	omajinai();

	char const *path = "./test.sam";
	aw_t *aw = aw_init(path, c->idx, AW_PARAMS(.clip = 'H'));
	aw_append_alignment(aw, c->idx, c->idx, (gaba_alignment_t const *const *)c->res, c->cnt);
	aw_clean(aw);

	char const *sam =
		"@HD\tVN:1.0\tSO:unsorted\n"
		"@SQ\tSN:sec0\tLN:4\n"
		"@SQ\tSN:sec1\tLN:4\n"
		"@SQ\tSN:sec2\tLN:8\n"
		"@RG\tID:1\n"
		"sec0\t0\tsec0\t1\t255\t4M\tsec1\t0\t0\tGGRA\t*\tRG:Z:1\n"
		"sec1\t0\tsec1\t1\t255\t4M\tsec2\t0\t0\tMGGG\t*\tRG:Z:1\n"
		"sec2\t0\tsec2\t1\t255\t8M\t*\t0\t0\tACVVGTGT\t*\tRG:Z:1\n"
		"sec2\t16\tsec0\t1\t255\t4M4H\tsec1\t0\t0\tBBGT\t*\tRG:Z:1\n"
		"sec1\t16\tsec1\t1\t255\t4M\tsec2\t0\t0\tCCCK\t*\tRG:Z:1\n"
		"sec0\t16\tsec2\t7\t255\t2H2M\t*\t0\t0\tTY\t*\tRG:Z:1\n"
		"sec0\t0\tsec0\t1\t255\t4M\tsec2\t0\t0\tGGRA\t*\tRG:Z:1\n"
		"sec2\t0\tsec2\t1\t255\t8M\t*\t0\t0\tACVVGTGT\t*\tRG:Z:1\n";
	char *rbuf = (char *)malloc(1024);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, 1024);

	assert(size == strlen(sam), "size(%lld, %lld)", size, strlen(sam));
	assert(memcmp(rbuf, sam, MIN2(size, strlen(sam))) == 0, "%s%s, %s, %s", dump(rbuf, size), dump(sam, strlen(sam)), rbuf, sam);

	zfclose(fp);
	free(rbuf);
	remove(path);
}


/* gpa format writer */
unittest()
{
	omajinai();

	char const *path = "./test.gpa";
	aw_t *aw = aw_init(path, c->idx, NULL);

	assert(aw != NULL, "%p", aw);

	aw_clean(aw);
	remove(path);
}

/* check header */
unittest()
{
	omajinai();

	char const *path = "./test.gpa";
	aw_t *aw = aw_init(path, c->idx, NULL);
	aw_clean(aw);

	char const *gpa = "H\tVN:Z:0.1\n";
	char *rbuf = (char *)malloc(strlen(gpa) + 1);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, strlen(gpa) + 1);

	assert(size == strlen(gpa), "size(%lld, %lld)", size, strlen(gpa));
	assert(memcmp(rbuf, gpa, MIN2(size, strlen(gpa))) == 0, "%s%s", dump(rbuf, size), dump(gpa, strlen(gpa)));

	zfclose(fp);
	free(rbuf);
	remove(path);
}

/* append alignment */
unittest()
{
	omajinai();

	char const *path = "./test.gpa";
	aw_t *aw = aw_init(path, c->idx, NULL);
	aw_append_alignment(aw, c->idx, c->idx, (gaba_alignment_t const *const *)c->res, c->cnt);
	aw_clean(aw);

	char const *gpa =
		"H\tVN:Z:0.1\n"
		"A\t0\tsec0\t0\t4\t+\tsec0\t0\t4\t+\t4M\t*\t1\tMQ:i:255\n"
		"A\t1\tsec1\t0\t4\t+\tsec1\t0\t4\t+\t4M\t0\t2\tMQ:i:255\n"
		"A\t2\tsec2\t0\t8\t+\tsec2\t0\t8\t+\t8M\t1\t*\tMQ:i:255\n"
		"A\t3\tsec0\t0\t4\t+\tsec2\t4\t4\t-\t4M\t*\t4\tMQ:i:255\n"
		"A\t4\tsec1\t0\t4\t+\tsec1\t4\t4\t-\t4M\t3\t5\tMQ:i:255\n"
		"A\t5\tsec2\t0\t2\t+\tsec0\t4\t2\t-\t2M\t4\t*\tMQ:i:255\n"
		"A\t6\tsec0\t0\t4\t+\tsec0\t0\t4\t+\t4M\t*\t7\tMQ:i:255\n"
		"A\t7\tsec2\t0\t8\t+\tsec2\t0\t8\t+\t8M\t6\t*\tMQ:i:255\n";
	char *rbuf = (char *)malloc(1024);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, 1024);

	assert(size == strlen(gpa), "size(%lld, %lld)", size, strlen(gpa));
	assert(memcmp(rbuf, gpa, MIN2(size, strlen(gpa))) == 0, "%s%s, %s, %s", dump(rbuf, size), dump(gpa, strlen(gpa)), rbuf, gpa);

	zfclose(fp);
	free(rbuf);
	remove(path);
}

/* gpa with prefix */
unittest()
{
	omajinai();

	char const *path = "./test.gpa";
	aw_t *aw = aw_init(path, c->idx, AW_PARAMS( .name_prefix = "aln" ));
	aw_append_alignment(aw, c->idx, c->idx, (gaba_alignment_t const *const *)c->res, c->cnt);
	aw_clean(aw);

	char const *gpa =
		"H\tVN:Z:0.1\n"
		"A\taln0\tsec0\t0\t4\t+\tsec0\t0\t4\t+\t4M\t*\taln1\tMQ:i:255\n"
		"A\taln1\tsec1\t0\t4\t+\tsec1\t0\t4\t+\t4M\taln0\taln2\tMQ:i:255\n"
		"A\taln2\tsec2\t0\t8\t+\tsec2\t0\t8\t+\t8M\taln1\t*\tMQ:i:255\n"
		"A\taln3\tsec0\t0\t4\t+\tsec2\t4\t4\t-\t4M\t*\taln4\tMQ:i:255\n"
		"A\taln4\tsec1\t0\t4\t+\tsec1\t4\t4\t-\t4M\taln3\taln5\tMQ:i:255\n"
		"A\taln5\tsec2\t0\t2\t+\tsec0\t4\t2\t-\t2M\taln4\t*\tMQ:i:255\n"
		"A\taln6\tsec0\t0\t4\t+\tsec0\t0\t4\t+\t4M\t*\taln7\tMQ:i:255\n"
		"A\taln7\tsec2\t0\t8\t+\tsec2\t0\t8\t+\t8M\taln6\t*\tMQ:i:255\n";
	char *rbuf = (char *)malloc(1024);

	zf_t *fp = zfopen(path, "r");
	int64_t size = zfread(fp, rbuf, 1024);

	assert(size == strlen(gpa), "size(%lld, %lld)", size, strlen(gpa));
	assert(memcmp(rbuf, gpa, MIN2(size, strlen(gpa))) == 0, "%s%s, %s, %s", dump(rbuf, size), dump(gpa, strlen(gpa)), rbuf, gpa);

	zfclose(fp);
	free(rbuf);
	remove(path);
}

/**
 * end of aw.c
 */
