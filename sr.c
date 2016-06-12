
/**
 * @file sr.c
 *
 * @brief sequence reader impl.
 */

#define UNITTEST_UNIQUE_ID			25
#include "unittest.h"

#include <stdint.h>
#include <stdlib.h>
#include "lmm.h"
#include "sr.h"
#include "fna.h"
#include "gref.h"
#include "log.h"
#include "sassert.h"


/* constants */
#define SR_SINGLE_READ_MEM_SIZE		( 4096 )

/* inline directive */
#define _force_inline				inline


/* assertions */
_static_assert((int32_t)FNA_UNKNOWN == (int32_t)SR_UNKNOWN);
_static_assert((int32_t)FNA_FASTA == (int32_t)SR_FASTA);
_static_assert((int32_t)FNA_FASTQ == (int32_t)SR_FASTQ);
_static_assert((int32_t)FNA_FAST5 == (int32_t)SR_FAST5);
_static_assert((int32_t)FNA_GFA == (int32_t)SR_GFA);

_static_assert((int32_t)GREF_FW_ONLY == (int32_t)SR_FW_ONLY);
_static_assert((int32_t)GREF_FW_RV == (int32_t)SR_FW_RV);


/**
 * @struct sr_s
 */
struct sr_s {
	char *path;
	fna_t *fna;
	gref_acv_t *acv;
	gref_idx_t *idx;
	struct sr_gref_s *(*iter_read)(
		sr_t *sr);
	struct sr_params_s params;
};

/**
 * @struct sr_gref_intl_s
 * @brief gref and iter container
 */
struct sr_gref_intl_s {
	lmm_t *lmm;
	char const *path;
	gref_t const *gref;
	gref_iter_t *iter;
	fna_seq_t *seq;
	uint8_t gref_need_free;
	uint8_t seq_need_free;
	uint8_t reserved2[6];
};

/**
 * @fn sr_dump_seq
 */
static _force_inline
void sr_dump_seq(
	sr_t *sr)
{
	/* init pool object */
	gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
		.k = sr->params.k,
		.seq_direction = sr->params.seq_direction,
		.seq_format = GREF_4BIT,
		.seq_head_margin = 32,
		.seq_tail_margin = 32,
		.copy_mode = GREF_COPY,
		.num_threads = sr->params.num_threads,
		.lmm = NULL));

	/* dump sequence */
	fna_seq_t *seq = NULL;
	while((seq = fna_read(sr->fna)) != NULL) {

		if(seq->type == FNA_SEGMENT) {
			/*
			for(int64_t i = 0; i < seq->s.segment.seq_len; i++) {
				fprintf(stderr, "%c", " AC G   T"[seq->s.segment.seq[i]]);
			}
			fprintf(stderr, "\n");
			*/

			gref_append_segment(pool,
				seq->s.segment.name,
				seq->s.segment.name_len,
				seq->s.segment.seq,
				seq->s.segment.seq_len);
		} else if(seq->type == FNA_LINK) {
			/* check cigar starts from '0' (indicating cigar is "0M") */
			if(seq->s.link.cigar[0] != '0') {
				debug("overlapping link is not supported.");
				break;
			} else {
				gref_append_link(pool,
					seq->s.link.from,
					seq->s.link.from_len,
					seq->s.link.from_ori,
					seq->s.link.to,
					seq->s.link.to_len,
					seq->s.link.to_ori);
			}
		} else {
			/* unknown type */
			debug("unknown sequence type appeared.");
		}

		fna_seq_free(seq);
	}
	
	/* freeze */
	sr->acv = gref_freeze_pool(pool);

	/* close fna context */
	fna_close(sr->fna); sr->fna = NULL;
	return;
}

/**
 * @fn sr_get_index
 */
struct sr_gref_s *sr_get_index(
	sr_t *sr)
{
	/* check if index is already built */
	if(sr->idx == NULL) {
		/* archive and build index */
		if(sr->acv == NULL) {
			sr_dump_seq(sr);
		}
		sr->idx = gref_build_index(sr->acv);

		if(sr->idx == NULL) {
			sr->acv = NULL;
			return(NULL);
		}
	}

	struct sr_gref_intl_s *r = (struct sr_gref_intl_s *)malloc(
		sizeof(struct sr_gref_intl_s));
	*r = (struct sr_gref_intl_s){
		.lmm = NULL,
		.path = sr->path,
		.gref = (gref_t const *)sr->idx,
		.iter = NULL,
		.gref_need_free = 0
	};
	return((struct sr_gref_s *)r);
}

/**
 * @fn sr_get_iter_graph
 */
static
struct sr_gref_s *sr_get_iter_graph(
	sr_t *sr)
{
	/* check if archive is already built */
	if(sr->acv == NULL) {
		sr_dump_seq(sr);
		if(sr->acv == NULL) {
			return(NULL);
		}
	} else {
		return(NULL);
	}

	struct sr_gref_intl_s *r = (struct sr_gref_intl_s *)malloc(
		sizeof(struct sr_gref_intl_s));
	*r = (struct sr_gref_intl_s){
		.lmm = NULL,
		.path = sr->path,
		.gref = sr->acv,
		.iter = gref_iter_init(sr->acv)
	};
	return((struct sr_gref_s *)r);
}

/**
 * @fn sr_get_iter_read
 */
static
struct sr_gref_s *sr_get_iter_read(
	sr_t *sr)
{
	if(sr->fna == NULL) {
		return(NULL);
	}

	/* gref objects */
	gref_acv_t *acv = NULL;

	/* init local memory */
	lmm_t *lmm_read = lmm_init(NULL, SR_SINGLE_READ_MEM_SIZE);
	struct sr_gref_intl_s *r = (struct sr_gref_intl_s *)lmm_malloc(
		lmm_read, sizeof(struct sr_gref_intl_s));

	/* read a sequence */
	fna_set_lmm(sr->fna, lmm_read);
	fna_seq_t *seq = NULL;
	while((seq = fna_read(sr->fna)) != NULL) {
		if(seq->type == FNA_SEGMENT) {
			/*
			for(int64_t i = 0; i < seq->s.segment.seq_len; i++) {
				fprintf(stderr, "%c", " AC G   T"[seq->s.segment.seq[i]]);
			}
			fprintf(stderr, "\n");
			*/

			gref_pool_t *pool = gref_init_pool(GREF_PARAMS(
				.k = sr->params.k,
				.seq_direction = sr->params.seq_direction,
				.seq_format = GREF_4BIT,
				.copy_mode = GREF_NOCOPY,
				.num_threads = sr->params.num_threads,
				.lmm = lmm_read));

			gref_append_segment(pool,
				seq->s.segment.name,
				seq->s.segment.name_len,
				seq->s.segment.seq,
				seq->s.segment.seq_len);
			acv = gref_freeze_pool(pool);
			break;
		}
		fna_seq_free(seq);
	}

	if(seq == NULL) {
		lmm_clean(lmm_read); lmm_read = NULL;
		fna_close(sr->fna); sr->fna = NULL;
		return(NULL);
	}

	*r = (struct sr_gref_intl_s){
		.lmm = lmm_read,
		.path = sr->path,
		.gref = acv,
		.iter = gref_iter_init(acv),
		.seq = seq,
		.gref_need_free = 1,
		.seq_need_free = 1
	};
	return((struct sr_gref_s *)r);
}

/**
 * @fn sr_get_iter
 */
struct sr_gref_s *sr_get_iter(
	sr_t *sr)
{
	return(sr->iter_read(sr));
}

/**
 * @fn sr_gref_free
 */
void sr_gref_free(
	struct sr_gref_s *_r)
{
	struct sr_gref_intl_s *r = (struct sr_gref_intl_s *)_r;
	if(r == NULL) { return; }

	gref_iter_clean(r->iter); r->iter = NULL;
	if(r->gref_need_free != 0) { gref_clean((gref_t *)r->gref); } r->gref = NULL;
	if(r->seq_need_free != 0) { fna_seq_free(r->seq); r->seq = NULL; }

	lmm_t *lmm = r->lmm; r->lmm = NULL;
	lmm_free(lmm, r);
	lmm_clean(lmm);
	return;
}

/**
 * @fn sr_init
 */
sr_t *sr_init(
	char const *path,
	sr_params_t const *params)
{
	/* check arguments sanity */
	if(path == NULL) { return(NULL); }

	struct sr_params_s const default_params = { 0 };
	params = (params == NULL) ? &default_params : params;

	/* malloc mem */
	struct sr_s *sr = (struct sr_s *)malloc(sizeof(struct sr_s));
	memset(sr, 0, sizeof(struct sr_s));
	sr->path = strdup(path);

	/* create fna object */
	sr->fna = fna_init(path, FNA_PARAMS(
		.file_format = params->format,
		.seq_encode = FNA_4BIT,
		.seq_head_margin = 32,
		.seq_tail_margin = 32));
	if(sr->fna == NULL) {
		goto _sr_init_error_handler;
	}

	static struct sr_gref_s *(*iter_read_table[])(sr_t *) = {
		[SR_FASTA] = sr_get_iter_read,
		[SR_FASTQ] = sr_get_iter_read,
		[SR_FAST5] = sr_get_iter_read,
		[SR_GFA] = sr_get_iter_graph
	};
	sr->iter_read = iter_read_table[sr->fna->file_format];

	/* copy params */
	sr->params = *params;
	return((sr_t *)sr);

_sr_init_error_handler:;
	free(sr->path); sr->path = NULL;
	fna_close(sr->fna); sr->fna = NULL;
	free(sr);
	return(NULL);
}

/**
 * @fn sr_clean
 */
void sr_clean(
	sr_t *sr)
{
	if(sr == NULL) { return; }

	gref_clean(sr->acv); sr->acv = NULL;
	free(sr->path); sr->path = NULL;
	fna_close(sr->fna); sr->fna = NULL;
	free(sr);
	return;
}


/* unittests */
#include <sys/stat.h>

unittest_config(
	.name = "sr",
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

unittest()
{
	char const *fasta_filename = "test.fa";

	remove(fasta_filename);
	sr_t *sr = sr_init(fasta_filename, NULL);
	assert(sr == NULL);
}

unittest()
{
	char const *fasta_filename = "test.fa";
	char const *fasta_content =
		">test1\n"
		"ACGTACGT\n"
		">test2\n"
		"TTTTGGGG\n"
		">test3\n"
		"AAAAAAAA";

	fdump(fasta_filename, fasta_content);
	sr_t *sr = sr_init(fasta_filename, NULL);

	assert(sr != NULL);

	sr_clean(sr);
	remove(fasta_filename);
}

unittest()
{
	char const *fasta_filename = "test.fa";
	char const *fasta_content =
		">test1\n"
		"ACGTACGT\n"
		">test2\n"
		"TTTTGGGG\n"
		">test3\n"
		"AAAAAAAA";

	fdump(fasta_filename, fasta_content);
	sr_t *sr = sr_init(fasta_filename, SR_PARAMS(
		.k = 3,
		.seq_direction = SR_FW_ONLY));
	assert(sr != NULL);

	struct sr_gref_s *idx = sr_get_index(sr);
	assert(idx != NULL);
	assert(idx->lmm == NULL);
	assert(strcmp(idx->path, fasta_filename) == 0);
	assert(idx->gref != NULL);
	assert(gref_get_total_len(idx->gref) == 24, "%lld", gref_get_total_len(idx->gref));
	assert(idx->iter == NULL);

	sr_gref_free(idx);
	sr_clean(sr);
	remove(fasta_filename);
}

unittest()
{
	char const *fasta_filename = "test.fa";
	char const *fasta_content =
		">test1\n"
		"ACGTACGT\n"
		">test2\n"
		"TTTTGGGG\n"
		">test3\n"
		"AAAAAAAA";

	fdump(fasta_filename, fasta_content);
	sr_t *sr = sr_init(fasta_filename, SR_PARAMS(
		.k = 3,
		.seq_direction = SR_FW_ONLY));
	assert(sr != NULL);

	/* test1 */
	struct sr_gref_s *iter = sr_get_iter(sr);
	assert(iter != NULL);
	assert(iter->lmm != NULL);
	assert(strcmp(iter->path, fasta_filename) == 0);
	assert(iter->gref != NULL);
	assert(iter->iter != NULL);
	sr_gref_free(iter);

	iter = sr_get_iter(sr);
	assert(iter != NULL);
	assert(iter->gref != NULL);
	assert(iter->iter != NULL);
	sr_gref_free(iter);

	iter = sr_get_iter(sr);
	assert(iter != NULL);
	assert(iter->gref != NULL);
	assert(iter->iter != NULL);
	sr_gref_free(iter);

	iter = sr_get_iter(sr);
	assert(iter == NULL);

	sr_clean(sr);
	remove(fasta_filename);
}


/**
 * end of sr.c
 */
