
/**
 * @file psort_radix_intl.c
 *
 * @brief write-combining radix sort implementation. included from psort.c
 *
 * @author Hajime Suzuki
 * @date 2016/3/20
 * @license MIT
 */

/* check if element type and read / write macros are defined */
#if !defined(elem_t) || !defined(rd) || !defined(wr) || !defined(ex)
#  error "elem_t, rd, wr, and ex macros must be defined."
#endif

/* constant */
#define WCR_BUF_ELEM_COUNT		( WCR_BUF_SIZE / sizeof(elem_t) )

/**
 * @fn psort_count_occ
 */
static
void join(psort_count_occ_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* initialize occ array */
	uint64_t *occ = (uint64_t *)ctx->occ->occ;
	memset(occ, 0, sizeof(uint64_t) * WCR_OCC_SIZE);	

	/* extract elem pointer */
	elem_t *src = (elem_t *)ctx->src;
	int64_t from = ctx->from;
	int64_t to = ctx->to;

	/* count occurrences */
	int64_t digit = ctx->digit;
	for(int64_t i = from; i < to; i++) {
		occ[ex(rd(src + i), digit)]++;
		debug("n(%llu), occ(%llu)", ex(rd(src + i), digit), occ[ex(rd(src + i), digit)]);
	}
	return;
}

/**
 * @fn psort_gather_occ
 * @brief gather occurrence table. item is ignored.
 */
static
void join(psort_gather_occ_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	struct psort_occ_s *occ = ctx->occ;
	struct psort_buffer_counter_s *cnt = ctx->cnt;

	int64_t threads = ctx->num_threads;
	uint64_t sum = 0;
	for(int64_t i = 0; i < WCR_OCC_SIZE; i++) {
		for(int64_t j = 0; j < threads; j++) {
			uint64_t curr_occ = occ[j].occ[i];

			/* store base index to occ[j].occ[i] */
			occ[j].occ[i] = sum & ~(WCR_BUF_ELEM_COUNT - 1);

			/* store initial buffer counter to cnt[j].cnt[i] */
			cnt[j].cnt[i] = sum & (WCR_BUF_ELEM_COUNT - 1);

			debug("sum(%llu), occ(%llu), cnt(%u)", sum, occ[j].occ[i], cnt[j].cnt[i]);

			/* update sum */
			sum += curr_occ;
		}
	}
	return;
}

/**
 * @fn psort_scatter
 */
static
void join(psort_scatter_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* extract pointers */
	uint64_t *occ = (uint64_t *)ctx->occ->occ;
	uint8_t *cnt = (uint8_t *)ctx->cnt->cnt;
	elem_t (*buf)[WCR_BUF_ELEM_COUNT] =
		(elem_t (*)[WCR_BUF_ELEM_COUNT])ctx->buf->buf;

	/* extract elem pointer */
	elem_t *src = (elem_t *)ctx->src;
	elem_t *dst = (elem_t *)ctx->dst;
	int64_t from = ctx->from;
	int64_t to = ctx->to;

	debug("from(%lld), to(%lld)", from, to);

	/* scatter */
	int64_t digit = ctx->digit;
	for(int64_t i = from; i < to; i++) {
		/* load an element, store to the buffer */
		elem_t register e = rd(src + i);
		uint64_t register n = ex(e, digit);
		uint8_t register c = cnt[n];

		wr(buf[n] + c, e);

		debug("e(%llu), n(%llu), cnt(%u), p(%p), e(%llu)",
			(uint64_t)e, n, cnt[n], buf[n] + (int64_t)cnt[n],
			(uint64_t)rd(buf[n] + (int64_t)cnt[n]));


		/** check if flush is needed */
		if(++c == WCR_BUF_ELEM_COUNT) {
			debug("bulk copy n(%llu)", n);
			/* bulk copy */
			memcpy_buf(dst + occ[n], buf[n]);

			/* update index */
			occ[n] += WCR_BUF_ELEM_COUNT;

			/* reset counter */
			c = 0;
		}

		/* write back cnt */
		cnt[n] = c;
	}
	return;
}

/**
 * @fn psort_flush
 */
static
void join(psort_flush_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	struct psort_occ_s *occ = ctx->occ;
	struct psort_buffer_counter_s *cnt = ctx->cnt;
	struct psort_buffer_s *buf = ctx->buf;
	elem_t *dst = (elem_t *)ctx->dst;

	/** flush the remaining content */
	int64_t threads = ctx->num_threads;
	int64_t prev_occ = 0;
	int64_t prev_cnt = 0;
	for(int64_t i = 0; i < WCR_OCC_SIZE; i++) {
		for(int64_t j = 0; j < threads; j++) {
			int64_t curr_occ = occ[j].occ[i];
			int64_t curr_cnt = cnt[j].cnt[i];

			debug("prev_occ(%lld), prev_cnt(%lld), curr_occ(%lld), curr_cnt(%lld)",
				prev_occ, prev_cnt, curr_occ, curr_cnt);

			/* reset prev_cnt if the current block differs from the previous one */
			prev_cnt = (prev_occ == curr_occ) ? prev_cnt : 0;
			debug("update prev_cnt(%lld)", prev_cnt);

			/* copy */
			for(int64_t k = prev_cnt; k < curr_cnt; k++) {
				debug("flush move k(%lld), e(%llu), p(%p)",
					k, (uint64_t)rd(((elem_t *)buf[j].buf[i]) + k), ((elem_t *)buf[j].buf[i]) + k);
				wr(dst + curr_occ + k,
					rd(((elem_t *)buf[j].buf[i]) + k));
			}

			/* save occ and cnt */
			prev_occ = curr_occ;
			prev_cnt = curr_cnt;
		}
	}
	return;
}

/**
 * @fn psort_copyback
 */
static
void join(psort_copyback_, SUFFIX)(
	struct psort_thread_context_s *ctx)
{
	/* extract elem pointer */
	elem_t *src = (elem_t *)ctx->src;
	elem_t *dst = (elem_t *)ctx->dst;
	int64_t from = ctx->from;
	int64_t to = ctx->to;

	/* copy back */
	memcpy(dst + from, src + from, sizeof(elem_t) * (to - from));
	return;
}

/**
 * @fn psort_partialsort_parallel
 */
static
void join(psort_partialsort_parallel_, SUFFIX)(
	void *src,
	int64_t len,
	int64_t num_threads,
	int64_t lower_digit,
	int64_t higher_digit)
{
	int64_t nt = (num_threads == 0) ? 1 : num_threads;

	/* malloc buffer */
	void *ptr = aligned_malloc(
		  nt * (
		  	  sizeof(struct psort_thread_context_s *)
		  	+ 3 * sizeof(void *)
			+ sizeof(struct psort_thread_context_s)
			+ sizeof(struct psort_occ_s)
			+ sizeof(struct psort_buffer_counter_s)
			+ sizeof(struct psort_buffer_s))
		+ sizeof(elem_t) * len,
		16);

	/* array of pointers to thread contexts */
	struct psort_thread_context_s **pth = (struct psort_thread_context_s **)ptr;
	
	/* array of pointers to functions */
	void **count_occ = (void **)&pth[nt];
	void **scatter = (void **)&count_occ[nt];
	void **copyback = (void **)&scatter[nt];

	/* thread contexts and working buffers */
	struct psort_thread_context_s *th = (struct psort_thread_context_s *)&copyback[nt];
	struct psort_occ_s *occ = (struct psort_occ_s *)&th[nt];
	struct psort_buffer_counter_s *cnt = (struct psort_buffer_counter_s *)&occ[nt];
	struct psort_buffer_s *buf = (struct psort_buffer_s *)&cnt[nt];
	void *dst = (void *)&buf[nt];

	/* initialize thread contexts */
	for(int64_t i = 0; i < nt; i++) {
		/* pointer to thread context */
		pth[i] = &th[i];

		/* pointer to functions */
		count_occ[i] = (void *)join(psort_count_occ_, SUFFIX);
		scatter[i] = (void *)join(psort_scatter_, SUFFIX);
		copyback[i] = (void *)join(psort_copyback_, SUFFIX);

		/* pointer to working buffers */
		th[i].occ = &occ[i];
		th[i].cnt = &cnt[i];
		th[i].buf = &buf[i];

		/* num_threads */
		th[i].num_threads = nt;

		/* array */
		th[i].from = i * len / nt;
		th[i].to = (i + 1) * len / nt;
	}

	/* initialize ptask object */
	ptask_t *pt = ptask_init(psort_dispatcher, (void **)pth, num_threads, 1024);

	/* LSB first radixsort */
	for(int64_t i = lower_digit; i < higher_digit; i++) {
		/* set digit and pointers */
		for(int64_t j = 0; j < nt; j++) {
			th[j].digit = i;
			th[j].src = src;
			th[j].dst = dst;
		}

		debug("start loop");
		for(int64_t j = 0; j < len; j++) {
			debug("%llu, ", ((elem_t *)src)[j]);
		}

		/* count occ */
		ptask_parallel(pt, count_occ, NULL);

		/* gather occ */
		join(psort_gather_occ_, SUFFIX)(th);

		/* scatter */
		ptask_parallel(pt, scatter, NULL);

		/* flush */
		join(psort_flush_, SUFFIX)(th);

		debug("end loop");
		for(int64_t j = 0; j < len; j++) {
			debug("%llu, ", ((elem_t *)dst)[j]);
		}

		/* swap */
		void *tmp = src; src = dst; dst = tmp;
	}

	/* copyback */
	if((higher_digit - lower_digit) & 0x01) {
		for(int64_t j = 0; j < nt; j++) {
			th[j].src = src;
			th[j].dst = dst;
		}
		ptask_parallel(pt, copyback, NULL);
	}

	/* cleanup ptask object */
	ptask_clean(pt);

	/* cleanup working memory */
	free(ptr);
	return;
}

/* unittests */
unittest_config(
	.name = "psort_radix_intl",
	.depends_on = { "ptask" }
);

/* small integer, single thread */
unittest()
{
	uint64_t raw[] =    { 1, 0, 2, 1, 0, 2, 0, 0, 1, 1 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(int64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 0, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* small integer, 4-thread */
unittest()
{
	uint64_t raw[] =    { 1, 0, 2, 1, 0, 2, 0, 0, 1, 1 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(int64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 4, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* middle integer, single thread */
unittest()
{
	uint64_t raw[] =    { 1000, 0, 2000, 1000, 0, 2000, 0, 0, 1000, 1000 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1000, 1000, 1000, 1000, 2000, 2000 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(int64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 0, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* middle integer, 4-thread */
unittest()
{
	uint64_t raw[] =    { 1000, 0, 2000, 1000, 0, 2000, 0, 0, 1000, 1000 };
	uint64_t sorted[] = { 0, 0, 0, 0, 1000, 1000, 1000, 1000, 2000, 2000 };

	elem_t *arr = (elem_t *)aligned_malloc(sizeof(elem_t) * 10, 16);
	for(int64_t i = 0; i < 10; i++) {
		wr(arr + i, p(raw[i]));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, 10, 4, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < 10; i++) {
		assert(e(rd(arr + i)) == sorted[i],
			"%llu, %llu", e(rd(arr + i)), (uint64_t)sorted[i]);
	}

	/* cleanup */
	free(arr);
}

/* inverse, long array, single thread */
unittest()
{
	int64_t const len = 10000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(int64_t i = 0; i < len; i++) {
		wr(arr + i, p((len - 1) - i));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, len, 0, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < len; i++) {
		assert(e(rd(arr + i)) == i, "%llu", e(rd(arr + i)));
	}

	/* cleanup */
	free(arr);
	return;
}

/* inverse, long array, 4-thread */
unittest()
{
	int64_t const len = 10000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(int64_t i = 0; i < len; i++) {
		wr(arr + i, p((len - 1) - i));
	}

	/* sort */
	join(psort_partialsort_parallel_, SUFFIX)(
		arr, len, 4, 0, sizeof(elem_t));

	/* check */
	for(int64_t i = 0; i < len; i++) {
		assert(e(rd(arr + i)) == i, "%llu", e(rd(arr + i)));
	}

	/* cleanup */
	free(arr);
	return;
}
#if 0
/* benchmark */
#include <sys/time.h>
unittest()
{
	int64_t const len = 200000000;
	elem_t *arr = (elem_t *)aligned_malloc(
		sizeof(elem_t) * len,
		16);

	for(int64_t i = 0; i < 5; i++) {
		/* init array */
		for(int64_t j = 0; j < len; j++) {
			wr(arr + j, p((len - 1) - j));
		}

		struct timeval ts, te;
		gettimeofday(&ts, NULL);

		/* sort */
		join(psort_partialsort_parallel_, SUFFIX)(
			arr, len, 4, 0, sizeof(elem_t));

		gettimeofday(&te, NULL);

		/* check */
		for(int64_t j = 1; j < len; j++) {
			assert(e(rd(arr + j - 1)) <= e(rd(arr + j)),
				"%llu, %llu",
				e(rd(arr + j - 1)), e(rd(arr + j)));
		}

		fprintf(stderr, "%lu us\n",
			(te.tv_sec - ts.tv_sec) * 1000000 + (te.tv_usec - ts.tv_usec));
	}

	/* cleanup */
	free(arr);
	return;
}
#endif

/* cleanup macros */
#undef WCR_BUF_ELEM_COUNT
#undef elem_t
#undef SUFFIX
#undef rd
#undef wr
#undef ex
#undef p
#undef e
#undef UNITTEST_UNIQUE_ID

/**
 * end of psort_radix_intl.c
 */
