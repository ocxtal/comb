
/**
 * @file psort.h
 *
 * @brief parallel integer sort (radixsort) library
 *
 * @author Hajime Suzuki
 * @date 2016/3/20
 * @license MIT
 */
#ifndef _PSORT_H_INCLUDED
#define _PSORT_H_INCLUDED

/**
 * @fn psort_full
 * @brief integer sort
 */
int psort_full(
	void *ptr,
	int64_t len,
	int64_t elem_size,
	int64_t num_threads);

/**
 * @fn psort_half
 * @brief key sort (sort the lower half of the element)
 */
int psort_half(
	void *ptr,
	int64_t len,
	int64_t elem_size,
	int64_t num_threads);

/**
 * @fn psort_partial
 */
int psort_partial(
	void *arr,
	int64_t len,
	int64_t elem_size,
	int64_t num_threads,
	int64_t from,
	int64_t to);

#endif /* _PSORT_H_INCLUDED */
/**
 * end of psort.h
 */
