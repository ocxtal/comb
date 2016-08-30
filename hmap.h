
/**
 * @file hmap.h
 *
 * @brief string to object hashmap
 */
#ifndef _HMAP_H_INCLUDED
#define _HMAP_H_INCLUDED

#include <stdint.h>

/**
 * @type hmap_header_t
 * @brief object must have a hmap_header_t field at the head.
 */
struct hmap_header_s {
	uint64_t reserved;
};
typedef struct hmap_header_s hmap_header_t;

/**
 * @struct hmap_params_s
 */
struct hmap_params_s {
	uint64_t hmap_size;
	void *lmm;
};
typedef struct hmap_params_s hmap_params_t;
#define HMAP_PARAMS(...)		( &((struct hmap_params_s const){ __VA_ARGS__ }) )

/**
 * @type hmap_t
 */
typedef struct hmap_s hmap_t;

/**
 * @struct hmap_key_s
 * @brief return value container for get_key function
 */
struct hmap_key_s {
	char const *ptr;
	uint32_t len;
};
typedef struct hmap_key_s hmap_key_t;

/**
 * @fn hmap_init
 */
hmap_t *hmap_init(
	uint64_t object_size,
	hmap_params_t const *params);

/**
 * @fn hmap_clean
 */
void hmap_clean(
	hmap_t *hmap);

/**
 * @fn hmap_flush
 */
void hmap_flush(
	hmap_t *hmap);

/**
 * @fn hmap_get_id
 * @brief returns index in the object array
 */
uint32_t hmap_get_id(
	hmap_t *hmap,
	char const *str,
	uint32_t len);

/**
 * @fn hmap_get_key
 */
struct hmap_key_s hmap_get_key(
	hmap_t *hmap,
	uint32_t id);

/**
 * @fn hmap_get_object
 */
void *hmap_get_object(
	hmap_t *hmap,
	uint32_t id);

/**
 * @fn hmap_get_count
 */
uint32_t hmap_get_count(
	hmap_t *hmap);

#endif /* _HMAP_H_INCLUDED */
/**
 * end of hmap.h
 */
