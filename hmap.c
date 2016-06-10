
/**
 * @file hmap.c
 *
 * @brief string to object hashmap
 */

#define UNITTEST_UNIQUE_ID			55
#include "unittest.h"

#include <string.h>
#include <stdint.h>
#include "hmap.h"
#include "lmm.h"
#include "log.h"
#include "sassert.h"


/* constants */
#define HMAP_DEFAULT_HASH_SIZE		( 128 )

/* inline directive */
#define _force_inline				inline

/* roundup */
// #define _roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )
#define _roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )

/* max, min */
#define MAX2(x, y)					( (x) < (y) ? (y) : (x) )
#define MIN2(x, y)					( (x) > (y) ? (y) : (x) )


/**
 * MurmurHash3 string hashing function,
 * extracted from https://github.com/aappleby/smhasher
 * modified to make the functions static and to return
 * hashed value directly.
 */
static _force_inline
uint32_t rotl32(uint32_t x,int8_t r)
{
	return((x << r) | (x >> (32 - r)));
}

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

static _force_inline
uint32_t getblock32(const uint32_t *p, int i)
{
	return(p[i]);
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static _force_inline
uint32_t fmix32(uint32_t h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return(h);
}

//-----------------------------------------------------------------------------
static _force_inline
uint32_t MurmurHash3_x86_32(
	const void *key,
	int32_t len,
	uint32_t seed)
{
	const uint8_t * data = (const uint8_t*)key;
	const int nblocks = len / 4;

	uint32_t h1 = seed;

	const uint32_t c1 = 0xcc9e2d51;
	const uint32_t c2 = 0x1b873593;

	//----------
	// body

	const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

	for(int i = -nblocks; i; i++)
	{
		uint32_t k1 = getblock32(blocks,i);

		k1 *= c1;
		k1 = rotl32(k1,15);
		k1 *= c2;

		h1 ^= k1;
		h1 = rotl32(h1,13); 
		h1 = h1*5+0xe6546b64;
	}

	//----------
	// tail

	const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

	uint32_t k1 = 0;

	switch(len & 3)
	{
	case 3: k1 ^= tail[2] << 16;
	case 2: k1 ^= tail[1] << 8;
	case 1: k1 ^= tail[0];
		k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
	};

	//----------
	// finalization

	h1 ^= len;

	h1 = fmix32(h1);

	return(h1);
}
/* end of MurmurHash3.cpp */



/**
 * MurmurHash3 wrapper functions
 */
static _force_inline
uint32_t hash_string(
	char const *str,
	int32_t len)
{
	return(MurmurHash3_x86_32(
		(void const *)str,
		(int)len,
		0xcafebabe));
}
static _force_inline
uint32_t hash_uint32(uint32_t val)
{
	return(MurmurHash3_x86_32(
		(void const *)&val,
		4,
		val));
}

/**
 * @struct hmap_header_intl_s
 */
struct hmap_header_intl_s {
	uint64_t key_base 	: 48;
	uint64_t key_len	: 16;
};
_static_assert(sizeof(struct hmap_header_intl_s) == sizeof(struct hmap_header_s));

/**
 * @struct hmap_pair_s
 */
struct hmap_pair_s {
	uint32_t id;
	uint32_t hash_val;
};
_static_assert(sizeof(struct hmap_pair_s) == 8);

/**
 * @struct hmap_s
 */
struct hmap_s {
	lmm_t *lmm;
	uint32_t mask;
	uint32_t object_size;
	lmm_kvec_t(uint8_t) key_arr;
	lmm_kvec_t(uint8_t) object_arr;
	uint32_t next_id;
	uint32_t reserved;
	struct hmap_pair_s *table;
};

/**
 * @fn hmap_init
 */
hmap_t *hmap_init(
	uint64_t object_size,
	hmap_params_t const *params)
{
	struct hmap_params_s const default_params = {
		.hmap_size = HMAP_DEFAULT_HASH_SIZE,
		.lmm = NULL
	};
	params = (params == NULL) ? &default_params : params;

	/* size must be power of 2 */
	uint64_t hmap_size = (params->hmap_size == 0)
		? HMAP_DEFAULT_HASH_SIZE
		: params->hmap_size;
	if((hmap_size & (hmap_size - 1)) != 0) {
		return(NULL);
	}

	/* malloc mem */
	lmm_t *lmm = (lmm_t *)params->lmm;
	struct hmap_s *hmap = lmm_malloc(lmm, sizeof(struct hmap_s));
	struct hmap_pair_s *table = lmm_malloc(lmm, sizeof(struct hmap_pair_s) * hmap_size);
	if(hmap == NULL || table == NULL) {
		goto _hmap_init_error_handler;
	}

	/* init context */
	hmap->lmm = lmm;
	hmap->mask = (uint32_t)hmap_size - 1;
	hmap->object_size = _roundup(object_size, 16);
	hmap->next_id = 0;
	hmap->table = table;
	lmm_kv_init(lmm, hmap->key_arr);
	lmm_kv_init(lmm, hmap->object_arr);

	/* init hashmap with invalid mark */
	memset(hmap->table, 0xff, sizeof(struct hmap_pair_s) * hmap_size);
	return((hmap_t *)hmap);

_hmap_init_error_handler:;
	lmm_free(lmm, hmap); hmap = NULL;
	lmm_free(lmm, table); table = NULL;
	return(NULL);
}

/**
 * @fn hmap_clean
 */
void hmap_clean(
	hmap_t *_hmap)
{
	struct hmap_s *hmap = (struct hmap_s *)_hmap;

	if(hmap != NULL) {
		lmm_kv_destroy(hmap->lmm, hmap->key_arr);
		lmm_kv_destroy(hmap->lmm, hmap->object_arr);
		lmm_free(hmap->lmm, hmap->table); hmap->table = NULL;
		lmm_free(hmap->lmm, hmap); hmap = NULL;
	}
	return;
}

/**
 * @fn hmap_flush
 */
void hmap_flush(
	hmap_t *_hmap)
{
	struct hmap_s *hmap = (struct hmap_s *)_hmap;

	if(hmap != NULL) {
		lmm_kv_clear(hmap->lmm, hmap->key_arr);
		lmm_kv_clear(hmap->lmm, hmap->object_arr);

		hmap->next_id = 0;
		memset(hmap->table, 0xff, sizeof(struct hmap_pair_s) * (hmap->mask + 1));
	}
	return;
}

/**
 * @fn hmap_object_get_ptr
 */
static _force_inline
struct hmap_header_intl_s *hmap_object_get_ptr(
	struct hmap_s *hmap,
	uint32_t id)
{
	return((struct hmap_header_intl_s *)(lmm_kv_ptr(hmap->object_arr) + (uint64_t)id * hmap->object_size));
}

/**
 * @fn hmap_object_get_key
 */
static _force_inline
struct hmap_key_s hmap_object_get_key(
	struct hmap_s *hmap,
	uint32_t id)
{
	struct hmap_header_intl_s *obj = hmap_object_get_ptr(hmap, id);
	return((struct hmap_key_s){
		.str = (char const *)lmm_kv_ptr(hmap->key_arr) + obj->key_base,
		.len = obj->key_len
	});
}

/**
 * @fn hmap_get_key
 */
struct hmap_key_s hmap_get_key(
	hmap_t *_hmap,
	uint32_t id)
{
	return(hmap_object_get_key((struct hmap_s *)_hmap, id));
}

/**
 * @fn hmap_expand
 */
static _force_inline
void hmap_expand(
	struct hmap_s *hmap)
{
	uint32_t prev_mask = hmap->mask;
	uint32_t prev_size = prev_mask + 1;
	uint32_t size = 2 * prev_size;
	uint32_t mask = size - 1;

	/* realloc */
	hmap->table = (struct hmap_pair_s *)lmm_realloc(hmap->lmm, hmap->table,
		sizeof(struct hmap_pair_s) * (uint64_t)size);

	/* update context */
	hmap->mask = mask;

	/* init expanded area with invalid mark */
	memset(&hmap->table[prev_size], 0xff, sizeof(struct hmap_pair_s) * prev_size);

	/* rehash */
	for(int64_t i = 0; i < prev_size; i++) {
		/* check id */
		uint32_t const invalid_id = (uint32_t)-1;
		uint32_t id = hmap->table[i].id;
		if(id == invalid_id) { continue; }

		/* check if move is needed */
		uint32_t base_hash_val = hmap->table[i].hash_val;
		if((base_hash_val & mask) == i) { continue; }

		/* move */
		uint32_t hash_val = base_hash_val;
		while(hmap->table[hmap->mask & hash_val].id != invalid_id) {
			hash_val = hash_uint32(hash_val);
		}
		hmap->table[i] = (struct hmap_pair_s){ (uint32_t)-1, (uint32_t)-1 };
		hmap->table[hmap->mask & hash_val] = (struct hmap_pair_s){
			.id = id,
			.hash_val = base_hash_val
		};
		debug("move id(%u) form %lld to %u", id, i, hmap->mask & hash_val);
	}
	debug("expanded, mask(%u)", hmap->mask);
	return;
}

/**
 * @fn hmap_get_id
 */
uint32_t hmap_get_id(
	hmap_t *_hmap,
	char const *str,
	int32_t len)
{
	struct hmap_s *hmap = (struct hmap_s *)_hmap;

	uint32_t const invalid_id = (uint32_t)-1;
	uint32_t id = invalid_id;
	uint32_t tmp_id = invalid_id;

	uint32_t base_hash_val = hash_string(str, len);
	uint32_t hash_val = base_hash_val;

	while((tmp_id = hmap->table[hmap->mask & hash_val].id) != invalid_id) {
		struct hmap_key_s ex_key = hmap_object_get_key(hmap, tmp_id);
		debug("ex_key at %u: str(%p), len(%d)", tmp_id, ex_key.str, ex_key.len);

		/* compare string */
		if(ex_key.len == len && strncmp(ex_key.str, str, MIN2(ex_key.len, len) + 1) == 0) {
			/* matched with existing string in the section array */
			id = tmp_id; break;
		}

		debug("collision at %u, key(%s), ex_key(%s), hash_val(%x)",
			hmap->mask & hash_val, str, ex_key.str, hash_uint32(hash_val));

		/* not matched, rehash */
		hash_val = hash_uint32(hash_val);
	}

	if(id == invalid_id) {
		/* rehash if occupancy exceeds 0.5 */
		if(hmap->next_id > (hmap->mask + 1) / 2) {
			debug("check size next_id(%u), size(%u)",
				hmap->next_id, (hmap->mask + 1) / 2);
			hmap_expand(hmap);
		}

		/* add (hash_val, id) pair to table */
		debug("id(%u), mask(%x), hash_val(%x), id(%u), base_hash_val(%x)",
			hmap->mask & hash_val, hmap->mask, hash_val, hmap->next_id, base_hash_val);
		hmap->table[hmap->mask & hash_val] = (struct hmap_pair_s){
			.id = (id = hmap->next_id++),
			.hash_val = base_hash_val
		};

		/* reserve working area */
		uint8_t tmp[hmap->object_size];
		// memset(tmp, 0, hmap->object_size);
		struct hmap_header_intl_s *h = (struct hmap_header_intl_s *)tmp;

		/* push key string to key_arr */
		h->key_len = len;
		h->key_base = lmm_kv_size(hmap->key_arr);
		lmm_kv_pushm(hmap->lmm, hmap->key_arr, str, len);
		lmm_kv_push(hmap->lmm, hmap->key_arr, '\0');

		/* add object to object array */
		lmm_kv_pushm(hmap->lmm, hmap->object_arr, tmp, hmap->object_size);
	}
	return(id);
}

/**
 * @fn hmap_get_object
 */
void *hmap_get_object(
	hmap_t *hmap,
	uint32_t id)
{
	return((void *)hmap_object_get_ptr((struct hmap_s *)hmap, id));
}

/**
 * @fn hmap_get_count
 */
uint32_t hmap_get_count(
	hmap_t *_hmap)
{
	struct hmap_s *hmap = (struct hmap_s *)_hmap;
	return(hmap->next_id);
}


/* unittests */
unittest_config(
	.name = "hmap",
);

#define make_string(x) ({ \
	char buf[128]; \
	sprintf(buf, "key-%" PRId64 "", (int64_t)(x)); \
	buf; \
})
#define make_args(x)	make_string(x), strlen(make_string(x))

#define UNITTEST_KEY_COUNT			( 32768 )

/* create context */
unittest()
{
	/* valid hash size */
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t), NULL);
	assert(hmap != NULL, "%p", hmap);
	hmap_clean(hmap);
}

/* with params */
unittest()
{
	/* valid hash size */
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t),
		HMAP_PARAMS(.hmap_size = 128));

	assert(hmap != NULL, "%p", hmap);
	hmap_clean(hmap);

	/* invalid hash size */
	hmap = hmap_init(16, HMAP_PARAMS(.hmap_size = 0));
	assert(hmap != NULL, "%p", hmap);
	hmap_clean(hmap);

	hmap = hmap_init(16, HMAP_PARAMS(.hmap_size = 127));
	assert(hmap == NULL, "%p", hmap);
	hmap_clean(hmap);
}

/* with lmm object */
unittest()
{
	lmm_t *lmm = lmm_init(NULL, 1024);
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t),
		HMAP_PARAMS(.lmm = lmm));

	assert(hmap != NULL, "%p", hmap);
	hmap_clean(hmap);
	lmm_clean(lmm);
}

/* append and get */
unittest()
{
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t), NULL);

	/* append key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));

		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* get key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);

		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* get key in reverse order */
	for(int64_t i = UNITTEST_KEY_COUNT - 1; i >= 0; i--) {
		struct hmap_key_s k = hmap_get_key(hmap, i);

		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	hmap_clean(hmap);
}

/* flush */
unittest()
{
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t), NULL);

	/* append key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));

		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* flush the content */
	hmap_flush(hmap);

	/* append again */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));

		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* get key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);

		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	hmap_clean(hmap);
}

/* with lmm object */
unittest()
{
	lmm_t *lmm = lmm_init(NULL, 1024);
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t),
		HMAP_PARAMS(.lmm = lmm));

	/* append key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		
		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* get key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);

		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	/* get key in reverse order */
	for(int64_t i = UNITTEST_KEY_COUNT - 1; i >= 0; i--) {
		struct hmap_key_s k = hmap_get_key(hmap, i);

		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	assert(hmap_get_count(hmap) == UNITTEST_KEY_COUNT, "count(%u)", hmap_get_count(hmap));

	hmap_clean(hmap);
	lmm_clean(lmm);
}

/* with different object size */
unittest()
{
	/* 32 */
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t) + 32, NULL);
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}

	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);
		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	hmap_clean(hmap);

	/* 36 */
	hmap = hmap_init(sizeof(hmap_header_t) + 36, NULL);
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}

	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);
		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	hmap_clean(hmap);

	/* 127 */
	hmap = hmap_init(sizeof(hmap_header_t) + 127, NULL);
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}

	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);
		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	hmap_clean(hmap);
}

/* with lmm */
unittest()
{
	lmm_t *lmm = lmm_init(NULL, 1024 * 1024);

	/* 127 */
	hmap_t *hmap = hmap_init(sizeof(hmap_header_t) + 127,
		HMAP_PARAMS(.lmm = lmm));
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		assert((int64_t)id == i, "i(%lld), id(%lld)", i, id);
	}

	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct hmap_key_s k = hmap_get_key(hmap, i);
		assert(k.len == strlen(make_string(i)), "a(%d), b(%d)", k.len, strlen(make_string(i)));
		assert(strcmp(k.str, make_string(i)) == 0, "a(%s), b(%s)", k.str, make_string(i));
	}
	hmap_clean(hmap);
	lmm_clean(lmm);
}


/* get_object */
unittest()
{
	struct str_cont_s {
		hmap_header_t header;
		char s[128];
	};
	hmap_t *hmap = hmap_init(sizeof(struct str_cont_s), NULL);

	/* append key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		struct str_cont_s *obj = hmap_get_object(hmap, id);

		assert(obj != NULL, "obj(%p)", obj);

		strcpy(obj->s, make_string(i));
	}

	/* get key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct str_cont_s *obj = hmap_get_object(hmap, i);

		assert(obj != NULL, "obj(%p)", obj);
		assert(strcmp(obj->s, make_string(i)) == 0, "%s, %s", obj->s, make_string(i));
	}

	/* get key in reverse order */
	for(int64_t i = UNITTEST_KEY_COUNT - 1; i >= 0; i--) {
		struct str_cont_s *obj = hmap_get_object(hmap, i);

		assert(obj != NULL, "obj(%p)", obj);
		assert(strcmp(obj->s, make_string(i)) == 0, "%s, %s", obj->s, make_string(i));
	}

	hmap_clean(hmap);
}

/* with lmm */
unittest()
{
	lmm_t *lmm = lmm_init(NULL, 1024 * 1024);

	struct str_cont_s {
		hmap_header_t header;
		char s[128];
	};
	hmap_t *hmap = hmap_init(sizeof(struct str_cont_s),
		HMAP_PARAMS(.lmm = lmm));

	/* append key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		uint32_t id = hmap_get_id(hmap, make_args(i));
		struct str_cont_s *obj = hmap_get_object(hmap, id);

		assert(obj != NULL, "obj(%p)", obj);

		strcpy(obj->s, make_string(i));
	}

	/* get key */
	for(int64_t i = 0; i < UNITTEST_KEY_COUNT; i++) {
		struct str_cont_s *obj = hmap_get_object(hmap, i);

		assert(obj != NULL, "obj(%p)", obj);
		assert(strcmp(obj->s, make_string(i)) == 0, "%s, %s", obj->s, make_string(i));
	}

	/* get key in reverse order */
	for(int64_t i = UNITTEST_KEY_COUNT - 1; i >= 0; i--) {
		struct str_cont_s *obj = hmap_get_object(hmap, i);

		assert(obj != NULL, "obj(%p)", obj);
		assert(strcmp(obj->s, make_string(i)) == 0, "%s, %s", obj->s, make_string(i));
	}

	hmap_clean(hmap);
	lmm_clean(lmm);
}


/**
 * end of hmap.c
 */
