
/**
 * @file zf.h
 *
 * @brief zlib file API compatible I/O wrapper library
 *
 * @author Hajime Suzuki
 * @date 2016/3/30
 * @license MIT
 */
#ifndef _ZF_H_INCLUDED
#define _ZF_H_INCLUDED

#include <stdint.h>

/* type aliases */
struct zf_s {
	char const *path;
	char const *mode;
	int reserved1[2];
	void *reserved2[10];
	int64_t reserved3[3];

};
typedef struct zf_s zf_t;


/**
 * @fn zfopen
 * @brief open file, similar to fopen / gzopen,
 * compression format can be explicitly specified adding an extension to `mode', e.g. "w+.bz2".
 */
zf_t *zfopen(
	char const *path,
	char const *mode);

/**
 * @fn zfclose
 * @brief close file, similar to fclose / gzclose
 */
int zfclose(
	zf_t *zf);

/**
 * @fn zfread
 * @brief read from file, similar to gzread
 */
size_t zfread(
	zf_t *zf,
	void *ptr,
	size_t len);

/**
 * @fn zfpeek
 * @brief read len (< 512k) without advancing pointer
 */
size_t zfpeek(
	zf_t *zf,
	void *ptr,
	size_t len);

/**
 * @fn zfwrite
 * @brief write to file, similar to gzwrite
 */
size_t zfwrite(
	zf_t *zf,
	void *ptr,
	size_t len);

/**
 * @fn zfgetc
 */
int zfgetc(
	zf_t *zf);

/**
 * @fn zfungetc
 */
int zfungetc(
	zf_t *zf,
	int c);

/**
 * @fn zfeof
 */
int zfeof(
	zf_t *zf);

/**
 * @fn zfputc
 */
int zfputc(
	zf_t *zf,
	int c);

/**
 * @fn zfputs
 */
int zfputs(
	zf_t *zf,
	char const *s);

/**
 * @fn zfprintf
 */
int zfprintf(
	zf_t *zf,
	char const *format,
	...);

#endif /* _ZF_H_INCLUDED */
/**
 * end of zf.h
 */
