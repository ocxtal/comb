
/**
 * @file zf.c
 *
 * @brief zlib-file API compatible I/O wrapper library
 */

#define UNITTEST_UNIQUE_ID			44
#define UNITTEST 					1

#include "unittest.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "kopen.h"
#include "sassert.h"
#include "zf.h"

#ifdef HAVE_Z
#include "zlib.h"
#endif

#ifdef HAVE_BZ2
#include "bzlib.h"
#endif


/* constants */
#define ZF_BUF_SIZE					( 512 * 1024 )		/* 512KB */
#define ZF_UNGETC_MARGIN_SIZE		( 32 )

/* function pointer type aliases */
typedef void *(*zf_dopen_t)(
	int fd,
	char const *mode);
typedef void *(*zf_open_t)(
	char const *path,
	char const *mode);
typedef int (*zf_init_t)(
	void *);
typedef void *(*zf_close_t)(
	void *fp);
typedef size_t (*zf_read_t)(
	void *fp,
	void *ptr,
	size_t len);
typedef size_t (*zf_write_t)(
	void *fp,
	void *ptr,
	size_t len);

/* wrapped functions */

/**
 * @fn zf_init_stdio
 * @brief wrapped setvbuf
 */
int zf_init_stdio(
	void *fp)
{
	return(setvbuf((FILE *)fp, NULL, _IOFBF, ZF_BUF_SIZE));
}

/**
 * @fn fread_wrap
 * @brief wrap fread to zlib compatible
 */
size_t fread_wrap(
	void *fp,
	void *ptr,
	size_t len)
{
	return(fread(ptr, 1, len, (FILE *)fp));
}

/**
 * @fn fwrite_wrap
 * @brief wrap fwrite to zlib compatible
 */
size_t fwrite_wrap(
	void *fp,
	void *ptr,
	size_t len)
{
	return(fwrite(ptr, 1, len, (FILE *)fp));
}

/* zlib-dependent functions */
#ifdef HAVE_Z
/**
 * @fn zf_init_gzip
 * @brief setup function, set buffer size to 512k
 */
int zf_init_gzip(
	void *fp)
{
	#if ZLIB_VERNUM >= 0x1240
		return(gzbuffer((gzFile)fp, ZF_BUF_SIZE));
	#else
		return(0);
	#endif
}
#endif

/**
 * @struct zf_functions_s
 * @brief function container
 */
struct zf_functions_s {
	char const *ext;
	zf_dopen_t dopen;
	zf_open_t open;
	zf_init_t init;
	zf_close_t close;
	zf_read_t read;
	zf_write_t write;	
};

/**
 * @struct zf_intl_s
 * @brief context container
 */
struct zf_intl_s {
	char *path;
	char *mode;
	int fd;
	int eof;		/* == 1 if fp reached EOF, == 2 if curr reached the end of buf */
	void *ko;
	void *fp;		/* one of {FILE * / gzFile / BZFILE *} */
	struct zf_functions_s fn;
	uint8_t *buf;
	int64_t size;
	int64_t curr, end;
	char ungetc_margin[ZF_UNGETC_MARGIN_SIZE];
};
_static_assert(offsetof(struct zf_intl_s, ungetc_margin) == sizeof(struct zf_s));
_static_assert_offset(struct zf_intl_s, path, struct zf_s, path, 0);
_static_assert_offset(struct zf_intl_s, mode, struct zf_s, mode, 0);

/**
 * @val fn_table
 * @brief extension and function table
 */
static
struct zf_functions_s const fn_table[] = {
	/* default */
	{
		.ext = "",
		.dopen = (zf_dopen_t)fdopen,
		.open = (zf_open_t)fopen,
		.init = (zf_init_t)zf_init_stdio,
		.close = (zf_close_t)fclose,
		.read = (zf_read_t)fread_wrap,
		.write = (zf_write_t)fwrite_wrap
	},
	/* gzip */
	{
		.ext = ".gz",
		#ifdef HAVE_Z
		.dopen = (zf_dopen_t)gzdopen,
		.open = (zf_open_t)gzopen,
		.init = (zf_init_t)zf_init_gzip,
		.close = (zf_close_t)gzclose,
		.read = (zf_read_t)gzread,
		.write = (zf_write_t)gzwrite
		#endif
	},
	/* bzip2 */
	{
		.ext = ".bz2",
		#ifdef HAVE_BZ2
		.dopen = (zf_dopen_t)BZ2_bzdopen,
		.open = (zf_open_t)BZ2_bzopen,
		.init = (zf_init_t)NULL,
		.close = (zf_close_t)BZ2_bzclose,
		.read = (zf_read_t)BZ2_bzread,
		.write = (zf_write_t)BZ2_bzwrite
		#endif
	},
	/* other unsupported formats */
	{ .ext = ".lz" },
	{ .ext = ".lzma" },
	{ .ext = ".xz" },
	{ .ext = ".z" }
};

/**
 * @fn zfopen
 * @brief open file, similar to fopen / gzopen,
 * compression format can be explicitly specified adding an extension to `mode', e.g. "w+.bz2".
 */
zf_t *zfopen(
	char const *path,
	char const *mode)
{
	if(path == NULL || mode == NULL) {
		return(NULL);
	}

	/* check length */
	uint64_t path_len = strlen(path);
	uint64_t mode_len = strlen(mode);
	if(path_len == 0 || mode_len == 0) {
		return(NULL);
	}

	char *path_dup = (char *)path;
	char *mode_dup = (char *)mode;

	char const *path_tail = path + path_len;
	char const *mode_tail = mode + mode_len;

	/* determine format */
	struct zf_functions_s const *fn = &fn_table[0];
	for(uint64_t i = 1; i < sizeof(fn_table) / sizeof(struct zf_functions_s); i++) {
		/* skip if ext string is longer than path string */
		if(path_len < strlen(fn_table[i].ext)) { continue; }

		/* check path */
		if(strncmp(path_tail - strlen(fn_table[i].ext),
			fn_table[i].ext,
			strlen(fn_table[i].ext)) == 0) {

			/* hit */
			fn = &fn_table[i];
			path_dup = strdup(path);
			path_dup[strlen(path) - strlen(fn_table[i].ext)] = '\0';
			break;
		}

		/* skip if ext string is longer than mode string */
		if(mode_len < strlen(fn_table[i].ext)) { continue; }

		/* check mode */
		if(strncmp(mode_tail - strlen(fn_table[i].ext),
			fn_table[i].ext,
			strlen(fn_table[i].ext)) == 0) {
			
			/* hit */
			fn = &fn_table[i];
			mode_dup = strdup(mode);
			mode_dup[strlen(mode) - strlen(fn_table[i].ext)] = '\0';
			break;
		}
	}

	/* check if functions are available */
	if(fn->open == NULL) {
		return(NULL);
	}

	/* malloc context */
	struct zf_intl_s *fio = (struct zf_intl_s *)malloc(
		sizeof(struct zf_intl_s) + ZF_BUF_SIZE);
	if(fio == NULL) {
		return(NULL);
	}
	memset(fio, 0, sizeof(struct zf_intl_s));
	fio->buf = (uint8_t *)(fio + 1);
	fio->size = ZF_BUF_SIZE;
	fio->fn = *fn;

	/* open file */
	if(mode[0] == 'r') {
		/* read mode, open file with kopen */
		fio->ko = kopen(path, &fio->fd);
		if(fio->ko == NULL) {
			goto _zfopen_finish;
		}
		fio->fp = fio->fn.dopen(fio->fd, mode);
	} else {
		/* write mode, check if stdout is specified */
		if(strncmp(path, "-", strlen("-")) == 0) {
			fio->path = "-";
			fio->fd = STDOUT_FILENO;
			fio->ko = NULL;
			fio->fp = stdout;
			goto _zfopen_finish;
		}

		/* open file */
		fio->fp = fio->fn.open(path, mode);
		fio->fd = -1;		/* fd is invalid in write mode */
		fio->ko = NULL;		/* ko is also invalid */
		goto _zfopen_finish;
	}

_zfopen_finish:;
	if(fio->fp == NULL) {
		/* something wrong occurred */
		if(fio->ko != NULL) {
			kclose(fio->ko); fio->ko = NULL;
		}
		free(fio); fio = NULL;
		return(NULL);
	}

	/* everything is going right */
	fio->path = (path_dup == path) ? strdup(path) : path_dup;
	fio->mode = (mode_dup == mode) ? strdup(mode) : mode_dup;
	fio->curr = fio->end = 0;
	if(fio->fn.init != NULL) {
		if(fio->fn.init(fio->fp) != 0) {
			zfclose((zf_t *)fio);
			return(NULL);
		}
	}

	return((zf_t *)fio);
}

/**
 * @fn zfclose
 * @brief close file, similar to fclose / gzclose
 */
int zfclose(
	zf_t *fp)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;

	if(fio == NULL) {
		return(1);
	}

	/* flush if write mode */
	if(fio->mode[0] != 'r') {
		fio->fn.write(fio->fp, (void *)fio->buf, fio->curr);
	}

	/* close file */
	if(fio->fp != NULL) {
		fio->fn.close(fio->fp); fio->fp = NULL;
		if(fio->ko != NULL) {
			kclose(fio->ko); fio->ko = NULL;
		}
	}
	free(fio->path); fio->path = NULL;
	free(fio->mode); fio->mode = NULL;
	free(fio); fio = NULL;
	return(0);
}

/**
 * @fn zfread
 * @brief read from file, similar to gzread
 */
size_t zfread(
	zf_t *fp,
	void *_ptr,
	size_t len)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	uint8_t *ptr = (uint8_t *)_ptr;
	size_t copied_size = 0;

	/* check eof */
	if(fio->eof == 2) {
		return(0);
	}

	/* if buffer is not empty, copy from buffer */
	if(fio->curr < fio->end) {
		uint64_t rem_size = fio->end - fio->curr;
		uint64_t buf_copy_size = (rem_size < len) ? rem_size : len;
		memcpy(ptr, (void *)&fio->buf[fio->curr], buf_copy_size);

		/* update */
		len -= buf_copy_size;
		ptr += buf_copy_size;
		fio->curr += buf_copy_size;
		copied_size += buf_copy_size;
	}

	/* if fp already reached EOF */
	if(fio->eof == 1) {
		fio->eof += (fio->curr == fio->end);
		return(copied_size);
	}

	/* issue fread */
	if(len > 0) {
		uint64_t read_size = fio->fn.read(fio->fp, ptr, len);
		fio->eof = 2 * (read_size < len);
		copied_size += read_size;
	}
	return(copied_size);
}

/**
 * @fn zfpeek
 * @brief read len (< 512k) without advancing pointer
 */
size_t zfpeek(
	zf_t *fp,
	void *_ptr,
	size_t len)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	uint8_t *ptr = (uint8_t *)_ptr;
	size_t copied_size = 0;

	/* check eof */
	if(fio->eof == 2) {
		return(0);
	}

	/* copy from buffer */
	if(fio->curr < fio->end) {
		uint64_t rem_size = fio->end - fio->curr;
		uint64_t buf_copy_size = (rem_size < len) ? rem_size : len;
		memcpy(ptr, (void *)&fio->buf[fio->curr], buf_copy_size);
		
		/* update length (without advancing pointer) */
		len -= buf_copy_size;
		ptr += buf_copy_size;
		copied_size += buf_copy_size;
	}

	if(len > 0) {
		/* move existing elements to the head of the buffer */
		if(fio->curr != 0 && fio->curr < fio->end) {
			memmove((void *)fio->buf, (void *)&fio->buf[fio->curr], fio->end - fio->curr);
			fio->end -= fio->curr;
			fio->curr = 0;
		}

		/* read */
		uint64_t read_size = fio->fn.read(fio->fp, (void *)&fio->buf[fio->end], fio->size - fio->end);

		/* copy */
		uint64_t copy_size = (read_size < len) ? read_size : len;
		memcpy(ptr, (void *)&fio->buf[fio->end], copy_size);

		/* update status */
		copied_size += copy_size;
		fio->eof = (read_size < (uint64_t)(fio->size - fio->end)) + (copied_size == 0);
		fio->end += read_size;
	}
	return(copied_size);
}

/**
 * @fn zfwrite
 * @brief write to file, similar to gzwrite
 */
size_t zfwrite(
	zf_t *fp,
	void *ptr,
	size_t len)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	return(fio->fn.write(fio->fp, ptr, len));
}

/**
 * @fn zfgetc
 */
int zfgetc(
	zf_t *fp)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;

	/* if the pointer reached the end, refill the buffer */
	if(fio->curr >= fio->end) {
		fio->curr = 0;
		fio->end = (fio->eof == 0)
			? fio->fn.read(fio->fp, fio->buf, fio->size)
			: 0;
		fio->eof = (fio->end < fio->size) + (fio->end == 0);
	}
	if(fio->eof == 2) {
		return(EOF);
	}
	return((int)fio->buf[fio->curr++]);
}

/**
 * @fn zfungetc
 */
int zfungetc(
	zf_t *fp,
	int c)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	if(fio->curr > -ZF_UNGETC_MARGIN_SIZE) {
		return(fio->buf[--fio->curr] = c);
	} else {
		return(-1);
	}
}

/**
 * @fn zfeof
 */
int zfeof(
	zf_t *fp)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	return((int)(fio->eof == 2));
}

/**
 * @fn zfputc
 */
int zfputc(
	zf_t *fp,
	int c)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	fio->buf[fio->curr++] = (uint8_t)c;
	
	/* flush if buffer is full */
	if(fio->curr == fio->size) {
		fio->curr = 0;
		uint64_t written = fio->fn.write(fio->fp, fio->buf, fio->size);
		if((int64_t)written != fio->size) {
			return(-1);
		}
	}
	return(c);
}

/**
 * @fn zfputs
 */
int zfputs(
	zf_t *fp,
	char const *s)
{
	while(*s != '\0') {
		zfputc(fp, (int)*s++);
	}
	zfputc(fp, '\n');
	return(0);
}

/**
 * @fn zfprintf
 */
int zfprintf(
	zf_t *fp,
	char const *format,
	...)
{
	struct zf_intl_s *fio = (struct zf_intl_s *)fp;
	va_list l;
	va_start(l, format);

	/* flush */
	if(fio->curr != 0) {
		uint64_t flush = fio->fn.write(fio->fp, fio->buf, fio->curr);
		/* something is wrong */
		if((int64_t)flush != fio->curr) {
			va_end(l);
			return(0);
		}
	}

	/* fprintf */
	uint64_t size = vsprintf((char *)fio->buf, format, l);

	/* flush */
	uint64_t written = fio->fn.write(fio->fp, fio->buf, size);
	fio->curr = size - written;

	va_end(l);
	return((int)written);
}

/* unittests */
#include <time.h>

/**
 * @fn init_rand
 */
void *init_rand(
	void *params)
{
	srand(time(NULL));
	return(NULL);
}
void clean_rand(
	void *gctx)
{
	return;		/* nothing to do */
}

unittest_config(
	.name = "zf",
	.init = init_rand,
	.clean = clean_rand
);

/**
 * @val ascii_table
 */
static
char const *ascii_table =
	" !\"#$%&\'()*+,-./"
	"0123456789:;<=>?"
	"@ABCDEFGHIJKLMNO"
	"PQRSTUVWXYZ[\\]^_"
	"`abcdefghijklmno"
	"pqrstuvwxyz{|}~\n";

/**
 * @fn make_random_array
 */
static inline
void *make_random_array(
	void *params)
{
	uint64_t len = (uint64_t)params;
	char *arr = (char *)malloc(len);
	for(uint64_t i = 0; i < len; i++) {
		arr[i] = ascii_table[rand() % strlen(ascii_table)];
	}
	return((void *)arr);
}

#define with(_len) \
	.init = make_random_array, \
	.clean = free, \
	.params = (void *)((uint64_t)(_len))

#define omajinai() \
	char *arr = (char *)ctx;

/* return NULL when failed to open */
unittest()
{
	char const *filename = "test.txt";

	remove(filename);
	zf_t *zf = zfopen(filename, "r");
	assert(zf == NULL);
}

/* redirect to stdout */
unittest()
{
	zf_t *zf = zfopen("-", "w");
	assert(zf != NULL);
	zfclose(zf);
}

/* test zfopen / zfclose with len 100M */
#define TEST_ARR_LEN 		1000000
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmp.txt", "w");
	assert(wfp != NULL, "%p", wfp);

	size_t written = zfwrite(wfp, arr, TEST_ARR_LEN);
	assert(written == TEST_ARR_LEN, "%llu", written);

	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmp.txt", "r");
	assert(rfp != NULL, "%p", rfp);
	assert(strcmp(rfp->path, "tmp.txt") == 0, "%s", rfp->path);

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	size_t read = zfread(rfp, rarr, TEST_ARR_LEN);
	assert(read == TEST_ARR_LEN, "%llu", read);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt");
}

/* getc / putc */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write with zfputc */
	zf_t *wfp = zfopen("tmp.txt", "w");

	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		zfputc(wfp, arr[i]);
	}

	zfclose(wfp);

	/* read with zfgetc */
	zf_t *rfp = zfopen("tmp.txt", "r");

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		rarr[i] = zfgetc(rfp);
	}

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt");
}

/* peek */
unittest(with(100000))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmp.txt", "w");
	zfwrite(wfp, arr, 100000);
	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmp.txt", "r");

	char *rarr1 = (char *)malloc(TEST_ARR_LEN);
	char *rarr2 = (char *)malloc(TEST_ARR_LEN);
	
	/* read the former half */
	size_t read1 = zfpeek(rfp, rarr1, 50000);
	assert(zfeof(rfp) == 0, "%d", zfeof(rfp));
	assert(read1 == 50000, "%llu", read1);


	size_t read2 = zfread(rfp, rarr2, 50000);
	assert(zfeof(rfp) == 0, "%d", zfeof(rfp));
	assert(read2 == 50000, "%llu", read2);

	assert(memcmp(arr, rarr1, 50000) == 0);
	assert(memcmp(arr, rarr2, 50000) == 0);

	/* read the latter half */
	read1 = zfpeek(rfp, rarr1, 50000);
	assert(read1 == 50000, "%llu", read1);

	read2 = zfread(rfp, rarr2, 50000);
	assert(read2 == 50000, "%llu", read2);

	assert(memcmp(&arr[50000], rarr1, 50000) == 0);
	assert(memcmp(&arr[50000], rarr2, 50000) == 0);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* cleanup */
	free(rarr1);
	free(rarr2);
	remove("tmp.txt");
}

/* ungetc */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmp.txt", "w");
	zfwrite(wfp, arr, TEST_ARR_LEN);
	zfclose(wfp);

	/* push pop */
	zf_t *rfp = zfopen("tmp.txt", "r");
	int save = 'a';
	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		zfungetc(rfp, save);
		int c = zfgetc(rfp);
		assert(c == save);
		save = zfgetc(rfp);
	}

	save = zfgetc(rfp);
	assert(save == EOF, "%d", save);

	zfclose(rfp);
	remove("tmp.txt");
}

/* zlib-dependent tests */
#ifdef HAVE_Z
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmp.txt.gz", "w");
	assert(wfp != NULL, "%p", wfp);

	size_t written = zfwrite(wfp, arr, TEST_ARR_LEN);
	assert(written == TEST_ARR_LEN, "%llu", written);

	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmp.txt.gz", "r");
	assert(rfp != NULL, "%p", rfp);
	assert(strcmp(rfp->path, "tmp.txt") == 0, "%s", rfp->path);

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	size_t read = zfread(rfp, rarr, TEST_ARR_LEN);
	assert(read == TEST_ARR_LEN, "%llu", read);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt.gz");
}

/* specify compression format with mode flag */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmpfile", "w.gz");
	assert(wfp != NULL, "%p", wfp);
	assert(strcmp(wfp->mode, "w") == 0, "%s", wfp->mode);

	size_t written = zfwrite(wfp, arr, TEST_ARR_LEN);
	assert(written == TEST_ARR_LEN, "%llu", written);

	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmpfile", "r.gz");
	assert(rfp != NULL, "%p", rfp);
	assert(strcmp(rfp->path, "tmpfile") == 0, "%s", rfp->path);
	assert(strcmp(rfp->mode, "r") == 0, "%s", rfp->mode);

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	size_t read = zfread(rfp, rarr, TEST_ARR_LEN);
	assert(read == TEST_ARR_LEN, "%llu", read);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmpfile");
}

/* getc / putc */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write with zfputc */
	zf_t *wfp = zfopen("tmp.txt.gz", "w");

	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		zfputc(wfp, arr[i]);
	}

	zfclose(wfp);

	/* read with zfgetc */
	zf_t *rfp = zfopen("tmp.txt.gz", "r");

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		rarr[i] = zfgetc(rfp);
	}

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt.gz");
}
#endif /* HAVE_Z */

/* bzip2-dependent tests */
#ifdef HAVE_BZ2
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmp.txt.bz2", "w");
	assert(wfp != NULL, "%p", wfp);

	size_t written = zfwrite(wfp, arr, TEST_ARR_LEN);
	assert(written == TEST_ARR_LEN, "%llu", written);

	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmp.txt.bz2", "r");
	assert(rfp != NULL, "%p", rfp);
	assert(strcmp(rfp->path, "tmp.txt") == 0, "%s", rfp->path);

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	size_t read = zfread(rfp, rarr, TEST_ARR_LEN);
	assert(read == TEST_ARR_LEN, "%llu", read);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt.bz2");
}

/* specify compression format with mode flag */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write */
	zf_t *wfp = zfopen("tmpfile", "w.bz2");
	assert(wfp != NULL, "%p", wfp);
	assert(strcmp(wfp->mode, "w") == 0, "%s", wfp->mode);

	size_t written = zfwrite(wfp, arr, TEST_ARR_LEN);
	assert(written == TEST_ARR_LEN, "%llu", written);

	zfclose(wfp);

	/* read */
	zf_t *rfp = zfopen("tmpfile", "r.bz2");
	assert(rfp != NULL, "%p", rfp);
	assert(strcmp(rfp->path, "tmpfile") == 0, "%s", rfp->path);
	assert(strcmp(rfp->mode, "r") == 0, "%s", rfp->mode);

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	size_t read = zfread(rfp, rarr, TEST_ARR_LEN);
	assert(read == TEST_ARR_LEN, "%llu", read);

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmpfile");
}

/* getc / putc */
unittest(with(TEST_ARR_LEN))
{
	omajinai();

	/* write with zfputc */
	zf_t *wfp = zfopen("tmp.txt.bz2", "w");

	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		zfputc(wfp, arr[i]);
	}

	zfclose(wfp);

	/* read with zfgetc */
	zf_t *rfp = zfopen("tmp.txt.bz2", "r");

	char *rarr = (char *)malloc(TEST_ARR_LEN);
	for(int64_t i = 0; i < TEST_ARR_LEN; i++) {
		rarr[i] = zfgetc(rfp);
	}

	/* EOF */
	assert(zfgetc(rfp) == EOF, "%d", zfgetc(rfp));
	assert(zfeof(rfp) != 0, "%d", zfeof(rfp));

	zfclose(rfp);

	/* compare */
	assert(memcmp(arr, rarr, TEST_ARR_LEN) == 0);

	/* cleanup */
	free(rarr);
	remove("tmp.txt.bz2");
}
#endif /* HAVE_BZ2 */

/**
 * end of zf.c
 */
