
/**
 * @file mem.h
 */
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

/**
 * @fn mem_dump_file
 */
static inline
char *mem_dump_file(
	FILE *fp)
{
	uint64_t used = 0, size = 1024;
	char *buf = malloc(size);
	if(buf == NULL) { return(NULL); }

	#define push(x) { \
		if(used == size) { \
			buf = realloc(buf, size *= 2); \
			if(buf == NULL) { return(NULL); } \
		} \
		buf[used++] = (x); \
	}

	int c;
	while((c = getc(fp)) != EOF) { push(c); }

	/* push terminator (double '\0') */
	push('\0');
	push('\0');

	#undef push
	return(buf);
}


#ifdef __linux__
/* dump /proc/meminfo, sum MemFree and Cached */

/**
 * @fn mem_dump_meminfo
 */
static inline
char *mem_dump_vm_stat(
	void)
{
	FILE *fp = NULL;
	char *res = NULL;

	/* open */
	if((fp = fopen("/proc/meminfo", "r")) == NULL) {
		goto _mem_dump_meminfo_error_handler;
	}

	/* dump */
	if((res = mem_dump_file(fp)) == NULL) {
		goto _mem_dump_meminfo_error_handler;
	}

	/* close file */
	if(fclose(fp) != 0) {
		goto _mem_dump_meminfo_error_handler;
	}
	return(res);

_mem_dump_meminfo_error_handler:
	if(fp != NULL) { fclose(fp); }
	if(res != NULL) { free(res); }
	return(NULL);
}

/**
 * @fn mem_estimate_free_size
 */
static inline
uint64_t mem_estimate_free_size(
	void)
{
	char *s = mem_dump_vm_stat();
	if(s == NULL) {
		return(0);
	}

	char const *labels[] = {
		"MemFree:",
		"Cached:",
		NULL
	};

	/* supposing s is terminated with double '\0's */
	uint64_t mem = 0;
	char *p = s;
	while(*p != '\0') {

		char *t = p;
		while(*t != '\n' && *t != '\0') { t++; }

		for(char const **l = labels; *l != NULL; l++) {
			if(strncmp(p, *l, strlen(*l)) == 0) {
				t[-strlen(" kB")] = '\0';
				mem += atoi(p + strlen(*l));
			}
		}

		p = t + 1;
	}

	free(s);
	return(mem * 1024);
}

#elif __APPLE__
/* dump vm_stat, sum free, inactive, speculative, and purgable */

/**
 * @fn mem_dump_vm_stat
 */
static inline
char *mem_dump_vm_stat(
	void)
{
	FILE *fp = NULL;
	char *res = NULL;

	/* open */
	if((fp = popen("vm_stat", "r")) == NULL) {
		goto _mem_dump_vm_stat_error_handler;
	}

	/* dump */
	if((res = mem_dump_file(fp)) == NULL) {
		goto _mem_dump_vm_stat_error_handler;
	}

	/* close file */
	if(pclose(fp) != 0) {
		goto _mem_dump_vm_stat_error_handler;
	}
	return(res);

_mem_dump_vm_stat_error_handler:
	if(fp != NULL) { pclose(fp); }
	if(res != NULL) { free(res); }
	return(NULL);
}

/**
 * @fn mem_estimate_free_size
 */
static inline
uint64_t mem_estimate_free_size(
	void)
{
	char *s = mem_dump_vm_stat();
	if(s == NULL) {
		return(0);
	}

	char const *labels[] = {
		"Pages free:",
		// "Pages active:",
		"Pages inactive:",
		"Pages speculative:",
		// "Pages throttled:",
		// "Pages wired down:",
		"Pages purgeable:",
		NULL
	};

	/* supposing s is terminated with double '\0's */
	uint64_t mem = 0;
	char *p = s;
	while(*p != '\0') {

		char *t = p;
		while(*t != '\n' && *t != '\0') { t++; }

		for(char const **l = labels; *l != NULL; l++) {
			if(strncmp(p, *l, strlen(*l)) == 0) {
				*t = '\0';
				mem += atoi(p + strlen(*l));
			}
		}

		p = t + 1;
	}

	free(s);
	return(mem * sysconf(_SC_PAGESIZE));
}
#endif

/**
 * end of mem.h
 */
