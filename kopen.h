
/**
 * @file kopen.h
 *
 * @brief a header for kopen and kclose
 */
#ifndef _KOPEN_H_INCLUDED
#define _KOPEN_H_INCLUDED

/**
 * @fn kopen
 */
void *kopen(const char *fn, int *_fd);

/**
 * @fn kclose
 */
int kclose(void *a);


#endif /** #ifndef _KOPEN_H_INCLUDED */
/**
 * end of kopen.h
 */
