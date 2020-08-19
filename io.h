#ifndef __IO_H__
#define __IO_H__

#include <stddef.h>	/* FILE, fprintf() */
#include <curses.h>	/* WINDOW, wprintw() */

/**
 * Seamlessly handle presence/absence of ncurses window.
 */
#define PRINT(wp, ...) if ((wp)) wprintw((wp), __VA_ARGS__); else fprintf(stderr, __VA_ARGS__);

/* print vectors */
void fprint_doubles(FILE *fp, double *v, size_t len, int precision, int newline);
void wprint_doubles(WINDOW *wp, double *v, size_t len, int precision, int newline);
void fprint_uints(FILE *fp, unsigned int *v, size_t len, int width, int newline);
void wprint_uints(WINDOW *fp, unsigned int *v, size_t len, int width, int newline);
void fprint_size_ts(FILE *fp, size_t *v, size_t len, int width, int newline);
void wprint_size_ts(WINDOW *fp, size_t *v, size_t len, int width, int newline);
void fprint_uchars(FILE *fp, unsigned char *v, size_t len, int width, int newline);
void wprint_uchars(WINDOW *fp, unsigned char *v, size_t len, int width, int newline);

/* read vectors */

/* print matrices */
void fprint_vectorized_sq_matrix(FILE *fp, double *mat, size_t n, int row);
void wprint_vectorized_sq_matrix(WINDOW *fp, double *mat, size_t n, int row);
void fprint_vectorized_matrix(FILE *fp, double *mat, size_t n, size_t l, int row);
void wprint_vectorized_matrix(WINDOW *wp, double *mat, size_t n, size_t l, int row);
void fprint_vectorized_uintmatrix(FILE *fp, unsigned int *mat, unsigned int n, unsigned int l,int row);

#endif
