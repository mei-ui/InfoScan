/*
 * Copyright (C) 2006 Sun Yat-Sen University
 * JianHua Yang(yjhua2110@yahoo.com.cn)
 *
 * This file is part of snoSeeker.
 *
 * snoSeeker is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * snoSeeker is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with snoSeeker; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef UTILS_H
#define UTILS_H

typedef int boolean;
#define FALSE 0
#define TRUE  1

#define EXIT_FAILURE 1
#define EXIT_SUCCESS 0

#define LINE_LEN 512

#ifndef MIN
#define MIN(x, y) (x)<(y)?(x):(y)
#endif

#ifndef MAX
#define MAX(x, y) (x)>(y)?(x):(y)
#endif

/* Some stuff to support large files in Linux. */
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE 1
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

extern char **splitWhitespace(char *s, int *entryNum);
/*split the line with whitesapce */

extern void freeWords(char **s, int entryNum);
/* free words */

extern void *safeMalloc(size_t bytes);
/* safety for malloc */

extern void *safeRealloc(void *block, size_t bytes);
/* safety for realloc */

extern void safeFree(void *block);
/* safety for free() */

extern char **splitString(char *s, char *delim, int *entryNum);
/*split the line with whitesapce */

extern boolean isDelim(char *delim, char c);
/* */

extern char *strClone(char *s);
/* clone string */

extern char *getLine(FILE *fp);
/* get whole line
 * modified from Vienna RNAfold utils.c 2008/10/4 11:08:01
*/

extern char *skipStartWhitespace(char *s);
/* skip leading whitespace */

extern boolean startStr(char *s, char *t);
/* s start with t */

extern int overlapLength(int qStart, int qEnd, int tStart, int tEnd);
/* caculate the overlap length */

/* follows function from Jim Kent jksrc/common.h  Thanks Jim */
extern void reverseBytes(char *bytes, long length);
/* reverse bytes */

extern void reverseComp(char *bytes);
/* reverse complement */

extern void complement(char *bytes);
/* complement sequence */

extern void toUpperStr(char *string);

#endif /* End UTILS_H */
