#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

/* taken from samtools' kstring.c */
typedef struct {
	uint64_t tab[4];
	int sep, finished;
	const char *p; // end of the current token
} ks_tokaux_t;
/* END - taken from samtools' kstring.c */

static inline void ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		s->m = size;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

static inline int kputw(int c, kstring_t *s)
{
	char buf[16];
	int l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	long l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

//WCH: unsigned long
//TODO: review code and test to make sure that we do not overflow!
static inline int kputul(unsigned long c, kstring_t *s)
{
	char buf[32];
	long l, i;
	unsigned long x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}
//END-WCH:

int ksprintf(kstring_t *s, const char *fmt, ...);
int ksplit_core(char *s, int delimiter, int *_max, int **_offsets);

/* taken from samtools' kstring.c */
static inline int *ksplit(kstring_t *s, int delimiter, int *n)
{
	int max = 0, *offsets = 0;
	*n = ksplit_core(s->s, delimiter, &max, &offsets);
	return offsets;
}
/* END - taken from samtools' kstring.c */

static inline int ks_set(const char *p, kstring_t *s)
{
	if (NULL==p) {
		s->l = 0;
		return 0;
	}
	
	size_t l = strlen(p);
	
	if (s->m < l) {
		s->m = l+1; // WCH: need to allocate for '\0'
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s, p, l);
	s->l = l;
	s->s[s->l] = 0;
	//TODO: unsigned and signed conversion!
	return (int)l;
}

#endif
