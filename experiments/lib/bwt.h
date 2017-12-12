#ifndef BWTH
#define BWTH

#include <string.h>
#include <stdio.h>
#include <limits.h>

#include "utils.h"

/*******************************************************************/

int rank(char *T, char c, int i);

int bwt_lcp_inplace(char *T, int n, int *LCP);

int bwt_inplace(char *T, int n);

char* bwt_reverse(char *T, int n);

//computes the permuted Lyndon array (suffix order) in quadratic time
int bwt_lyndon_inplace(char *T, uint_t *LA, int n);

/*******************************************************************/

#endif
