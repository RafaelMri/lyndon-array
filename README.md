# lyndon-array

This algorithm computes the Lyndon-array (LA) of string s[0, n-1] during the BWT inversion.

## Introduction

First, we compute SA \[1\] and BWT of s[0, n-1], then during the BWT-inversion we compute LA.

## Build requirements

An ANSI C Compiler (e.g. GNU GCC)


## API

```c
/** @brief computes the lyndon array of string s[0, n-1] during the BWT inversion 
 *
 *  @param s	input string with s[n-1]=0
 *  @param A 	lyndon array 
 *  @param n	string length
 *  @return -1 if an error occured, otherwise the depth of the recursive calls.
 */
int compute_lyndon_bwt(unsigned char *s, uint_t *A, uint_t n);
```

## Example

**Compilation:**

```sh
gcc -c external/sacak-lcp.c
gcc -c external/malloc_count/malloc_count.c
gcc -c lyndon-array.c
gcc test.c -o test sacak-lcp.o lyndon-array.o malloc_count.o -ldl
```

**Run a test:**

```c
./test graindraining
```

**Output:**

```c
sizeof(int_t) = 4 bytes
Text = graindraining$
TOTAL:
CLOCK = 0.000033 TIME = 0.000000
i	SA	pLA	BWT	suffixes
0	13	1	g	$
1	2	11	r	aindraining$
2	7	6	r	aining$
3	5	2	n	draining$
4	12	1	n	g$
5	0	2	$	graindraining$
6	3	2	a	indraining$
7	10	2	n	ing$
8	8	2	a	ining$
9	4	1	i	ndraining$
10	11	1	i	ng$
11	9	1	i	ning$
12	1	1	g	raindraining$
13	6	1	d	raining$
isLyndonArray!!
malloc_count ### exiting, total: 14,626, peak: 11,386, current: 0
```

## References

\[1\] Nong, G., Practical linear-time O(1)-workspace suffix sorting for constant alphabets, ACM Trans. Inform. Syst., vol. 31, no. 3, pp. 1-15, 2013