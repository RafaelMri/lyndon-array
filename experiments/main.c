/*
 * Optimal suffix sorting and LCP array construction for constant alphabets
 *
 * Authors: Felipe A. Louza, Simon Gog, Guilherme P. Telles
 * contact: louza@ic.unicamp.br
 * 15/01/2016
 *
 */

#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <limits.h>

#include "lib/utils.h"
#include "lib/file.h"
#include "../external/malloc_count/malloc_count.h"
#include "../external/sacak-lcp.h"
#include "gsaca_cl/gsaca.h"
#include "../lyndon-array.h"

#ifndef DEBUG
	#define DEBUG 0 
#endif

/*******************************************************************/
unsigned char* cat_char(unsigned char** R, int k, int_t *n){

	(*n)++; //add 0 at the end

	int_t i, j;
	int_t l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		for(j=0; j<m; j++){
			if(R[i][j]<255) str[l++] = R[i][j]+1;
//			str[l++] = R[i][j];
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
	if(*n>l){
	  str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char)); 
	  //printf("str_realloc(%" PRIdN ", %" PRIdN ")\n", *n, l);
	}

	*n = l;

return str;
}

/*******************************************************************/

int main(int argc, char** argv){

time_t t_start=0;
clock_t c_start=0;

int VALIDATE=0, MODE=0;

	if(argc!=6){
		dies(__func__,"argc!=4");
	}

	unsigned char **R;
	int_t i, n=0;
	int   k;

	char* c_dir = argv[1];
	char* c_file = argv[2];

	sscanf(argv[3], "%d", &k);
	sscanf(argv[4], "%u", &MODE);
	sscanf(argv[5], "%u", &VALIDATE);

	file_chdir(c_dir);

	//disk access
	R = (unsigned char**) file_load_multiple(c_file, k, &n);
	if(!R){
		fprintf(stderr, "Error: less than %d strings in %s\n", k, c_file);
		return 0;
	}

	//concatenate strings
	unsigned char *str = NULL;
	str = cat_char(R, k, &n);

	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));

	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	//free memory
	for(i=0; i<k; i++)
		free(R[i]);
	free(R);

	//sorted array
	uint_t *LA = (uint_t*) malloc(n*sizeof(int_t));

	switch(MODE){

		case 1: printf("## LYNDON_BWT ##\n"); 
			compute_lyndon_bwt(str, (uint_t*)LA, n);
			break;

		case 2:	printf("## LYNDON_NSV ##\n"); 
			compute_lyndon_nsv(str, (uint_t*)LA, n);
			break;

		case 3:	printf("## GSACA-lyndon ##\n"); 
			time_start(&t_start, &c_start);
			gsaca_cl(str, (uint_t*)LA, n);
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;
		
		default: break;
	}

	// validate	
	if(VALIDATE==1){
		if(!lyndon_check(str, LA, n, 0)) printf("isNOTLyndonArray!!\n");
		else printf("isLyndonArray!!\n");
	}


	#if DEBUG
//		suffix_array_print(SA, (unsigned char*)str, min(10,n), sizeof(char));
	#endif

	free(LA);
	free(str);

return 0;
}

