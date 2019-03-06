// vim: noai:ts=2:sw=2
/*
 * Lyndon Array Construction during Burrows-Wheeler Inversion
 *
 * Authors: Felipe A. Louza, W. F. Smyth, Guilherme P. Telles
 * contact: louza@ic.unicamp.br
 * 15/04/2016
 *
 */

#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <limits.h>

#include "lib/utils.h"
#include "lib/file.h"
#include "lib/bwt.h"
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

int VALIDATE=0, MODE=0, COUNT=0;

	if(argc!=7){
		dies(__func__,"argc!=7");
	}

	unsigned char **R;
	int_t i, n=0;
	int   k;

	char* c_dir = argv[1];
	char* c_file = argv[2];

	sscanf(argv[3], "%d", &k);
	sscanf(argv[4], "%u", &MODE);
	sscanf(argv[5], "%u", &VALIDATE);
	sscanf(argv[6], "%u", &COUNT); //count Lyndon factors 

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

  char* copy = NULL;
  if(VALIDATE==1 && MODE==5){
		copy = (char*) malloc((n+1)*sizeof(char));
		strcpy(copy, (char*)str);
  }

	//sorted array
	uint_t *LA = (uint_t*) malloc(n*sizeof(int_t));
	for(i=0; i<n; i++) LA[i]=0;

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
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;
		
		case 4:	printf("## MAX_LYN ##\n"); 
			compute_lyndon_max_lyn(str, (uint_t*)LA, n);
			break;

		case 5:	printf("## BWT_INPLACE_LYN ##\n"); 
			bwt_lyndon_inplace((char*)str, (uint_t*)LA, n);
			break;

		default: break;
	}

	// validate	
	if(VALIDATE==1){
		
		if(MODE==5){
			free(str);
			str = (unsigned char*) copy;
		}
		
		if(!lyndon_check(str, LA, n, 0)){
			printf("isNOTLyndonArray!!\n");
			fprintf(stderr, "ERROR\n");
		}
		else {
			printf("isLyndonArray!!\n");
		}
	}

	#if DEBUG
	for(i=0; i<min(n,20); i++){

		printf("%" PRIdN ") %" PRIdN "\t", i, LA[i]);
  	printf("%c\t", str[i]-1);
		int_t j=i;
		for(j=i; j<(int_t) min(i+10,i+LA[i]+1); j++)
			printf("%c", str[j]-1);
		printf("\n");
	}
	#endif

	if(COUNT==1){
		i=0;
		int_t count=0, max=0;
		while(i<n){
			#if DEBUG
				printf("%d\t%d\n", i, LA[i]);
			#endif
			count++;
			if(LA[i]>max) max=LA[i];
			i+=LA[i];
		}

		printf("##\n");
		printf("Number of Lyndon factors: % "PRIdN"\n", count);
		printf("Average length: %.2lf\n", (double)n/(double)count);
		printf("Maximum length: % "PRIdN"\n", max);
		printf("##\n");
	}

	free(LA);
	free(str);

return 0;
}

