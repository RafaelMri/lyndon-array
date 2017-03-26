#include "lyndon-array.h"
#include "external/sacak-lcp.h"

#ifndef TIME
	#define TIME 1
#endif

#ifndef STEP_TIME
	#define STEP_TIME 1
#endif

#define SIGMA 255

//0 = Lyndon-array (LA)
//1 = permuted LA (in suffix-ordering)
#ifndef PERMUTED
	#define PERMUTED 0
#endif

#ifndef PRINT
	#define PRINT 0
#endif

#define STACK_SIZE 10240/(sizeof(uint_t)*2) //10KB

/**********************************************************************/

typedef struct _pair{
        uint_t i;
        uint_t j;
} s_pair;

typedef struct{
        s_pair *array;
        uint_t top;
        uint_t size;
} stack;

int stack_init(stack *S, uint_t size){

	S->array = (s_pair*) malloc(size*sizeof(s_pair));	
	if(S->array==NULL) return 0;

	S->size = size;
	S->top = 0;

return 0;
}

void stack_push(stack *S, uint_t i, uint_t j){

	if(S->top==S->size){
		S->size += STACK_SIZE;
		S->array = (s_pair*) realloc(S->array, S->size*sizeof(s_pair));
	} 

	S->array[S->top].i = i;
	S->array[S->top].j = j;

	(S->top)++;
}

s_pair stack_top(stack *S){

	s_pair aux;
	aux.i = S->array[S->top-1].i;
	aux.j = S->array[S->top-1].j; 

return aux;
}

void stack_pop(stack *S){
	(S->top)--;
}

void stack_print(stack *S){

	uint_t i;

	for(i=S->top-1; i>0; i--)
		printf("%d\t(%d, %d)\n", i, S->array[i].i, S->array[i].j);
	
	printf("%d\t(%d, %d)\n", i, S->array[i].i, S->array[i].j);

//	while(S->top>0){
//		printf("%d\t(%d, %d)\n", i, stack_top(S).i, stack_top(S).j);
//		stack_pop(S);
//	}

}

/****************************************************************************/

void tstart(time_t *t_time, clock_t *c_clock){
        *t_time = time(NULL);
        *c_clock =  clock();
}

double tstop(time_t t_time, clock_t c_clock){
        double aux1 = (clock() - c_clock) / (double)(CLOCKS_PER_SEC);
        double aux2 = difftime (time(NULL),t_time);
        printf("CLOCK = %lf TIME = %lf\n", aux1, aux2);
        return aux1;
}

/**********************************************************************/

int compute_lyndon_bwt(unsigned char *s, uint_t *A, uint_t n){

int_t i;

	#if TIME
		time_t t_total= 0;clock_t c_total= 0;
		tstart(&t_total, &c_total); 
	#endif

	#if STEP_TIME
		time_t t_time = 0;clock_t c_time = 0;
		tstart(&t_time, &c_time); 
		printf("1. SA construction:\n");
	#endif
	
	// 1. sort
	#if PERMUTED  	//1. uses 1 array for SA -> LF -> pLA
		uint_t *SA = A;
	#else		//0. uses 2 arrays, SA -> LF, and LA
		uint_t *SA = (uint_t*) malloc(n*sizeof(uint_t));
	#endif

	sacak(s, SA, n);

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("2. Obtain BWT:\n"); 
	#endif

	//2. compute BWT
	unsigned char* bwt = malloc(n*sizeof(unsigned char));
	for(i=0; i<n; i++){
		bwt[i] = (SA[i])?s[SA[i]-1]:'$';
	}	

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("3. Compute LF:\n"); 
	#endif

	//3. compute LF-array (in the space of SA[0,n-1])
	uint_t *LF = SA;
	uint_t *C = (uint_t*) malloc(SIGMA*sizeof(uint_t));

	//counter
	for(i=0; i<SIGMA; i++) C[i]=0;
	for(i=0; i<n; i++) C[bwt[i]]++;

	uint_t sum=0;
	for(i=0; i<SIGMA; i++){
		sum+=C[i]; C[i]=sum-C[i];
	}

	for(i=0; i<n; i++){
		LF[i] = C[bwt[i]]++;
	}

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("4. Compute Lyndon:\n");
	#endif

	//BWT-inversion and Lyndon array construction
	stack S;
	S.top=0; S.size=STACK_SIZE;
	stack_init(&S,STACK_SIZE);
	stack_push(&S, 0, 0);//(0, 0)

	#if PRINT
		printf("step\tpos\tT^{rev}\tLyndon\n");
	#endif

	uint_t *LA = A;

	uint_t pos = 0;
	uint_t step = 1;//(n-1)-i

	for(i=n-1; i >= 0; i--){	

		while(stack_top(&S).i > pos) stack_pop(&S);

		uint_t next = LF[pos];

		#if PERMUTED
			LA[pos] = step-stack_top(&S).j;
		#else
			LA[i] = step-stack_top(&S).j;
		#endif

		#if PRINT
			printf("%d\t%d\t%c\t%d\n", step, pos, bwt[pos], step-stack_top(&S).j);
		#endif

		stack_push(&S, pos, step++);

		pos = next; //pos = LF(pos)
	}

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
	#endif

	#if TIME
		printf("TOTAL:\n");
		fprintf(stderr,"%.6lf\n", tstop(t_total, c_total)); 
	#endif


	free(S.array);
	free(C);
	free(bwt);

	#if PERMUTED == 0
		free(SA);
	#endif

return 0;
}

/**********************************************************************/

int lyndon_check(unsigned char *s, uint_t *LA, uint_t n, int print){

int check=1;
uint_t i;

	uint_t *SA = (uint_t*) malloc(n*sizeof(uint_t));
	sacak(s, SA, n);

	uint_t *ISA = (uint_t*) malloc(n*sizeof(uint_t));
	for(i=0; i<n; i++) ISA[SA[i]]=i;

	// print output
	if(print){
		#if PERMUTED
			printf("i\tSA\tpLA\tBWT\tsuffixes\n");
		#else
			printf("i\tSA\tLA\tBWT\tsuffixes\n");
		#endif
	
		for(i = 0; i < n; ++i) {
	
			#if PERMUTED
				char j = (SA[i])? s[SA[i]-1]:'$';
			#else
				char j = (i)? s[i-1]:'$';
			#endif
			printf("%d\t%d\t%d\t%c\t",i, SA[i], LA[i],j);
			
			#if PERMUTED
			    int start=SA[i];
			#else
			    int start=i;
			#endif
	
			uint_t k;	
			for(k = start; k < n; ++k) {
				printf("%c", s[k]);
			}
			printf("$\n");
		}
	}

	//simple checker
	for(i=0; i<n; i++){
		
		#if PERMUTED
			int start=SA[i];
			int pos = i;
		#else
			int start=i;
			int pos = ISA[i];
		#endif

		if(start==n-1){
			if(LA[i]!=1){
				check=0;
				printf("*%d) %d != %d\n", i, LA[i], 1);
				break;
			}
			continue;			
		}

		int size=1;
		while(ISA[start+size]>pos)size++; //ISA[n-1]=0

		if(LA[i]!=size){
			check=0;
			printf("%d) %d != %d\n", i, LA[i], size);
			break;
		}
	}

	free(SA);
	free(ISA);
	
	return check;
}

/**********************************************************************/

int compute_nsv(uint_t* NSV, uint_t *array, int_t n){

	stack S;
	S.top=0; S.size=STACK_SIZE;
	stack_init(&S,STACK_SIZE);
	stack_push(&S, 0, -1);//(idx, value)=(pos, step)

        for(int_t i=n-1; i>0; i--){

                NSV[i] = stack_top(&S).i;
		
                if(array[i-1]>array[i]){
			stack_push(&S, i, array[i]);
                }
                else{
			while(stack_top(&S).j >= array[i-1]) stack_pop(&S);
                }
        }

        NSV[0] = stack_top(&S).i;

return 0;
}

/**********************************************************************/

int compute_lyndon_nsv(unsigned char *s, uint_t *A, uint_t n){

int_t i;

	#if TIME
		time_t t_total= 0;clock_t c_total= 0;
		tstart(&t_total, &c_total); 
	#endif

	#if STEP_TIME
		time_t t_time = 0;clock_t c_time = 0;
		tstart(&t_time, &c_time); 
		printf("1. SA construction:\n");
	#endif
	
	// 1. sort
	//uint_t *SA = A;
	uint_t *SA = (uint_t*) malloc(n*sizeof(uint_t));
	sacak(s, SA, n);

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("2. Compute ISA:\n"); 
	#endif

	//2. compute ISA
	uint_t *ISA = (uint_t*) malloc(n*sizeof(uint_t));
	for(i=0; i<n; i++) ISA[SA[i]]=i;

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("3. Compute NSV-array:\n"); 
	#endif

	//3. compute NSV

	uint_t *NSV = (uint_t*) malloc(n*sizeof(uint_t));
	//uint_t *NSV = A;
	compute_nsv(NSV, ISA, n);
	NSV[n-1]=n;

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
		tstart(&t_time, &c_time);
		printf("4. Compute Lyndon:\n");
	#endif

	//for(i=0; i<n; i++) printf("%d\t%d\t%d\t%d\t%d\n", i, SA[i], ISA[i], NSV[i], 0);

	#if PRINT
		printf("step\tpos\tT^{rev}\tLyndon\n");
	#endif

	uint_t *LA = A;

	uint_t pos = 0;
	for(i=0; i < n; i++){	

		pos = ISA[i];
		LA[pos] = NSV[pos]-pos;
	}

	#if STEP_TIME
		fprintf(stderr,"%.6lf\n", tstop(t_time, c_time));
	#endif

	#if TIME
		printf("TOTAL:\n");
		fprintf(stderr,"%.6lf\n", tstop(t_total, c_total)); 
	#endif

	free(SA);
	free(ISA);

return 0;
}

/**********************************************************************/

