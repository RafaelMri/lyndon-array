CC = gcc
CFLAGS += -Wall 
#CFLAGS += -g -O0
CFLAGS += -D_FILE_OFFSET_BITS=64 -m64 -O3 -fomit-frame-pointer -Wno-char-subscripts 

LFLAGS = -lm -lrt -ldl

DIR = dataset/
INPUT = input.100.txt

K	= 1
MODE 	= 1
CHECK	= 1
OUTPUT	= 0
COUNT = 1

LIBOBJ = \
	lib/utils.o\
	lib/file.o\
	lib/bwt.o\
	gsaca_cl/gsaca.o\
	../external/malloc_count/malloc_count.o\
	../external/sacak-lcp.o\
	../lyndon-array.o
	
##
M64 = 0
DEBUG = 0
PERMUTED = 0
STEP_TIME = 1
CHECK = 1
SAVE_SPACE = 1
COMPUTE_BWT = 1
##

DEFINES = -DDEBUG=$(DEBUG) -DM64=$(M64) -DPERMUTED=$(PERMUTED) -DSTEP_TIME=$(STEP_TIME)  -DSAVE_SPACE=$(SAVE_SPACE) -DCOMPUTE_BWT=$(COMPUTE_BWT)

CFLAGS += $(DEFINES)

all: compile

clean:
	\rm -f *.o ../*.o ../external/*.o lib/*o main 

compile: main.c ${LIBOBJ} 
	$(CC) $(CFLAGS) $(LFLAGS) -o main main.c ${LIBOBJ}

run:
	./main $(DIR) $(INPUT) $(K) $(MODE) $(CHECK) $(COUNT)

valgrind:
	valgrind --tool=memcheck --leak-check=full --track-origins=yes ./main $(DIR) $(INPUT) $(K) $(MODE) $(CHECK) $(COUNT)
