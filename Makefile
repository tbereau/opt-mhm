CC=gcc
#CC_MP=icc -openmp -lm -wd981 -DOPENMP=1
CC_MP=gcc -fopenmp -lm -DOPENMP=1
CFLAGS=-c -Wall -O3
LINKER=-lm

# For a serial version of the code, type 'make serial'.
# For a threaded version, type 'make openmp'.
# Right now, the serial version uses g++ and the threaded version
# icpc, but one can try other combinations... 

all: serial openmp

openmp: opt-mhm_mp.o
	$(CC_MP) opt-mhm_mp.o -o opt-mhm_mp

serial: opt-mhm.o
	$(CC) $(LINKER) opt-mhm.o -o opt-mhm

opt-mhm.o: opt-mhm.c
	$(CC) $(CFLAGS) opt-mhm.c -o opt-mhm.o

opt-mhm_mp.o: opt-mhm.c
	$(CC_MP) $(CFLAGS) opt-mhm.c -o opt-mhm_mp.o

clean:
	rm -rf  *o

