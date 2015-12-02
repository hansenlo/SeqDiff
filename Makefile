CC=g++
CFLAGS=-g -O3 -std=c++11 -I -Wno-deprecated -fopenmp
DEPS=ReadCluster.h spooky.h kmerAnalysis.h utilities.h

#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)


ReadCluster.o:ReadCluster.cpp ReadCluster.h
	$(CC) -c ReadCluster.cpp $(CFLAGS)

#variantFinder: variantFinder.o readCluster.o spooky.o kmerAnalysis.o utilities.o
#	$(CC) -o variantFinder variantFinder.o readCluster.o spooky.o kmerAnalysis.o utilities.o $(CFLAGS)