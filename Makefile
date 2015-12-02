CC=g++
CFLAGS=-O3 -std=c++11 -I -Wno-deprecated -fopenmp $(DEBUG)
OBJS=variantFinder.o ReadCluster.o spooky.o kmerAnalysis.o utilities.o
DEBUG=-g

#DEPS=ReadCluster.h spooky.h kmerAnalysis.h utilities.h

variantFinder : $(OBJS)
	$(CC) -o variantFinder $(OBJS) $(CFLAGS) 

variantFinder.o: variantFinder.cpp ReadCluster.h kmerAnalysis.h utilities.h 
	$(CC) $(CFLAGS) -c variantFinder.cpp

ReadCluster.o: ReadCluster.cpp ReadCluster.h spooky.h
	$(CC) $(CFLAGS) -c ReadCluster.cpp

spooky.o: spooky.h spooky.cpp
	$(CC) $(CFLAGS) -c spooky.cpp

kmerAnalysis.o: kmerAnalysis.cpp kmerAnalysis.h ReadCluster.h spooky.h
	$(CC) $(CFLAGS) -c kmerAnalysis.cpp

utilities.o: utilities.cpp utilities.h
	$(CC) $(CFLAGS) -c utilities.cpp


#ReadCluster.o:ReadCluster.cpp ReadCluster.h
#	$(CC) -c ReadCluster.cpp $(CFLAGS)

#variantFinder: variantFinder.o readCluster.o spooky.o kmerAnalysis.o utilities.o
#	$(CC) -o variantFinder variantFinder.o readCluster.o spooky.o kmerAnalysis.o utilities.o $(CFLAGS)