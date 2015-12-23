#ifndef READCLUSTER_H
#define READCLUSTER_H

#include<sparsehash/sparse_hash_map>
#include<sparsehash/dense_hash_map>
//#include<bitset>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<boost/dynamic_bitset.hpp>
//#include<algorithm>
//#include <boost/algorithm/string.hpp>
//#include<math.h>
#include<stdint.h>
//#include <stddef.h>
//#include<stdio.h>
//#include <iterator>
#include<sstream>
//#include<memory>
#include<ctime>
#include <omp.h>
#include <iostream>
#include<list>

#include "utilities.h"



class ReadCluster{
  
private:
  std::vector<std::string> readSeqs; //array of reads that all have a set of unique kmers in common and hence are clustered together
  std::string contig; //assembled contig from the collection of reads in readSeq
  std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> > kmerPositions; //a vector of hash tables each hash table contains the position of each kmer in a read the 
  //key is the kmer compressed as a 64 bit integer the value is the position on the read of
  //the start of the kmer 
  std::vector<std::pair<uint_fast64_t, long>> kmerCounts; //is a vector of pairs the first value  in the pair is the kmer the second value is the count of how many reads contain that kmer 
  int clusterKmerSize; //variable to hold the size of the kmers used to assemble the reads 
  std::vector<int> startPositions; //startPosition the kmer used to align the read with at least one other read in the cluster 
  //std::vector<bool> usedReads; //vector of boolean values each element of the vector represents one of the reads in the cluster. If a read has already been used in assembling the cluster its value is set to true

public:
  ReadCluster()=default;
 ReadCluster(long numReads):startPositions(numReads, -1){contig="0";} //parameter numReads is the number of reads that will go into this cluster
  void addSeq(std::string &); //functon to add a single read to the vector of reads
  std::vector<std::string> getSeqs();
  std::string getContig();
  void printReads();
  std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> >& getPositions(){ return kmerPositions;}; 
  std::string mergeReads(int kmerSize, int cutoffMinNuc);
  void addSequences(std::vector<std::string> &); //function to add a vector of reads
  void getKmers(); //function that will calculate the set of kmers for each read and store them in the class variables kmerPositions and kmerCounts 
  void setKmerSize(int size); 
  int getNumReads();
  void printKmerPositions();
  int getKmerSize(){return clusterKmerSize;};
  std::pair<uint_fast64_t, long> getPair(long index); //given an index return the pair kmer count pair corresponding to that index
  void setStartPositions(uint_fast64_t kmer); //given a kmer set the start position of that kmer for each read that contains the given kmer usually only called once for the kmer that occurs most often among the set of reads
  void printStartPositions();
  //void setUsedReadsFalse(); //will set the vector of usedReads all to false using the member variable clusterKmerSize if clusterKmerSize is not defined will throw an error

};

//very simple compare function to pass to sort in order to sort pairs by the second element
bool comparePairs(const std::pair<uint_fast64_t, long>&i, const std::pair<uint_fast64_t, long>&j);

#endif
