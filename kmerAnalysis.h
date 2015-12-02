#ifndef KMERANALYSIS_H
#define KMERANALYSIS_H

#include "ReadCluster.h"


//node of a linked list will hold a set of sequences each thread will operate on
struct node {
  std::vector<std::string> chunk; //vector of strings contains the read sequence
  std::vector<std::string> qualityScores; //vector of strings contains the quality scores for each read use this for debugging purposes
};


inline void countKmers(google::sparse_hash_map<uint_fast64_t, int, customHash> &controlKmers, char nextLineFlag, std::string inputFile, int kmerSize);

//void countUniqueKmers(std::vector <uint_fast64_t> &controlCtr1, std::vector <uint_fast64_t> &controlCtr2,   sparse_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, std::string inputFile, int kmerSize)
//void countUniqueKmers(std::vector <uint_fast64_t> &controlCtr1, std::vector <uint_fast64_t> &controlCtr2,   dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers,  char nextLineFlag, std::string inputFile, int kmerSize)
void countUniqueKmers(google::sparse_hash_map<uint_fast64_t, int, customHash>  &controlKmers, google::sparse_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, std::string inputFile, int kmerSize);

void readUniqueKmers(google::dense_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, std::string inputFile, int kmerSize, int ctrCutoff);

//debugging function for a given cluster ID will print all kmers that point to that cluster
void printSingleCluster( google::dense_hash_map<uint_fast64_t, long, customHash> &clusterKmers, uint_fast64_t clusterID, int kmerSize);

	     	      
//function to print out the clusters of reads and all unique reads
void printClusters(std::vector< std::vector<std::string> > &clusterBuffer, google::dense_hash_map<uint_fast64_t, int, customHash> &clusterFiles, std::vector< std::string > &uniqueReadsBuffer, std::ofstream &uniqueOut, std::vector<std::shared_ptr<std::ofstream> > &files, int tid);

void assignClusters(node * workNodePtr, //a pointer pointing to a chunk of work 
		    google::dense_hash_map<uint_fast64_t, long, customHash> &clusterKmers, //hash table key is kmer value is the cluster that kmer belongs to
		    google::dense_hash_map<uint_fast64_t, int, customHash> &clusterFiles,  //hash table key is a cluster id value is the index into the vector of output file names
		    std::vector<std::shared_ptr<std::ofstream> > &files, //set of opened file handlers that clusters will be written to
		    google::dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers, //hash table key is a kmer that is unqiue to the exp read library value is the number of times kmer occurs in experiment read library
		    std::ofstream &uniqueOut, //ofstream object file handler contains location of output file for all unique reads
		    int kmerSize, //kmer size
		    int numFiles, //number of files to print clusters to
		    uint_fast64_t &clusterCtr, //count of how many clusters there currently are will be shared among threads
		    int tid //thread id number passing this for debugging purposes
		    );


//Given a hash table containing unique kmers and their counts get reads that have unique kmers above a certain threshold

std::vector<std::string> getReads(google::dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers, int numFiles, char nextLineFlag, std::string inputFile, int kmerSize);

//function to read in Clusters do some filtering and pass the filtered Clusters to be assembled
void readInClusters(std::vector<std::string> &fileNames);

#endif
