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

void readUniqueKmers(google::dense_hash_map<std::bitset<bitSetSize>, int, stdHash>  &uniqueKmers,  char nextLineFlag, std::string inputFile, int kmerSize, int ctrCutoff);

//debugging function for a given cluster ID will print all kmers that point to that cluster
void printSingleCluster( google::dense_hash_map<uint_fast64_t, long, customHash> &clusterKmers, uint_fast64_t clusterID, int kmerSize);

	     	      
//function to print out the clusters of reads and all unique reads
void printClusters(std::vector< std::vector<std::string> > &clusterBuffer, google::dense_hash_map<uint_fast32_t, int, stdHash> &clusterFiles, std::vector<int> &qualityBuffer, std::vector< std::string > &qualityStringBuffer,  std::vector<std::shared_ptr<std::ofstream> > &files);

//function to print out unique reads only header of unique reads is (int read id)_(int clusterID)_(int quality flag)
void printUniqueReads(std::vector< std::string > &uniqueReadsBuffer, std::vector<int> &qualityBuffer, std::ofstream &uniqueOut, int tid);



void assignClusters(node * workNodePtr, //a pointer pointing to a chunk of work 
		    google::dense_hash_map<std::bitset<bitSetSize>, uint_fast32_t, stdHash> &clusterKmers, //hash table key is kmer value is the cluster that kmer belongs to
		    google::dense_hash_map<uint_fast32_t, int, stdHash> &clusterFiles,  //hash table key is a cluster id value is the index into the vector of output file names
		    google::dense_hash_map<std::bitset<bitSetSize>, int, stdHash>  &uniqueKmers, //hash table key is a kmer that is unqiue to the exp read library value is the number of times kmer occurs in experiment read library
		    std::ofstream &uniqueOut, //ofstream object file handler contains location of output file for all unique reads
		    int kmerSize, //kmer size
		    uint_fast64_t &clusterCtr, //count of how many clusters there currently are will be shared among threads
		    int tid, //thread id number passing this for debugging purposes
		    int numFiles, //number of files to print clusters to
		    std::unordered_map<double, double > &linkClusters //hash table key is a cluster id the value is the cluster id that the key cluster should be merged with  
		    );

void sendClustersToFile(std::string &uniqueReadFile, //the file location where all the reads containing novel words were placed 
			google::dense_hash_map<uint_fast32_t, int, stdHash> &clusterFiles,  //hash table key is a cluster id value is the index into the vector of output file names
			std::vector<std::shared_ptr<std::ofstream> > &files, //set of opened file handlers that clusters will be written to
			std::unordered_map<double, double > &linkClusters,  //hash table key is a cluster id the value is the cluster id that the key cluster should be merged with 
			int numFiles ); //number of files to print clusters to

//hash table key is a cluster id the value is the cluster id that the key cluster should be merged with 
//function will take a hash table of links to clusters and flag any clusters which have to many links to other clusters by putting a -1 for the value
void filterClusters(std::unordered_map<double, double > &linkClusters, //hash table key is a cluster id the value is the cluster id that the key cluster should be merged with 
		    int maxConnectivity); //max number of clusters that can be connected any more connections and the cluster is flagged as a bad cluster




//Given a hash table containing unique kmers and their counts get reads that have unique kmers above a certain threshold

//int numberCalled This variable holds the number of times this function was called per input sequencing file

std::vector<std::string> getReads(google::dense_hash_map<std::bitset<bitSetSize>, int, stdHash>  &uniqueKmers, int numFiles, char nextLineFlag, std::string inputFile, int kmerSize);

//function to read in Clusters do some filtering and pass the filtered Clusters to be assembled
void readInCluster(std::string &fileName, int cutoffClusterSize, int clusterKmerSize, int tid,  std::ofstream &contigOut, long &clusterNumber, std::ofstream &debuggingMatrix); //tid is for debugging purposes variable stores the thread id

void mergeClusters(std::vector<std::string> &fileNames);

#endif
