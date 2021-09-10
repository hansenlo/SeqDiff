#ifndef READCLUSTER_H
#define READCLUSTER_H

#include<sparsehash/sparse_hash_map>
#include<sparsehash/dense_hash_map>
#include<bitset>
//#include<fstream>
#include<string>
#include<stdlib.h>
//#include<boost/dynamic_bitset.hpp>
//#include<algorithm>
//#include <boost/algorithm/string.hpp>
#include<math.h>
#include<cmath>
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
#include <functional>
#include <tuple>
#include<sstream>
#include<numeric>

#include "utilities.h"

#include <tclap/CmdLine.h>




// C++ program to print DFS traversal from a given vertex in a  given graph
 
// Graph class represents a directed graph using adjacency list representation
class Graph
{
  //int V;    // No. of vertices
    //std::list<int> *adj;    // Pointer to an array containing adjacency lists
    std::vector<std::list<int>> graph;

    void DFSUtil(int v, std::vector<bool> &visited, std::vector<int> &cluster);  // A function used by DFS
public:
 Graph(int V):graph(V){}// Constructor build a graph with a known number of vertices
    void addEdge(int v, int w);   // function to add an edge to graph
    void DFS(int v, std::vector<int> &cluster, std::vector<bool> &visited);    // DFS traversal of the vertices reachable from v
    void addVertex() //adding a new vertex if necassary
    {      
      std::list<int> temp; //adding a new vertex
      graph.push_back(temp);
    }
};
 


class ReadCluster{
  
private:
  std::vector<std::string> readSeqs; //array of reads that all have a set of unique kmers in common and hence are clustered together
  std::vector<std::string> qualityStrings; //array of quality strings corresponding to the reads in the cluster
  std::string contig; //assembled contig from the collection of reads in readSeq
  std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> > kmerPositions; //a vector of hash tables each hash table contains the position of each kmer in a read the 
  //key is the kmer compressed as a 64 bit integer the value is the position on the read of
  //the start of the kmer 
  google::dense_hash_map<uint_fast64_t, long, customHash> kmerInReadCounts; //is a hash table the key is the kmer the valur is the count of how many reads contain that kmer
  
  std::vector<bool> sameStrand;  //vector to keep track of whether a read is on the right strand if it is the value will be set to true; 

  //std::vector<std::pair<uint_fast64_t, long>> kmerCounts; //is a vector of pairs the first value  in the pair is the kmer the second value is the count of how many reads contain that kmer 
  int clusterKmerSize; //variable to hold the size of the kmers used to assemble the reads 
  std::vector<int> startPositions; //startPosition the kmer used to align the read with at least one other read in the cluster 
  //std::vector<bool> usedReads; //vector of boolean values each element of the vector represents one of the reads in the cluster. If a read has already been used in assembling the cluster its value is set to true
  
  spp::sparse_hash_map<uint_fast64_t, bool, customHash> presentMultipleTimes; //hash table to keep track of which kmers are found at least twice in any read if a kmer is found at least twice in a read 

public:
  ReadCluster()=default;
 ReadCluster(long numReads):startPositions(numReads, -1){contig="0"; kmerInReadCounts.set_empty_key(-1); } //parameter numReads is the number of reads that will go into this cluster
  void addSeq(std::string &); //functon to add a single read to the vector of reads
  std::vector<std::string> getSeqs();
  std::string getContig();
  void printReads();
  std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> >& getPositions(){ return kmerPositions;}; 
  void mergeReads(std::vector<std::string> &contigs, int kmerSize, int cutoffMinNuc, std::ofstream &debugging, std::vector<std::string> &clusterID); //kmerSize is the size of the kmers cutoffMinNuc is the minimum number of nucleotides in a column to call a base debugging and cluster ID are purely for debugging purposes.  

  //void addSequences(google::dense_hash_map<std::string, int, customHash> *); //function to add a vector of reads

  void addSequences(std::unordered_map<std::string, std::string> *allReads); //function to add the set of reads that belong to a cluster

  //need to add this subroutine. Function will take as input the most common kmer present in the most reads and will ensure all reads are on the same strand as the most common kmer
  void revCompCluster(uint_fast64_t maxKmer);

  uint_fast64_t getKmers(); //function that will calculate the set of kmers for each read and store them in the class variable kmerPositions functions returns the kmer that occurs most often in the reads that belong to the cluster 
  void setKmerSize(int size); 
  int getNumReads();
  void printKmerPositions();
  int getKmerSize(){return clusterKmerSize;};
  //`std::pair<uint_fast64_t, long> getPair(long index); //given an index return the pair kmer count pair corresponding to that index
  void setStartPositions(uint_fast64_t kmer); //given a kmer set the start position of that kmer for each read that contains the given kmer usually only called once for the kmer that occurs most often among the set of reads
  void printStartPositions();
  //void setUsedReadsFalse(); //will set the vector of usedReads all to false using the member variable clusterKmerSize if clusterKmerSize is not defined will throw an error
  
  uint_fast64_t numberDiff(std::vector<std::vector<char>> &alignmentMatrix, int row1, int row2, int start1, int start2, int &sizeAligned, int readSize); //Given the alignment matrix of reads and the row index and start positition of the 2 reads calculate percent difference of the 2 reads for the aligned regions sizeAligned is the number of aligned bases between the 2 reads the aligned bases may not necassarily match

   //given a matrix of aligned reads merge the columns into a contig also provided is a vector of rows that should 
  //be used in the assembly rows indexes not listed in the vector will be ignored 
  //also returned is the positions of the Ns in the contig
  //startContig contains the starting column in the alignmentMatrix of where the contig is starting to be assembled
  //longestDistN contains the longest distance between Ns for the contig
  //subClusterCtr contains a count of how many columns have a high second highest nuclotide count 
  //problematicCols is a list of columns which have a large number of counts for the second highest nucleotide
  std::string assembleContig(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, double &nucCtr, double &badColCtr, int cutoffMinNuc, double &Nctr, double &startContig, double &longestDistN,
			     int &subClusterCtr, std::vector<int> &problematicCols);
  
  //std::vector<double> &Npositions

  void getDistanceGraph(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &startMatrix, Graph &numMatches, int readSize); //given an alignment matrix return the matrix of distances between pairs of reads startMatrix is the start position of each read in the alignment matrix

  void getSubClusters(std::vector<std::vector<int>> &numMatches, std::vector<std::vector<int>> &clusters);


  //debugging function only will print the alignment Matrix for the set of rows passed to it 
  void printMatrix(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, double &percentBadCol, std::string &clusterID, std::ofstream &debugging, double &percentNs); 
 
  //checks to see if columns that are Ns are approximately 50% to different bases if so will return contigs for both cases if successfull function will return a 1 and a vector of new contigs
  bool extractHetro(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &newContigs, int startContig, int cutoffMinNuc);
 

  //trys to split a contig into two different contigs representing two different variants 
  bool extractVariants(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &newContigs, double &startContig, int cutoffMinNuc, int readSize);
 
  bool assembleSubClusters(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &newContigs, double &startContig, int cutoffMinNuc, int readSize, std::vector<int> &problematicCols);
 



//function to check the quality of the assembled contig returns through the parameter 
  //list several different metrics of how well the assembly worked
  bool checkContig(std::string &combinedNuc, int readSize, double &badColCtr, double &Nctr, double &percentBadCol, double &percentNs) 
  {                                                                               

    percentBadCol=(badColCtr/combinedNuc.length());
    percentNs=(Nctr/combinedNuc.length());

    //return; //debugging must remove!!

 
    if(combinedNuc.length() < (readSize/2)) //if contig to short throw it out
      {
	return(false);
	//combinedNuc="0";
      }

    if(percentBadCol > 0.1) //if to many bad columns throw it out
      {

	return(false);
	//combinedNuc="0";
      }

  //if number of not confident nucleotide calls is greater than 10% of the contig throw it out
    if(percentNs > 0.05)
      {
	return(false);
	//combinedNuc="0";
      }
  

    return(true);
  }



};

//very simple compare function to pass to sort in order to sort pairs by the second element
bool comparePairs(const std::pair<uint_fast64_t, long>&i, const std::pair<uint_fast64_t, long>&j);


#endif


