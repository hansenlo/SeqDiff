#ifndef UTILITIES_H
#define UTILITIES_H

#include<string>
#include<algorithm>

#include "spooky.h"

#include<iostream>

#include "sparsepp.h"

#include<sparsehash/sparse_hash_map>
#include<sparsehash/dense_hash_map>

#include <cmath>

//#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include<bitset>
#include <functional>

#include <unordered_map>

//#include<boost/dynamic_bitset.hpp>
//#include <boost/functional/hash.hpp>

#include<fstream>


const int bitSetSize=192; //the size of bitset objects for the entire program



struct stdHash{

  
  size_t operator()(const std::bitset<bitSetSize> &b) const{

    std::hash<std::bitset<bitSetSize>> hashObject;

    //return(1);

    return size_t(hashObject(b));
    

  }

};



struct customHash{

  
  size_t operator()(const unsigned long value ) const{


    return size_t(SpookyHash::Hash64(&value, 8, 20));
	       //return size_t(value);
  }

};


struct customHashVer2{

  
  //size_t operator()( boost::dynamic_bitset<> value) const{

  size_t operator()( const unsigned long value) const{


    size_t x;

    x=SpookyHash::Hash64(&value, sizeof(value), 20);

    std::cout<<"x is "<<x<<std::endl;
	       //return size_t(value);
  }

};




void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr);


std::string bit2String(std::bitset<bitSetSize> &number, int kmerSize);


void revComplement(std::string &seq);


void printBitString(std::vector <uint_fast64_t> &controlCtr1, std::vector <uint_fast64_t> &controlCtr2, std::string outFile);


//function to create a bit rev lookup table for use 
//in reverse complementing a nucleotide string
//arguments are the word size you want to reverse i.e. 8 bits 16 bits etc 
//and the lookup table which will be populated
//*******IMPORTANT******* bit reversal is for nucleotide strings!!! don't use bit reversal as universal reversal table*********
void createBitRevTable(int word, std::vector< std::bitset<bitSetSize> > &bitTable);

//function to create a bit rev lookup table for use 
//in reverse complementing a nucleotide string
//arguments are the word size you want to reverse i.e. 8 bits 16 bits etc 
//and the lookup table which will be populated
//*******IMPORTANT******* bit reversal is for nucleotide strings!!! don't use bit reversal as universal reversal table*********
//function to create a bit reversal table in the context of a machine word of 64 bits only 
void createBitRevTableMachineWord(int word, std::vector<uint_fast32_t> &bitTable);

void revComplementBitString(std::bitset<bitSetSize> &reversedKey, std::bitset<bitSetSize> &key, std::bitset<bitSetSize> &clearbitWord, std::vector< std::bitset<bitSetSize> > &bitTable, int bitWordSize, int kmerSize);


void revComplementMachineWord(uint_fast64_t &key, uint_fast64_t &revKey, int clusterKmerSize, std::vector<uint_fast32_t> &bitTable);

void readInFasta(std::unordered_map<std::string, std::string> &genome, std::string fileName); //function to read in a multi fasta file and store it in a hash table

//given a hash table of unique words 
//a file containing assembled contigs 
//return a hash table where the key is the contig id and value is the four words whose counts need to be looked up the first two words are the two unique words coming from either side of the contig
//the next two words are the corresponding words that match the reference i.e. one bp further out from the unique words
void getUniqueWordsContigs(google::dense_hash_map<std::bitset<bitSetSize>, int, stdHash>  &uniqueKmers, std::string contigFile, std::unordered_map<uint_fast64_t, std::vector< std::bitset<bitSetSize> >, stdHash> &pairedUniqueKmers);

//function will search the string source find all occurances of the parameter find and replace them with whatever is in the variable replace 
void find_and_replace(std::string& source, std::string const& find, std::string const& replace);


//function taks as input an integer and will return the kmer coded for by the integer
std::string int2String(uint_fast64_t number, int kmerSize);

#endif
