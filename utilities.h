#ifndef UTILITIES_H
#define UTILITIES_H

#include<string>
#include<algorithm>

#include "spooky.h"

#include<iostream>


//#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include<bitset>
#include <functional>

#include<boost/dynamic_bitset.hpp>
//#include <boost/functional/hash.hpp>

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



void revComplementBitString(std::bitset<bitSetSize> &reversedKey, std::bitset<bitSetSize> &key, std::bitset<bitSetSize> &clearbitWord, std::vector< std::bitset<bitSetSize> > &bitTable, int bitWordSize, int kmerSize);



#endif
