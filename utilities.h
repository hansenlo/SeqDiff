#ifndef UTILITIES_H
#define UTILITIES_H

#include<string>
#include<algorithm>


void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr);


std::string bit2String( uint_fast64_t number, int kmerSize);


void revComplement(std::string &seq);


void printBitString(std::vector <uint_fast64_t> &controlCtr1, std::vector <uint_fast64_t> &controlCtr2, std::string outFile);


//function to create a bit rev lookup table for use 
//in reverse complementing a nucleotide string
//arguments are the word size you want to reverse i.e. 8 bits 16 bits etc 
//and the lookup table which will be populated
//*******IMPORTANT******* bit reversal is for nucleotide strings!!! don't use bit reversal as universal reversal table*********
void createBitRevTable(int word, std::vector<uint_fast32_t> &bitTable);

#endif
