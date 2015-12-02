#include "utilities.h"

#include <iostream>
#include<fstream>
#include<boost/dynamic_bitset.hpp>

using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::ofstream;
using std::ios;
using std::cerr;
using boost::dynamic_bitset;


struct eqKey{

  bool operator()(const int a, const int b) const
  {
    if(a==b)
      {
	return true;
      }else
      {

	return false;
      }

  }

};


void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr)
{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
}


string bit2String( uint_fast64_t number, int kmerSize)
{
  
  int i, numberNuc;
  string bitString, nucString, nucBit;

  nucString="";
  
  dynamic_bitset<> bitRep(64, number);
  
  
  to_string(bitRep, bitString);



  for(i=64-(kmerSize*2); i<64; i+=2)
    {
      nucBit=bitString.substr(i, 2);
      

      if(nucBit=="00")
	{

	  nucString=nucString+"A";

	}else if(nucBit=="11")
	{
	  nucString=nucString+"T";

	}else if (nucBit=="10")
	{
	  nucString=nucString+"G";

	}else if (nucBit=="01")
	{
	  nucString=nucString+"C";
	}

      
    
    }

  //cout<<"nuc String is "<<nucString<<endl;
  //cout<<"bit string is "<<bitString<<endl;

  
  return(nucString);
  


  //return bitString;
}








void revComplement(string &seq){

  int i;

  //reverse the string
  std::reverse(seq.begin(), seq.end());

  //complementing the string
  for(i=0; i<seq.length(); i++)
    {
      
      //complementing the string
      switch(seq[i])
		    {
		    case 'A' :
		      seq[i]='T';
		      break;
		    
		    case 'a' :
		      seq[i]='T';
		      break;
		   
		    case 'T' :
		      seq[i]='A';
		      break;
		    
		    case 't' :
		      seq[i]='A';
		      break;
		    
		    case 'G' :
		      seq[i]='C';
		      break;
		    
		    case 'g' :
		      seq[i]='C';
		      break;

		    case 'C' :
		      seq[i]='G';
		      break;

		    case 'c' :
		      seq[i]='G';
		      break;

		    case 'N' :
		      seq[i]='N';
		      break;

		    case 'n' :
		      seq[i]='N';
		      break;

		      
		    default:
		      cerr<<"line is "<<seq<<endl;
		      cerr<<"Base is not a nucleotide exiting in revComplement function "<<endl;
		      exit(EXIT_FAILURE);

		    }




    }

    


}



void printBitString(vector <uint_fast64_t> &controlCtr1, vector <uint_fast64_t> &controlCtr2, string outFile)
{
  ofstream fout(outFile.c_str());
  uint_fast64_t numberKmers=0;
  uint_fast64_t i, size, kmer;
  int bit, j;

  size=controlCtr1.size();

  //cout<<"size is "<<size<<endl;

  cerr<<"Printing Kmers "<<endl;

  for(i=0; i<size; i++)
    {
      /*
      if(i>10)
	return;

      dynamic_bitset<> bitRep2(8, controlCtr1[i]);
      cout<<" index is  "<<i<<endl;
      cout<<bitRep2<<endl;
      */

      for(j=0; j<64; j++)
	{
	  if((controlCtr1[i] & (1 << j))>0)
	    {
	      kmer=(i*64)+j;
	      fout<<kmer<<"\t"<<1<<"\n";
	      //cout<<kmer<<"\t"<<1<<"\n";
	      

	    }
	}

      
    }
  



}


//function to create a bit rev lookup table for use 
//in reverse complementing a nucleotide string
//arguments are the word size you want to reverse i.e. 8 bits 16 bits etc 
//and the lookup table which will be populated
//*******IMPORTANT******* bit reversal is for nucleotide strings!!! don't use bit reversal as universal reversal table*********
void createBitRevTable(int word, std::vector<uint_fast32_t> &bitTable){

  unsigned long long int size=pow(2, word);
  unsigned long long int i, j, revBitStr;
  unsigned long long int temp=0;
  

  //create appropriate bit string to reset 
  unsigned long long int reset=pow(2,word)-1;

  //std::vector<uint_fast64_t> bitTable(size);

  //loop though all numbers that can be created with word size of bits
  for(i=0; i<size; i++)
    {
      revBitStr=0;

      /*
      if(i % 10000000==0)
	{
	  cerr<<"number of integers reversed  is "<<i<<endl;
	}
      */

      //dynamic_bitset<> bitRepKey(64, i);
      //cout<<bitRep<<endl;
      //cout<<"bit string key is "<<"  "<<i<<" "<<bitRepKey<<endl;

      for(j=0; j<word; j+=2)
	{
	  //shift two bits into first two bit positions
	  //clear all other bits
	  temp=(i>>j)&3;

	  //move the two bits being examined into the proper reverse positions
	  temp<<=(word-j-2);
  
	  //set the two bits in the rev string to the proper pair of bits
	  revBitStr|=temp;
	 
	  //dynamic_bitset<> bitRepKeyRev(64, revBitStr);
	  //cout<<bitRep<<endl;
	  //cout<<"bit string key is "<<"  "<<revBitStr<<" "<<bitRepKeyRev<<endl;

	}
      
      revBitStr=~revBitStr;
      revBitStr&=reset;

      //std::bitset<64> x(revBitStr);
      //cout<<x<<endl;
      bitTable[i]=revBitStr;
    }

  //return(bitTable);
}

