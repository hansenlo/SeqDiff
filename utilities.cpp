#include "utilities.h"

#include <iostream>
//#include<fstream>
//#include<boost/dynamic_bitset.hpp>

using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::ofstream;
using std::ios;
using std::cerr;
//using boost::dynamic_bitset;
using std::ifstream;

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



string int2String( uint_fast64_t number, int kmerSize)
{
  
  int i, numberNuc;
  string bitString, nucString, nucBit;

  nucString="";
  
  std::bitset<64> bitRep(number);
  
  
  bitString=bitRep.to_string();



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




string bit2String( std::bitset<bitSetSize> &number, int kmerSize)
{
  
  int i, numberNuc;
  string bitString, nucString, nucBit;

  nucString="";
  
  //dynamic_bitset<> bitRep(64, number);
  
  
  bitString=number.to_string();

  //to_string(bitRep, bitString);



  for(i=bitSetSize-(kmerSize*2); i<bitSetSize; i+=2)
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
//*********Also Key word size cannot be bigger than 64 bits!!###########
void createBitRevTable(int word, std::vector< std::bitset<bitSetSize> > &bitTable){

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


//BAD programing practice hard code in the reverse complement must change this if I decide to keep 
void revComplementBitString( std::bitset<bitSetSize> &reversedKey, std::bitset<bitSetSize> &key, std::bitset<bitSetSize> &clearbitWord, std::vector< std::bitset<bitSetSize> > &bitTable, int bitWordSize, int kmerSize)
{

  /*
  cout<<"key is "<<key<<endl;
  cout<<"clearbitWord is "<<clearbitWord<<endl;

  cout<<"integer representation is "<<(key&clearbitWord).to_ulong()<<endl;


  cout<<"key now is "<<(key&clearbitWord)<<endl;
  */

  /*
  int i, shiftRight, shiftLeft;

  shiftLeft=0;

  for(i=bitSetSize-bitWordSize; i>0; i-=bitWordSize)
    {
      shiftRight=i;
      
      reversedKey=
      

      shiftLeft+=bitWordSize;
    }
  */

  //cout<<"size of bit table is "<<sizeof(bitTable[0])<<endl;


  
  //super confusing but fast the bit string size is 192 I find the revComplement of the 
  //first 16 bits then move them to the back of the bit string I then find the reveComplement of
  //bitst 16 through 32 and then move it to its appropriate place at the back of the bit string 
  //and so on and so forth
  reversedKey = (bitTable[(key & clearbitWord).to_ulong()] << 176) | 
    ( bitTable[( (key >> 16) & clearbitWord ).to_ulong()] << 160) | 
    (bitTable[ ((key >> 32) & clearbitWord).to_ulong()] << 144) |
    ( bitTable[( (key >> 48) & clearbitWord ).to_ulong()] << 128) | 
    (bitTable[ ((key >> 64) & clearbitWord).to_ulong()] << 112) |
    ( bitTable[( (key >> 80) & clearbitWord ).to_ulong()] << 96) | 
    (bitTable[ ((key >> 96) & clearbitWord).to_ulong()] << 80) |
    ( bitTable[( (key >> 112) & clearbitWord ).to_ulong()] << 64) | 
    (bitTable[ ((key >> 128) & clearbitWord).to_ulong()] << 48) |
    ( bitTable[( (key >> 144) & clearbitWord ).to_ulong()] << 32) | 
    (bitTable[ ((key >> 160) & clearbitWord).to_ulong()] << 16) |
    (bitTable[ ((key >> 176) & clearbitWord).to_ulong()]); 

  reversedKey=reversedKey>>(bitSetSize-kmerSize*2); //need to move the bits that represent nucleotides back into the right side of the bit string

  

  /*
  reversedKey = (bitTable[key & 0xffff] << 48) | 
    (bitTable[(key >> 16) & 0xffff] << 32) | 
    (bitTable[(key >> 32) & 0xffff] << 16) |
    (bitTable[(key >> 48) & 0xffff]);  

  reversedKey=reversedKey>>(64-kmerSize*2);

  */


  /*orginal 64 bit bit reversal

  reversedKey = (bitTable[(key & clearbitWord).to_ulong()] << 48) | 
    ( bitTable[( (key >> 16) & clearbitWord ).to_ulong()] << 32) | 
    (bitTable[ ((key >> 32) & clearbitWord).to_ulong()] << 16) |
    (bitTable[ ((key >> 48) & clearbitWord).to_ulong()]);  
  */


  //  reversedKey=reversedKey>>(64-kmerSize*2);
  
  

} 
void readInFasta(std::unordered_map<std::string, std::string> &genome, string fileName)
{
  ifstream fasta;
  string line;
  string chr;

  fasta.open(fileName.c_str(), ifstream::in);
  if(!fasta.is_open())
    {
      cerr<<"could not open file "<<fileName<< " check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }

      

      while(fasta.good())  //read in all the clusters
	{
	  getline(fasta, line);
      
	  if(line[0]=='>')
	    {


	      chr=line.substr(1); 

	      cerr<<"chr is "<<chr<<endl;

	      if(chr=="GL000220.1")
		{
		  string temp="1";
		}

	      genome[chr]="";
	    }else
	    {
	      
	      
	      genome[chr]=genome[chr].append(line);

	    }
	  

	}

}

//function to create a bit rev lookup table for use 
//in reverse complementing a nucleotide string
//arguments are the word size you want to reverse i.e. 8 bits 16 bits etc 
//and the lookup table which will be populated
//*******IMPORTANT******* bit reversal is for nucleotide strings!!! don't use bit reversal as universal reversal table*********
//function to create a bit reversal table in the context of a machine word of 64 bits only 
void createBitRevTableMachineWord(int word, std::vector<uint_fast32_t> &bitTable)
{

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

void revComplementMachineWord(uint_fast64_t &key, uint_fast64_t &revKey, int clusterKmerSize, std::vector<uint_fast32_t> &bitTable)
{

  
  revKey = (bitTable[key & 0xffff] << 48) | 
    (bitTable[(key >> 16) & 0xffff] << 32) | 
    (bitTable[(key >> 32) & 0xffff] << 16) |
    (bitTable[(key >> 48) & 0xffff]);  

  revKey=revKey>>(64-clusterKmerSize*2);

}
