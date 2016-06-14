#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include <stdlib.h> 

#include "ReadCluster.h"
#include "utilities.h"
#include "kmerAnalysis.h"



using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::ifstream;
using google::dense_hash_map;
using google::sparse_hash_map;
using std::vector;
using std::ofstream;
//using boost::dynamic_bitset;
using std::bitset;
using std::hash;
using std::unordered_map;
using std::stringstream;
using std::ofstream;


//given a sam file parse it and return it as a matrix every row of the matrix will be an alignment file
//mapCutoff is the mapping quality cutoff
void parseSam(vector<vector<string>> &sam, string samFile, int mapCutoff=0);
void split(vector<string> &tokenizedString,  string str, char delimiter);

//given a vcf file parse it and return it as a matrix every row of the matrix will be a variant call
void parseVcf(vector<vector<string>> &vcf, string vcfFile);


//given a genome a chromosome a position on the genome and a kmer size return through string kmer the sequence extracted from the genome
void getKmer(unordered_map<string, string> &genome, string chr, long position, string &kmer, int kmerSize); 

		    //given a  sequence and a kmerSize convert the sequence to a bit string if this cannot be done then return an empty 
void convertBitString(string &sequence, int kmerSize, bitset<bitSetSize> &bitString);

//will return 0 for heterzygote and 1 for homozygote will return 2 if cannot decide
int callZygosity(dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int indexVcfRecord, vector<vector<string>> &vcf, int kmerSize, unordered_map<string, string> &genome, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, int zygosityCutoff);

//int callZygosity(dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, vector<string> &vcfRecord, int kmerSize, unordered_map<string, string> &genome, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, int zygosityCutoff);



//given a vcf matrix print it to the given file
void printVcf(vector<vector<string>> &vcf, string outputFile);

//this function will treat each variant as an island not considering nearby variants
int checkUniqueness(unordered_map<string, string> &genome, string &altSeq, dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int kmerSize, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, int zygosityCutoff, vector<string> &vcfRecord, bool forward); //bool forward variable means do you start from the begining or end of the alternative Seq
//forward set to true means start from the beginning bool set to false means start from the end 

//this function will group variants together when considering them  
//bool forward variable means do you start from the begining or end of the alternative Seq
//forward set to true means start from the beginning bool set to false means start from the end 
int checkUniquenessGroup(string &altSeq, string &mixedSeq, string &wildType, dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int kmerSize, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, vector<int> &positionsVariant, int zygosityCutoff, bool forward);


//given a vcf record and the set of all vcf records return the alterative sequence that incorportes the variant in the vcf record and all variants within kmerSize 
 //function returns the starting positions in the sequence of the variant of interest order is altSeq, mixedSeq, wildType
vector<int> getAltSeq(string &altSeq, string &mixedSeq, string &wildType, int indexVcfRecord, vector<vector<string>> &vcf, int kmerSize, unordered_map<string, string> &genome);


//given a collection of vcf records a specific index into that collection and a flag indicating whether to look upstream or downstream of the given variant
//return a vector in indexes representing variants withint kmerSize of the given variant
void collectVcf(vector<vector<string>> &vcf, int indexVcfRecord, vector<int> &indexToInclude, bool backward, int kmerSize); 


int main(int argc, char *argv[] )
{

  //vector<vector<string>> samResults;
  
  int kmerSize=45;


  vector<vector<string>> vcfResults;
  

    //cerr<<"starting to parse vcf file "<<endl;
					  //parseSam(samResults, "/data5/SeqDiffResults/Results/Alignment/unique_platinumChr21_plusUnmapped.sam", 20);

				     //parseVcf(vcfResults, "/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/indelCalls.vcf"); 

  //parseVcf(vcfResults, "/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/indelCalls_BWA.vcf"); 

  parseVcf(vcfResults, "/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/allCalls_merged_BWA.vcf"); 

  //parseVcf(vcfResults, "/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/temp.vcf"); 



    cerr<<"Starting to read in genome \n";

	  unordered_map<string, string> genome;

					  
	  readInFasta(genome, "/data/Genomes/human19/allChrhg19InOrder.fa"); //function to read in a multi fasta file and store it in a hash table

	  /*
	  cerr<<"reached this point started alternate Seq generation "<<endl;

	  string altSeq, mixedSeq, wildType;
	  vector<int> positionsVariant;

	  //alt Seq is with all variants included with a mix of reference sequence
	  //mixed Seq is a mix of reference seq and variants the variant of interest is converted to reference not included
	  //wildType is just the reference seq
	  positionsVariant=getAltSeq(altSeq, mixedSeq, wildType, 2922, vcfResults, kmerSize, genome);
	  
	  cerr<<"alternate sequence is "<<altSeq<<" \n";
	  cerr<<"mixedSeq is "<<mixedSeq<<endl;
	  cerr<<"wildTypeSeq is "<<wildType<<endl;

	  for(int z=0; z<positionsVariant.size(); z++)
	    {
	      cerr<<positionsVariant[z]<<"\t";
	    }

	  cerr<<"\n";

	  return(0);
	  */



int i;


dense_hash_map<bitset<bitSetSize>, int, stdHash> uniqueKmers;
bitset<bitSetSize> emptyKey;


   //################NEED TO ADD CODE MAKING SURE kmer size is small enough this never happens##############
emptyKey.set(bitSetSize-1);
uniqueKmers.set_empty_key(emptyKey);


dense_hash_map<bitset<bitSetSize>, int, stdHash> allKmers;
allKmers.set_empty_key(emptyKey);


uniqueKmers.resize(100000000);
allKmers.resize(400000000);   

char continueFlag='0';

 cerr<<"reading in unique kmers "<<endl;

				     //reading in unique kmers
 readUniqueKmers(uniqueKmers, continueFlag, "/data7/PlatinumAlignments/platinumChr21_plusUnmapped_kmerCounts_RefSubtracted_45k_database_human_Readable.dat", kmerSize, 8); //Need to uncomment for code to work

 cerr<<"reading in all kmers "<<endl;

				     //reading in all kmers
  readUniqueKmers(allKmers, continueFlag, "/data7/PlatinumAlignments/platinumChr21_plusUnmapped_kmerCounts_45k_database_human_Readable.dat", kmerSize, 2); //Need to uncomment for code to work




				     /*
string kmer="";

getKmer(genome, "chr21", 15663704-kmerSize-1, kmer, kmerSize); 

cerr<<"kmer now is "<<endl;
cerr<<kmer<<endl;
				     */


bitset<bitSetSize> bitString;


				     cerr<<"starting to process vcf "<<endl;
	    
int bitWord=16;
std::vector< bitset<bitSetSize> > bitTable(pow(2,bitWord));
bitset<bitSetSize> clearbitWord(pow(2,bitWord)-1);
createBitRevTable(bitWord, bitTable);
	    
int call, genotypeIndex, samePositionCtr;
bool samePosition=false;
int sameIndex, j;
samePositionCtr=0;

 string longestRefSeq;
			
	  for(i=0; i<vcfResults.size(); i++)
	      {
		if(vcfResults[i][0][0]=='#') //do not run header vcf lines
		  {
		    continue;
		  }

		if(i % 1000==0)
		  {
		    
		    cerr<<"number of records processed is "<<i<<endl;

		  }
		

		genotypeIndex=vcfResults[i].size()-1;


		if(i<(vcfResults.size()-1))
		  {
		    

		    //if variants are at the same location in the genome call them both heterozygout
		    if( (vcfResults[i][0]==vcfResults[i+1][0] ) && (vcfResults[i][1]==vcfResults[i+1][1]) )
		      {

			//if(vcfResults[i][1]=="9427905")
			//{
	
			    //}
		
			
			samePositionCtr++;
			samePosition=true;
			continue;

		      }else
			 {
			   //if variants all have the same position then determine appropriate 
			   if(samePosition)
			     {

			       sameIndex=i-samePositionCtr;
		
			       int max=0;
			       //finding the longest reference sequence and storing it
			       for(j=sameIndex; j<=(i); j++)
				 {				   
		
				   if(vcfResults[j][3].length()>max)
				     {
				       longestRefSeq=vcfResults[j][3];
				       
				     }

				 }

			 
			       for(j=sameIndex; j<=(i); j++)
				 {				   
				  
				     
				       string altSeq=vcfResults[j][4];
				       string newAltSeq;

				       if(vcfResults[j][3].length() > vcfResults[j][4].length()) //then is a deletion
					 {
					   newAltSeq=longestRefSeq;
					   
					   newAltSeq.erase(1, (vcfResults[j][3].length()-1)); //removing the deletion from the reference sequence
					   
					   //					   if(vcfResults[j][1]=="15663777")
					   //{
					   //cerr<<"newAltSeq in Deletion is "<<newAltSeq<<" alt Seq is "<<vcfResults[j][4]<<" ref Seq is "<<vcfResults[j][3]<<" longestRefSeq is "<<longestRefSeq<<" size to delete is "<<vcfResults[j][3].length()-1<<endl;
					   
					   //}


					 }else if(vcfResults[j][3].length() < vcfResults[j][4].length()) //then is an insertion
					 {
					 
					   newAltSeq=longestRefSeq;
					  
					   string insertedSeq=vcfResults[j][4].substr(1); //getting the sequence that was inserted
					   newAltSeq.insert(1, insertedSeq); //adding inserted sequence
			   
					   //if(vcfResults[j][1]=="15663777")
					   //{
					   //cerr<<"newAltSeq in Insertion is "<<newAltSeq<<" alt Seq is "<<vcfResults[j][4]<<" ref Seq is "<<vcfResults[j][3]<<" longestRefSeq is "<<longestRefSeq<<endl;
					   //}



					 }else if(vcfResults[j][4].length()==vcfResults[j][3].length()) //then is a snp
					 {
					   newAltSeq=longestRefSeq;
					   newAltSeq[1]=vcfResults[j][4][0];

					 }

				       if(j==sameIndex)
					 {
					    //clearing the alternative Sequence
					   vcfResults[sameIndex][4]=newAltSeq;

					 }else
					 {

					   vcfResults[sameIndex][4]=vcfResults[sameIndex][4]+","+newAltSeq; //adding all the alternative allels for variants that map to the same location
					   vcfResults[j][0]="-";
				       
					 }
				     
				 
				 }

			       //setting the genotype and reference sequence for merged variants 
			       vcfResults[sameIndex][genotypeIndex]=std::to_string(samePositionCtr+1)+"|"+std::to_string(samePositionCtr); //setting the right genotype for multiple variants that map to the same location
			       vcfResults[sameIndex][3]=longestRefSeq;

			       samePosition=false;
			       samePositionCtr=0;
			       continue;

			     }else
			     {

			       samePosition=false;
			       samePositionCtr=0;

			     }


			 }	      
		  }

	      

		//cerr<<" i is "<<i<<" "<<vcfResults[i][1]<<endl;
	    //if more than 4 reads contain the reference sequence than do not call it a homozygoty	    
		//call=callZygosity(uniqueKmers, allKmers, vcfResults[i], kmerSize, genome, bitTable, clearbitWord, 4);


		//cerr<<" i is "<<i<<" "<<vcfResults[480][1]<<endl;
	    

		//call=callZygosity(uniqueKmers, allKmers, 480, vcfResults, kmerSize, genome, bitTable, clearbitWord, 4);
		call=callZygosity(uniqueKmers, allKmers, i, vcfResults, kmerSize, genome, bitTable, clearbitWord, 4);

		//return(0);

		//call=callZygosity(uniqueKmers, allKmers, vcfResults[28], kmerSize, genome, bitTable, clearbitWord, 4);

		    
		    if(call==0) //heterozygote
		      {
			
			//cerr<<"variant is a heterozygot! "<<endl;
			
			vcfResults[i][genotypeIndex]="0/1";
		      }else if(call==1)
		      {
		
			//cerr<<"variant is a homozygout! "<<endl;

			vcfResults[i][genotypeIndex]="1/1";
		      }

		    //if the call is 2 than do not do anything keep the orginal call

		    //return(0);
		
		  
	      }

  cerr<<"done processing vcf starting to print vcf file "<<endl;

  printVcf(vcfResults, "test.vcf");



	    /*
convertBitString(kmer, kmerSize, bitString);
cerr<<"converted bit string is "<<bit2String(bitString, kmerSize)<<endl;
	    */

					 /*    
  cerr<<"line of vcf file "<<endl;
  cerr<<"number of vcf records is "<<vcfResults.size()<<endl;

  for(i=0; i<vcfResults[10].size(); i++)
      {
	cerr<<vcfResults[10][i]<<"\t";
      }
    
    cerr<<endl;
    
					 */

return(0);

}


void parseSam(vector<vector<string>> &sam, string samFile, int mapCutoff)
{

  ifstream samRecords;
  string line;
  string chr;
  int ctr;

  samRecords.open(samFile.c_str(), ifstream::in);
  if(!samRecords.is_open())
    {
      cerr<<"could not open file "<<samFile<< " check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }

      
  ctr=0;

      while(samRecords.good())  //read in all the clusters
	{
	  getline(samRecords, line);

	  if(ctr % 1000000==0)
	    {
	      cerr<<"number of sam records is "<<sam.size()<<endl;
	    }

	  ctr++;

	  vector<string> tokenizedString;
	  split(tokenizedString, line, '\t');	  

	  if(tokenizedString.size()>=11) //remove header lines
	    {	  
	  //if sequence maps with mapping quality better than the cutoff than return it
	      if(stoull(tokenizedString[4], NULL, 10) >= mapCutoff)
		{
	
		  sam.push_back(tokenizedString);
	      
		}
	    }


	}

  

}

void split(vector<string> &tokenizedString, string str, char delimiter) {
  //vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  
  while(getline(ss, tok, delimiter)) {
    tokenizedString.push_back(tok);
  }


}


void parseVcf(vector<vector<string>> &vcf, string vcfFile)
{

  ifstream vcfRecords;
  string line;
  string chr;
  int ctr;


  vcfRecords.open(vcfFile.c_str(), ifstream::in);
  if(!vcfRecords.is_open())
    {
      cerr<<"could not open file "<<vcfFile<< " check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }

      
  ctr=0;

      while(vcfRecords.good())  //read in all the clusters
	{
	  getline(vcfRecords, line);

	  if(ctr % 100000==0)
	    {
	      cerr<<"number of vcf records is "<<vcf.size()<<endl;
	    }

	  ctr++;

	  vector<string> tokenizedString;
	  if(!line.empty())
	    {
	      split(tokenizedString, line, '\t');	  
	      vcf.push_back(tokenizedString);

	    }

	  
	}

  
}
		  

void getKmer(unordered_map<string, string> &genome, string chr, long position, string &kmer, int kmerSize)
{
  
  //cerr<<"position is "<<position<<endl;
  //cerr<<"chr is "<<chr<<" length is "<<genome[chr].length()<<" kmer size is "<<kmerSize<<endl;
  //if string is empty grab the sequence from the reference
  //if(kmer.empty())
  //{
      kmer=genome[chr].substr(position, kmerSize);
      //}else
      //{
      //kmer=kmer.append(genome[chr].substr(position, 1)); //append the next base onto the already existing kmer 

      //}

}


 void convertBitString(string &sequence, int kmerSize, bitset<bitSetSize> &bitString)
{
  
  string line, qualityScores, temp, foo, word, kmer;
  char continueFlag, start;
  int ctr, i, revControlValue, controlValue, quality, kmerCtr;
  long totalCtr, limit;
  bool fastq, flag, secondTime, usedRead, qualityReadIn;
  uint_fast64_t index, positionCtr;

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  bitset<bitSetSize> key, AbitString, CbitString, GbitString, TbitString, reversedKey, firstKey, secondKey;
    
  vector<string> splitLine;

  //sparse_hash_map<uint_fast64_t, int, customHash> allKmers;
  //dense_hash_map<uint_fast64_t, int, customHash> denseKmers;

  //denseKmers.set_empty_key(-20);
  //denseKmers.resize(3000000000);

  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
 

  //11 bitstring represents T place into right side of bit string
  TbitString.set((2*kmerSize)-1);
  TbitString.set((2*kmerSize)-2);
 
  //10 bitstring represents G places into right side of bit string
  GbitString.set((2*kmerSize)-1);

  //01 bitstring represents C places into right side
  CbitString.set((2*kmerSize)-2);


  //cerr<<"Sequence in convertBitSTring is "<<sequence<<endl;


 if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      exit(EXIT_FAILURE);
    }


 // stringstream ss(inputFile);

  
  /*
  //building the look up table to reverse nucleotide bit strings
  int bitWord=16;
  std::vector< bitset<bitSetSize> > bitTable(pow(2,bitWord));
  bitset<bitSetSize>clearbitWord(pow(2,bitWord)-1);
  createBitRevTable(bitWord, bitTable);
  */

   //revKey=0;
   bitString.reset();
   positionCtr=0;
   totalCtr=0;
   flag=true;
   ctr=3;
  //arrayIndex=0;
   flag=false;
 	
       //ensuring the line is large enough to fit a kmer in it
      if(sequence.length() < kmerSize)
       	{
	  cerr<<"sequence is to short it is  "<<sequence.length()<<" should be at least "<<kmerSize<<endl;
	  exit(1);
	}
      
  

      //starting a new sequence
      for(i=sequence.length()-1; i>=0; i--)
	{

	   //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(sequence[i]) == std::string::npos) {
		
	    //revKey=0;
	    bitString.reset();
	    sequence=""; //if cannot convert to a bitstring because string contains non nucleotide characters return an empty string
		  //flag=true;
		  continue;
	      }else
		{
		  //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
		  bitString=bitString>>2;
	      
		  switch(sequence[i])
		    {
		    case 'A' :
		      bitString=bitString|AbitString;
		      break;
		    
		    case 'a' :
		      bitString=bitString|AbitString;
		      break;
		   
		    case 'T' :
		      bitString=bitString|TbitString;
		      break;
		    
		    case 't' :
		      bitString=bitString|TbitString;
		      break;
 		    
		    case 'G' :
		      bitString=bitString|GbitString;
		      break;
		    
		    case 'g' :
		      bitString=bitString|GbitString;
		      break;

		    case 'C' :
		      bitString=bitString|CbitString;
		      break;

		    case 'c' :
		      bitString=bitString|CbitString;
		      break;
		      
		    default:
		      cerr<<sequence<<endl;
		      cerr<<"Base is not a nucleotide exiting "<<endl;
		      exit(EXIT_FAILURE);

		    }

		  
		}
  	 

	}
	     
}

int callZygosity(dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int indexVcfRecord, vector<vector<string>> &vcf, int kmerSize, unordered_map<string, string> &genome, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, int zygosityCutoff)
{

  int i, lengthRef, lengthAlt, startRefSeq, endRefSeq;
  string  endPart;

  uint_fast64_t position;

  int bitWord=16;

  /*
  lengthRef=vcfRecord[3].length(); //getting the length of the reference seq. 
  lengthAlt=vcfRecord[4].length(); //getting the length of the alternative seq. 

  position=stoull(vcfRecord[1], NULL, 10);

  
  //cerr<<"position of variant is "<<position<<endl;

  //getting the initial kmer subtracting 1 becase vcf coordinates are 1 based
  getKmer(genome, vcfRecord[0], position-kmerSize-1, altSeq, kmerSize); 


  
  altSeq.append(vcfRecord[4]);
  
  //cerr<<"got alt seq first part "<<endl;

  getKmer(genome, vcfRecord[0], position+lengthRef-1, endPart, kmerSize*2); 

  altSeq.append(endPart);

  //cerr<<"got alt seq end part "<<endl;

  //cerr<<"altSeq is "<<altSeq<<endl;

  //getting the start in the genome of the reference sequence segment that corresponds to the altenative sequence
  startRefSeq=(position+lengthRef-1+2*kmerSize)-altSeq.length();
  endRefSeq=startRefSeq+altSeq.length();
  */



	  string altSeq, mixedSeq, wildType;
	  vector<int> positionsVariant;

	  //alt Seq is with all variants included with a mix of reference sequence
	  //mixed Seq is a mix of reference seq and variants the variant of interest is converted to reference not included
	  //wildType is just the reference seq

	  //cerr<<"before get alternative sequence "<<vcf[indexVcfRecord][1]<<endl;

	  positionsVariant=getAltSeq(altSeq, mixedSeq, wildType, indexVcfRecord, vcf, kmerSize, genome);

	  //return(0);

	  //converting to uppercase
	  //std::transform(altSeq.begin(), altSeq.end(), altSeq.begin(), ::toupper);
	  //std::transform(mixedSeq.begin(), mixedSeq.end(), mixedSeq.begin(), ::toupper);
	  // std::transform(wildType.begin(), wildType.end(), wildType.begin(), ::toupper);
	  
	  //cerr<<"lenght of alt seq is "<<altSeq.length()<<endl;
	  //cerr<<vcf[indexVcfRecord][1]<<endl;
	  //cerr<<"after get AltSeq "<<endl;
	 


	  /*
	  cerr<<"alternate sequence is "<<altSeq<<" \n";
	  cerr<<"mixedSeq is "<<mixedSeq<<endl;
	  cerr<<"wildTypeSeq is "<<wildType<<endl;

	  for(int z=0; z<positionsVariant.size(); z++)
	    {
	      cerr<<positionsVariant[z]<<"\t";
	    }

	  cerr<<"\n";
	  */






  //getKmer(genome, vcfRecord[0], startRefSeq, refSeq, altSeq.length()); 

  //looking for the first unique word







  bool isUnique=false;
  int seqPosition=altSeq.length()-kmerSize;




  //cerr<<"starting from end of alt seqeunce "<<endl;
  //bool forward variable means do you start from the begining or end of the alternative Seq
  bool forward=false;
// int zygosityFlagBackward=checkUniqueness(genome, altSeq, uniqueKmers, allKmers, kmerSize, bitTable, clearbitWord, zygosityCutoff, vcfRecord, forward);


  //cerr<<"before backward check Uniqueness "<<endl;
  int zygosityFlagBackward=checkUniquenessGroup(altSeq, mixedSeq, wildType, uniqueKmers, allKmers, kmerSize, bitTable, clearbitWord, positionsVariant, zygosityCutoff, forward);
  //cerr<<"After backward check Uniqueness "<<endl;


  //cerr<<"end of alt sequence  backward flag is "<<zygosityFlagBackward<<endl;

 forward=true;
  // int zygosityFlagForward=checkUniqueness(genome, altSeq, uniqueKmers, allKmers, kmerSize, bitTable, clearbitWord, zygosityCutoff, vcfRecord, forward);
 int zygosityFlagForward=checkUniquenessGroup(altSeq, mixedSeq, wildType, uniqueKmers, allKmers, kmerSize, bitTable, clearbitWord, positionsVariant, zygosityCutoff, forward);
  
  //
 //cerr<<"starting from beginning of alt sequence forward flag is "<<zygosityFlagForward<<endl;


  if(zygosityFlagBackward > 1 && zygosityFlagForward > 1)
    {
      return(2);
    }

  //if ref kmer is present coming from either the back or front of the variant then call it a heterozyote
  if(zygosityFlagForward==0 || zygosityFlagBackward==0)
    {
      return(0);
    }

if(zygosityFlagForward==1 && zygosityFlagBackward==1)
    {
      return(1);
    }

 

  //cerr<<"starting while loop "<<endl;

  //return(2);

}

void printVcf(vector<vector<string>> &vcf, string outputFile)
{
  ofstream vcfOut(outputFile); //location of where to put file containing contigs

  int i, j;

  //cerr<<"size of matrix is "<<vcf.size()<<endl;

  //int temp=vcf.size()-2;

  //cerr<<"printing last record in file "<<vcf[temp][0]<<endl;
  
  for(i=0; i<vcf.size(); i++)
    {
      if(vcf[i][0]=="-") //if vcf record stars with this symbol means it is one of many that map to same location only print one of those and skip the rest
	{
	  continue; 
	}

      for(j=0; j<(vcf[i].size()-1); j++)
	{
	  vcfOut<<vcf[i][j]<<"\t";
	  
	  //cerr<<vcf[i][j]<<"\t";
	  
	}

      //cerr<<"\n";
      int lastColumn=vcf[i].size()-1;
      vcfOut<<vcf[i][lastColumn]<<"\n";
      //cerr<<vcf[i][lastColumn]<<"\n";

    }
  

}


int checkUniqueness(unordered_map<string, string> &genome, string &altSeq, dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int kmerSize, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, int zygosityCutoff, vector<string> &vcfRecord, bool forward) //bool forward variable means do you start from the begining or end of the alternative Seq
{

  int bitWord=16;

  int i, lengthRef, lengthAlt, startRefSeq, endRefSeq;
  uint_fast64_t position;
  int seqPosition, seqCtr;

  lengthRef=vcfRecord[3].length(); //getting the length of the reference seq. 
  lengthAlt=vcfRecord[4].length(); //getting the length of the alternative seq. 

  position=stoull(vcfRecord[1], NULL, 10);

  if(!forward) //starting from the end of the alternative sequence
    {
    
      //getting the start of the reference sequence segment that corresponds to the altenative sequence
      startRefSeq=(position+lengthRef-1+2*kmerSize)-altSeq.length();
      endRefSeq=startRefSeq+altSeq.length();
      
      seqPosition=altSeq.length()-kmerSize;

    }else
    {
      startRefSeq=position-kmerSize-1;
      seqPosition=0;

    }

    bool isUnique=false;
    seqCtr=altSeq.length()-kmerSize;

    //    cerr<<"seqCtr is "<<seqCtr<<endl;

  while(!isUnique && seqCtr>0)
    {
     
      bitset<bitSetSize> bitString;
      bitset<bitSetSize> reversedString; //bit string to hold the reverse complement

      //cerr<<"seq Pos is "<<seqPosition<<endl;
      string kmerToTest=altSeq.substr(seqPosition, kmerSize);

    
      convertBitString(kmerToTest, kmerSize, bitString); //converting to a bitstring
      
      if(!forward)
	{
	  seqPosition--;
	}else
	{
	  seqPosition++;
	}

      seqCtr--; //count of how many kmers I have already tested

      if(kmerToTest.empty()) //if failed to convert to a bit string then do not test in hash table
	{
	  continue;
	}

      revComplementBitString(reversedString, bitString, clearbitWord, bitTable, bitWord, kmerSize);
      
      if(uniqueKmers.count(bitString)>0 || uniqueKmers.count(reversedString)>0)
	{
	  isUnique=true;
	  
	  //cerr<<"unique bit string is "<<bit2String(bitString, kmerSize)<<" with a count of "<<uniqueKmers[bitString]<<endl;


	  //getting the reference kmer corresponding to the variant
	  string refSeq;

	  if(!forward)
	    {
	      getKmer(genome, vcfRecord[0], startRefSeq+seqPosition+1, refSeq, kmerSize); 

	    }else
	    {
	      getKmer(genome, vcfRecord[0], startRefSeq+seqPosition-1, refSeq, kmerSize); 
	    }

	 
	  //cerr<<"ref equnivalent of unique word is "<<refSeq<<" seqPosition is "<<seqPosition<<" startRefSeq is "<<startRefSeq<<endl;

	  

	  //get the reference sequence corresponding to the variant
	  bitString.reset();
	  reversedString.reset();
	  convertBitString(refSeq, kmerSize, bitString); //converting to bit string the reference word

	  //reverseComplementing the reference kmer
	  revComplementBitString(reversedString, bitString, clearbitWord, bitTable, bitWord, kmerSize);
	  
	  // int count=allKmers.count(bitString)+allKmers.count(reversedString);
	  //	    cerr<<"reference bit string is "<<bit2String(bitString, kmerSize)<<" with a count of for reverse "<<allKmers[reversedString]<<" with a count of for + strand "<<allKmers[bitString]<<" sum is "<<temp<<endl;
	  //	    cerr<<"zygosityCutoff is "<<zygosityCutoff<<endl;

	  
	  //if( (allKmers.count(bitString)+allKmers.count(reversedString)) > 0 )
	    
	      int sum=0;
	      if(allKmers.count(bitString)>0 && allKmers.count(reversedString)>0)
		{
		  sum=allKmers[bitString]+allKmers[reversedString];
		
		}else if(allKmers.count(bitString)>0)
		{

		  sum=allKmers[bitString];

		}else if (allKmers.count(reversedString)>0)
		{
		  sum=allKmers[reversedString];
		
		}


	      if(sum > zygosityCutoff)
		{

		  //cerr<<"variant is a heterozygoute!! "<<endl;
		  //cerr<<"reference bit string is "<<bit2String(bitString, kmerSize)<<" with a count of "<<allKmers[bitString]<<endl;

		  return(0);
		}else
		{
	      
		  //cerr<<"variant is a homozygoute! "<<endl;
	    
		  return(1);
		}

	    
	}

    }

  return(2);



}      	  

//function returns the starting positions in the sequence of the variant of interest order is altSeq, mixedSeq, wildType
vector<int> getAltSeq(string &altSeq, string &mixedSeq, string &wildType, int indexVcfRecord, vector<vector<string>> &vcf, int kmerSize, unordered_map<string, string> &genome)
{
  
  int lengthRef; //=vcf[indexVcfRecord][3].length(); //getting the length of the reference seq. 
  int lengthAlt; //=vcf[indexVcfRecord][4].length(); //getting the length of the alternative seq. 
  vector<int> postionsVariant; //starting position in the sequence of the variant of interest
  int posAlt, posMixed, posWild;

  vector<int> indexToInclude; //vector to hold the indexes of the vcf records to incoporate into the alternate sequence


  //cerr<<"before forward collectVcf "<<endl;

  //get vcf records within kmerSize upstream of the highlighted variant
  collectVcf(vcf,  indexVcfRecord, indexToInclude, true, kmerSize);


  //cerr<<"before backward collectVcf "<<endl;
  //get vcf records within kmerSize downstream of the highlighted variant
  collectVcf(vcf,  indexVcfRecord, indexToInclude, false, kmerSize);
  
  //cerr<<"after collectVcf "<<endl;

  //sorting the vector
  sort(indexToInclude.begin(), indexToInclude.end());


  /*
  int j;

  for(j=0; j<indexToInclude.size(); j++)
    {
      cerr<<indexToInclude[j]<<endl;
    }
  */
  
  //cerr<<"Position is "<<vcf[indexToInclude[0]][1]<<endl;  


/*
  //beginning


//getting the initial kmer subtracting 1 becase vcf coordinates are 1 based
  getKmer(genome, vcfRecord[0], position-kmerSize-1, altSeq, kmerSize); 

  
  altSeq.append(vcfRecord[4]);
  
  //cerr<<"got alt seq first part "<<endl;

  getKmer(genome, vcfRecord[0], position+lengthRef-1, endPart, kmerSize*2); 

  altSeq.append(endPart);

  //cerr<<"got alt seq end part "<<endl;

  //cerr<<"altSeq is "<<altSeq<<endl;

  //getting the start in the genome of the reference sequence segment that corresponds to the altenative sequence
  startRefSeq=(position+lengthRef-1+2*kmerSize)-altSeq.length();
  endRefSeq=startRefSeq+altSeq.length();

  //end 
  */






  //cerr<<indexToInclude[0]<<"\t"<<indexToInclude[1]<<endl;
  
  //cerr<<"position of variant is "<<position<<endl;

  int i;
  //string altSeq;
  long positionPrev, positionCurrent, lengthPrevRef, lengthCurrRef;
  string endPart;

  positionCurrent=stoull(vcf[indexToInclude[0]][1], NULL, 10);

  //getting the initial sequence
  getKmer(genome, vcf[indexToInclude[0]][0], positionCurrent-(3*kmerSize)-1, altSeq, 3*kmerSize);
  mixedSeq=altSeq; //getting the wild type sequence only 
  wildType=altSeq; //getting the witd type sequence only


  if(indexToInclude[0]!=indexVcfRecord)
    {
     
      mixedSeq.append(vcf[indexToInclude[0]][4]);   

    }else
    {

      posAlt=altSeq.length();
      posMixed=altSeq.length();
      posWild=altSeq.length();
      
      mixedSeq.append(vcf[indexToInclude[0]][3]); //append the reference sequence
     
    }

  altSeq.append(vcf[indexToInclude[0]][4]); //adding the alternant sequence from the variant
  wildType.append(vcf[indexToInclude[0]][3]);
 

   lengthCurrRef=vcf[indexToInclude[0]][3].length();

   //cerr<<"adding to alt seq this  "<<vcf[indexToInclude[0]][4]<<endl;
 
   //cerr<<"alt seq at this point is "<<altSeq<<endl;
 
  for(i=1; i<indexToInclude.size(); i++)
    {

      positionPrev=stoull(vcf[indexToInclude[i-1]][1], NULL, 10);
      positionCurrent=stoull(vcf[indexToInclude[i]][1], NULL, 10);

     
      lengthPrevRef=vcf[indexToInclude[i-1]][3].length();
      
      if((positionCurrent-positionPrev-lengthPrevRef)<0) //if two variants have the same position ignore the second and use the first
	{
	  continue;
	
	}



      lengthCurrRef=vcf[indexToInclude[i]][3].length();
      

      //cerr<<"length is "<<positionCurrent-positionPrev-lengthPrevRef<<endl;

      
      string section;
      getKmer(genome, vcf[indexToInclude[i]][0], positionPrev+lengthPrevRef-1, section, positionCurrent-positionPrev-lengthPrevRef); 
      
      altSeq.append(section);
      mixedSeq.append(section);
      wildType.append(section);

      //cerr<<"after appending the section on "<<altSeq<<" section is "<<section<<endl;

      if(indexToInclude[i]!=indexVcfRecord)
	{

	  mixedSeq.append(vcf[indexToInclude[i]][4]); //append alternate sequence if not the variant of interest 

	}else
	{    

	  posAlt=altSeq.length();
	  posMixed=mixedSeq.length();
	  posWild=wildType.length();
	  mixedSeq.append(vcf[indexToInclude[i]][3]); //append the reference sequence    
	}

      
      wildType.append(vcf[indexToInclude[i]][3]);      
      altSeq.append(vcf[indexToInclude[i]][4]);

      //cerr<<"just appended "<<vcf[indexToInclude[i]][4]<<endl;
      
 
    }

  //attaching the last little bit of reference sequence

  //cerr<<"wild sequence is "<<wildType<<endl;


  getKmer(genome, vcf[indexToInclude[0]][0], positionCurrent+lengthCurrRef-1, endPart, kmerSize*3); 

  //cerr<<"now appending "<<endPart<<endl;
  
  altSeq.append(endPart);
  mixedSeq.append(endPart);
  wildType.append(endPart);

  postionsVariant.push_back(posAlt);
  postionsVariant.push_back(posMixed);
  postionsVariant.push_back(posWild);



  //cerr<<"alt seq is "<<altSeq<<endl;

      /*
  //cerr<<"got alt seq first part "<<endl;

  getKmer(genome, vcfRecord[0], position+lengthRef-1, endPart, kmerSize*2); 

  altSeq.append(endPart);

  //cerr<<"got alt seq end part "<<endl;

  //cerr<<"altSeq is "<<altSeq<<endl;

  //getting the start in the genome of the reference sequence segment that corresponds to the altenative sequence
  startRefSeq=(position+lengthRef-1+2*kmerSize)-altSeq.length();
  endRefSeq=startRefSeq+altSeq.length();
      */


  return(postionsVariant);  

}
void collectVcf(vector<vector<string>> &vcf, int indexVcfRecord, vector<int> &indexToInclude, bool backward, int kmerSize)
{

  long position=stoull(vcf[indexVcfRecord][1], NULL, 10);

  int currentIndex, currentPosition;
  
  double sizeVcf=vcf.size(); //number of vcf records

  //cerr<<"before if statement collectVcf "<<endl;

  if(backward){

    currentPosition=stoull(vcf[indexVcfRecord][1], NULL, 10);
    currentIndex=indexVcfRecord;

  }else
    {
      if((indexVcfRecord+1)<sizeVcf)
	{

	  currentPosition=stoull(vcf[indexVcfRecord+1][1], NULL, 10);
	  currentIndex=indexVcfRecord+1;
	}else
	{

	  return;

	}
    
    }


  //cerr<<"after if statement collectVcf "<<endl;
    
  
  //getting all vcf records with kmerSize of the variant of interest that occur before the variant
  while( (vcf[indexVcfRecord][0]==vcf[currentIndex][0]) && (abs(position-currentPosition) < kmerSize))
    {      
      //cerr<<"adding index "<<currentIndex<<endl;
      ///cerr<<"size of vcf file is "<<sizeVcf<<endl;

      indexToInclude.push_back(currentIndex);
      
      
      if((currentIndex-1)>=0)
	{
	  if(backward)
	    {
	      currentIndex--;
	    }else
	    {
	      currentIndex++;
	    }


	}else
	{
	  break;
	}
     
      if(currentIndex==sizeVcf)
	{
	  break;
	}


      if(vcf[currentIndex][0][0]!='#')//if the first character is a pound symbol than have reached the comments part of vcf file 
	{

	  currentPosition=stoull(vcf[currentIndex][1], NULL, 10);
	}else
	{
	  break;
	}
      
	  //cerr<<"currentPosition is "<<currentPosition<<endl;
    }




}

//forward set to true means start from the beginning bool set to false means start from the end 

int checkUniquenessGroup(string &altSeq, string &mixedSeq, string &wildType, dense_hash_map<bitset<bitSetSize>, int, stdHash> &uniqueKmers, dense_hash_map<bitset<bitSetSize>, int, stdHash> &allKmers, int kmerSize, std::vector< bitset<bitSetSize> > &bitTable, bitset<bitSetSize> &clearbitWord, vector<int> &positionsVariant, int zygosityCutoff, bool forward)
{

  int bitWord=16;

  int i, lengthRef, lengthAlt, startRefSeq, endRefSeq;
  uint_fast64_t position;
  int seqPosition, seqCtr;

  bool isUnique=false;
    
  seqCtr=0;


    //    cerr<<"seqCtr is "<<seqCtr<<endl;

    if(forward)
      {
	seqPosition=positionsVariant[0]-kmerSize;
      }else
      {
	seqPosition=positionsVariant[0]+1;
      }


    
    while(!isUnique && seqPosition>=0 && (seqPosition <= (altSeq.length()-kmerSize)) )
      {
     
	bitset<bitSetSize> bitString;
	bitset<bitSetSize> reversedString; //bit string to hold the reverse complement for the wild type
	bitset<bitSetSize> bitStringMixed;
	bitset<bitSetSize> reversedStringMixed; //bit string to hold the reverse complement for the mixed

	
	//cerr<<"seqPosition is "<<seqPosition<<" alt Seq is "<<altSeq<<endl;
	//return(0);

	//cerr<<"seq Pos in alt seq is  "<<seqPosition<<endl;
	string kmerToTest=altSeq.substr(seqPosition, kmerSize);
    
	//cerr<<"kmerTo Test is "<<kmerToTest<<endl;

	convertBitString(kmerToTest, kmerSize, bitString); //converting to a bitstring
	
	//cerr<<"after convert to Bit String "<<endl;

	if(kmerToTest.empty()) //if failed to convert to a bit string then do not test in hash table
	  {
	    //cerr<<"inside railed to convert ***** "<<endl;

	    if(!forward)
	      {
		seqPosition--;
		seqCtr--; //offset from the starting position of the variant

	      }else
	      {
		seqPosition++;
		seqCtr++; //offset from the starting position of the variant
	      }



	    continue;
	  }

	revComplementBitString(reversedString, bitString, clearbitWord, bitTable, bitWord, kmerSize);
             
	if(uniqueKmers.count(bitString)>0 || uniqueKmers.count(reversedString)>0)
	  {
	    isUnique=true;
	  
	  //cerr<<"unique bit string is "<<bit2String(bitString, kmerSize)<<" with a count of "<<uniqueKmers[bitString]<<endl;

	  //getting the reference kmer corresponding to the variant and the mixed kmer
	    string refKmer, mixedKmer;

	    //cerr<<"In unique Condition before kmer extraction from mixed and ref seqCtr is "<<seqCtr<<endl;

	    //cerr<<"unique sequence is "<<kmerToTest<<endl;

	    if(forward)
	      {
		refKmer=wildType.substr(positionsVariant[2]-kmerSize+seqCtr+1, kmerSize); 
		mixedKmer=mixedSeq.substr(positionsVariant[1]-kmerSize+seqCtr+1, kmerSize); 
	      

	      }else
	      {

		//cerr<<"Position in wild type is "<<positionsVariant[2]+seqCtr+1<<" size of wild type is "<<wildType.length()<<endl;
		//cerr<<"Position in mixed type is "<<positionsVariant[1]+seqCtr+1<<" size of mixed is "<<mixedSeq.length()<<endl;
		
		refKmer=wildType.substr(positionsVariant[2]+seqCtr-1, kmerSize); 
		mixedKmer=mixedSeq.substr(positionsVariant[1]+seqCtr-1, kmerSize); 
       
	      }

	    //cerr<<"wild Type is "<<refKmer<<endl;
	    //cerr<<"mixed is "<<mixedKmer<<endl;


	    // cerr<<"In unique Condition after kmer extraction from mixed and ref"<<endl;
	 
	    
	    //cerr<<"positionsVariant is "<<positionsVariant[2]<<" seq Ctr is "<<seqCtr<<" size of sequence is "<<refKmer.size()<<endl;



	  //cerr<<"ref equnivalent of unique word is "<<refSeq<<" seqPosition is "<<seqPosition<<" startRefSeq is "<<startRefSeq<<endl;

	  //get the reference sequence corresponding to the variant
	  bitString.reset();
	  reversedString.reset();
	  bitStringMixed.reset();
	  reversedStringMixed.reset();

	  //cerr<<"refKmer is "<<refKmer<<endl;
	  convertBitString(refKmer, kmerSize, bitString); //converting to bit string the reference word

	  //reverseComplementing the reference kmer
	  revComplementBitString(reversedString, bitString, clearbitWord, bitTable, bitWord, kmerSize);

	  //cerr<<"mixedKmer is "<<mixedKmer<<endl;
	  convertBitString(mixedKmer, kmerSize, bitStringMixed); //converting to bit string the mixed word
	  //reverseComplementing the mixed kmer
	  revComplementBitString(reversedStringMixed, bitString, clearbitWord, bitTable, bitWord, kmerSize);


	  
	  // int count=allKmers.count(bitString)+allKmers.count(reversedString);
	  //	    cerr<<"reference bit string is "<<bit2String(bitString, kmerSize)<<" with a count of for reverse "<<allKmers[reversedString]<<" with a count of for + strand "<<allKmers[bitString]<<" sum is "<<temp<<endl;
	  //	    cerr<<"zygosityCutoff is "<<zygosityCutoff<<endl;

	  
	  //if( (allKmers.count(bitString)+allKmers.count(reversedString)) > 0 )
	    
	      int sum=0;
	      if(allKmers.count(bitString)>0 && allKmers.count(reversedString)>0 && allKmers.count(bitStringMixed)>0 && allKmers.count(reversedStringMixed)>0)
		{
		  sum=allKmers[bitString]+allKmers[reversedString]+allKmers[bitStringMixed]+allKmers[reversedStringMixed];
		
		}else if(allKmers.count(bitString)>0)
		{

		  sum=allKmers[bitString];

		}else if (allKmers.count(reversedString)>0)
		{
		  sum=allKmers[reversedString];
		
		}else if (allKmers.count(bitStringMixed)>0)
		{
		  sum=allKmers[bitStringMixed];
		
		}else if (allKmers.count(reversedStringMixed)>0)
		{
		  sum=allKmers[reversedStringMixed];
		}


	      if(sum > zygosityCutoff)
		{

		  //cerr<<"variant is a heterozygoute!! "<<endl;
		  //cerr<<"reference bit string is "<<bit2String(bitString, kmerSize)<<" with a count of "<<allKmers[bitString]<<endl;

		  return(0);
		}else
		{
	      
		  //cerr<<"variant is a homozygoute! "<<endl;
	    
		  return(1);
		}

	    
	  }


	if(!forward)
	  {
	    seqPosition--;
	    seqCtr--; //offset from the starting position of the variant

	  }else
	  {
	    seqPosition++;
	    seqCtr++; //offset from the starting position of the variant
	  }
	



      }

  return(2);



}










