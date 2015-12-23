#include "kmerAnalysis.h"
#include "utilities.h"


using google::dense_hash_map;
using google::sparse_hash_map;
using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::ofstream;
using std::ios;
using std::ifstream;
using std::cerr;
using std::istringstream;
using std::to_string;
using std::list;
using std::make_shared;

inline void countKmers(sparse_hash_map<uint_fast64_t, int, customHash> &controlKmers, char nextLineFlag, string inputFile, int kmerSize)
{
  ifstream seqFile;
  string line;
  char continueFlag, start;
  int ctr, i;
  long totalCtr;
  bool fastq, flag;
  uint_fast64_t index;

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  uint_fast64_t key, positionCtr, AbitString, CbitString, GbitString, TbitString; 

  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));
  
  //places into right side of bit string
 GbitString=0*pow(2, (2*kmerSize))+1*pow(2, ((2*kmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

 //place into right side of bitstring
 CbitString=1*pow(2, (2*kmerSize)-2)+0*pow(2, ((2*kmerSize)));


 continueFlag=nextLineFlag;

 if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      exit(EXIT_FAILURE);
    }


  seqFile.open(inputFile.c_str());
   

  if(!seqFile.is_open())
    {
      cerr<<"could not open file check to see if exists "<<" file is "<<inputFile<<endl;
      exit(EXIT_FAILURE);
    }
  
  
   getline(seqFile, line);
   start=toupper(line[0]);
   
   fastq=false;
   if(start!='>')
     {
       fastq=true; 
     }
     

   


   //###############################code blocking dealing with sequncing libraries##################

   //revKey=0;
  key=0;
  positionCtr=0;
  totalCtr=0;
  //flag=true;
  ctr=3;
  //arrayIndex=0;
  //flag=false;

  
  while(seqFile.good())
    {
     
      totalCtr++;
      ctr++;
      //read in every line of the sequence file
      getline(seqFile, line);
       
      //cout<<totalCtr<<endl;
      if(totalCtr % 1000000==0)
	{
	  
	  cerr<<"lines processed from input file "<<totalCtr<<endl;
	
	  //return;

	  //return(0);
	 
	}
   	
      if(fastq)
	{
	 
	  //very very odd if ctr is set to 0 
	  //code slows down dramatically I have no 
	  //idea why
	  if(ctr<3)
	    {
	      continue;
	    }else
	    {
	      ctr=-1;
	    }

	  
	}


      if(!fastq)
	{
	  //checking to see if line starts with multiple Ns if so probably
	  //repetive region and needs to be skipped
	  if(line[0]=='N' &&  line[1]=='N')
	    {
	      continue;
	    }
	  
	}

      //cout<<line<<endl;
	
       //ensuring the line is large enough to fit a kmer in it
      if(!(line.length()>= kmerSize))
	{
	  continue;
	}
      
      
      //new way of comparision looking to make sure line starts with a valid
      //sequence char
      if (validChar.find(line[0]) == std::string::npos) { 
	continue;
      }

      /* old way of doing comparision
      //check to make sure reading in line that has sequence info
      //cout<<"first char line is"<<line[0]<<endl;
      //start=toupper(line[0]);
      if(line[0]!=('A') && line[0]!=('T') && line[0]!=('C') && line[0]!=('G') && line[0]!=('N')&&line[0]!=('a') && line[0]!=('t') && line[0]!=('c') && line[0]!=('g'))
	{
	  //cout<<"testing "<<line<<endl;
	  continue;
	}
      */

      //if sequence does not continue on the next line reset all values to zero 
      //since sequence on next line is unrelated to current sequence
      if(continueFlag=='0')
	{
	  //revKey=0;
	  key=0;
	  positionCtr=0;
	  //flag=true;
	  
	}
      

      //starting a new sequence
      for(i=line.length()-1; i>=0; i--)
	{
	     
	  /*
	      //next character is something besides sequence 
	      if(nuc.find_first_not_of("ACGTacgt", 0)!=string::npos)
		{
		  //reset bit string to zero
		  revKey=0;
		  key=0;
		  positionCtr=0;
		  flag=true;
		  continue;
		 
		}
	  */


	  /*
	  if(line[i]!=('A') && line[i]!=('T') && line[i]!=('C') && line[i]!=('G') && line[i]!=('a') && line[i]!=('t') && line[i]!=('c') && line[i]!=('g'))
	    {
	  //cout<<"testing "<<line<<endl;
	      revKey=0;
	      key=0;
	      positionCtr=0;
	      flag=true;
	      continue;
	    }
	  */

	  
	  //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(line[i]) == std::string::npos) {
		
	    //revKey=0;
		  key=0;
		  positionCtr=0;
		  //flag=true;
		  continue;
	      }else
		{
		  //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
		  key=key>>2;
	      
		  switch(line[i])
		    {
		    case 'A' :
		      key=key|0;
		      break;
		    
		    case 'a' :
		      key=key|0;
		      break;
		   
		    case 'T' :
		      key=key|TbitString;
		      break;
		    
		    case 't' :
		      key=key|TbitString;
		      break;
		    
		    case 'G' :
		      key=key|GbitString;
		      break;
		    
		    case 'g' :
		      key=key|GbitString;
		      break;

		    case 'C' :
		      key=key|CbitString;
		      break;

		    case 'c' :
		      key=key|CbitString;
		      break;
		      
		    default:
		      cerr<<"Base is not a nucleotide exiting "<<endl;
		      exit(EXIT_FAILURE);

		    }


		  //dynamic_bitset<> bitRep2(64, key);
		  //cout<<" key is "<<key<<endl;
		  //cout<<bitRep2<<endl;
		  //cout<<"kmer is "<<bit2String(key, kmerSize)<<endl;
		  //cout<<"line substring is "<<line.substr(i-kmerSize+1, kmerSize)<<endl;
		  //cout<<"line is "<<line<<endl;


		  positionCtr++;		  
		  
		}
	  /*
	  dynamic_bitset<> bitRep2(64, GbitString);
	  cout<<" key is "<<key<<endl;
	  cout<<bitRep2<<endl;
	  cout<<"nucleotide is "<<line[i]<<endl;
	  */
	      //not yet reached a full kmer in size
	      if((positionCtr<kmerSize))
		{
		  continue;  
		}
	      
	      /*
	      if(flag)
		{
		  flag=false;
		  //key=key<<(64-(kmerSize*2));

		}
	      */




	      //key/64 is integer division finding which char is going to be modified
	      //key%64 tells me which bit in the given character is getting changed
	      //controlCtr1[key/64]|=(1ull<<(key%64));
	      controlKmers[key]=controlKmers[key]+1;

	      //cout<<"key is "<<key<<endl;



	      //dynamic_bitset<> bitRep2(64, controlCtr1[key/64]);
	      //dynamic_bitset<> bitRep(64, (1ull<<(key%64)));
	      

	      //cout<<" key is "<<key<<endl;
	      //cout<<bitRep2<<endl;
	      //cout<<"kmer is "<<bit2String(key, kmerSize)<<endl;
	      //cout<<"line substring is "<<line.substr(i, kmerSize)<<endl;
	      //cout<<"line is "<<line<<endl;

	      //cout<<"key in count control kmers  is "<<key<<" "<<bit2String(key, kmerSize)<<" bit to modify is  "<<key%64<<" "<<bitRep2<<endl;
	      //cout<<"left shifted "<<bitRep<<endl;

	      /*
	      dynamic_bitset<> bitRep2(8, key/8);
	      cout<<" key is "<<key<<" index is "<<key/8<<" bit to modify is "<<key%8<<endl;
	      cout<<bitRep2<<endl;
	      */

//return(1);
				       
	      //exit(1);

	}

     
    }


  seqFile.close();

  //return(controlCtr1);
}


//void countUniqueKmers(vector <uint_fast64_t> &controlCtr1, vector <uint_fast64_t> &controlCtr2,   sparse_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, string inputFile, int kmerSize)
//void countUniqueKmers(vector <uint_fast64_t> &controlCtr1, vector <uint_fast64_t> &controlCtr2,   dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers,  char nextLineFlag, string inputFile, int kmerSize)
void countUniqueKmers(sparse_hash_map<uint_fast64_t, int, customHash>  &controlKmers, sparse_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, string inputFile, int kmerSize)
{
  ifstream seqFile;
  string line, qualityScores, temp, foo;
  char continueFlag, start;
  int ctr, i, revControlValue, controlValue, quality;
  long totalCtr;
  bool fastq, flag, secondTime, usedRead, qualityReadIn;
  uint_fast64_t index, reversedKey, firstKey, secondKey;

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  uint_fast64_t key, positionCtr, AbitString, CbitString, GbitString, TbitString; 

  //sparse_hash_map<uint_fast64_t, int, customHash> allKmers;
  //dense_hash_map<uint_fast64_t, int, customHash> denseKmers;

  //denseKmers.set_empty_key(-1);
  //denseKmers.resize(3000000000);

  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));
  
  //places into right side of bit string
 GbitString=0*pow(2, (2*kmerSize))+1*pow(2, ((2*kmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

 //place into right side of bitstring
 CbitString=1*pow(2, (2*kmerSize)-2)+0*pow(2, ((2*kmerSize)));


 continueFlag=nextLineFlag;

 if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      exit(EXIT_FAILURE);
    }


  seqFile.open(inputFile.c_str());
   

  if(!seqFile.is_open())
    {
      cerr<<"could not open file check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }
  
  
   getline(seqFile, line);
   start=toupper(line[0]);
   
   fastq=false;
   if(start!='>')
     {
       fastq=true; 
     }
     


  //building the look up table to reverse nucleotide bit strings
   int bitWord=16;
   std::vector<uint_fast32_t> bitTable(pow(2,bitWord));

   createBitRevTable(bitWord, bitTable);


   //revKey=0;
  key=0;
  positionCtr=0;
  totalCtr=0;
  flag=true;
  ctr=3;
  //arrayIndex=0;
  flag=false;
  firstKey=0;
  secondKey=0;
  while(seqFile.good())
    {

      totalCtr++;
      ctr++;
      //read in every line of the sequence file
      getline(seqFile, line);
       
      //cout<<totalCtr<<endl;
      if(totalCtr % 10000000==0)
	{
	  
	  cerr<<"lines processed from input file "<<totalCtr<<" size of hash table is "<<uniqueKmers.size()<<endl;
	  
	
	  //return;

	  //return(0);
	 
	}
   	
      if(fastq)
	{
	  
	  //very very odd if ctr is set to 0 
	  //code slows down dramatically I have no 
	  //idea why
	  if(ctr<3)
	    {
	      continue;
	    }else
	    {
	      ctr=-1;
	    }
	  
	  
	}
	
       //ensuring the line is large enough to fit a kmer in it
      if(!(line.length()>= kmerSize))
	{
	  continue;
	}
      
      
      //new way of comparision looking to make sure line starts with a valid
      //sequence char
      if (validChar.find(line[0]) == std::string::npos) { 
	continue;
      }

   



       //if sequence does not continue on the next line reset all values to zero 
      //since sequence on next line is unrelated to current sequence
      if(continueFlag=='0')
	{
	  //revKey=0;
	  key=0;
	  positionCtr=0;
	  flag=true;
	  
	}
      
      //cout<<line<<endl;

      qualityReadIn=false;
      secondTime=false;
      usedRead=false;
      firstKey=0; 
      secondKey=0;
      //starting a new sequence
      for(i=line.length()-1; i>=0; i--)
	{

	   //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(line[i]) == std::string::npos) {
		
	    //revKey=0;
		  key=0;
		  positionCtr=0;
		  //flag=true;
		  continue;
	      }else
		{
		  //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
		  key=key>>2;
	      
		  switch(line[i])
		    {
		    case 'A' :
		      key=key|0;
		      break;
		    
		    case 'a' :
		      key=key|0;
		      break;
		   
		    case 'T' :
		      key=key|TbitString;
		      break;
		    
		    case 't' :
		      key=key|TbitString;
		      break;
 		    
		    case 'G' :
		      key=key|GbitString;
		      break;
		    
		    case 'g' :
		      key=key|GbitString;
		      break;

		    case 'C' :
		      key=key|CbitString;
		      break;

		    case 'c' :
		      key=key|CbitString;
		      break;
		      
		    default:
		      cerr<<line<<endl;
		      cerr<<"Base is not a nucleotide exiting "<<endl;
		      exit(EXIT_FAILURE);

		    }

		  positionCtr++;		  
		  
		}

	   //not yet reached a full kmer in size
	      if((positionCtr<kmerSize))
		{
		  continue;  
		}
	      
	
	      //cout.setf(ios::fixed);
	      //cout<<line.substr(i-kmerSize+1, kmerSize)<<endl;
	      //dynamic_bitset<> bitRep(64, key);
	      //cout<<bitRep<<endl;
	      //cout<<line.substr(i, kmerSize)<<"  "<<key<<" "<<bitRep<<endl;


	      //dynamic_bitset<> bitRep(64, key);
	      //cout<<bitRep<<endl;
	      //cout<<line.substr(i, kmerSize)<<"  before reversal "<<key<<" "<<bitRep<<endl;
	      //cout<<bit2String(key, kmerSize)<<endl;
	    
	      
	      //revKey=revKey&revSubString;

	      //16 bit reverse complement 
	      reversedKey = (bitTable[key & 0xffff] << 48) | 
		(bitTable[(key >> 16) & 0xffff] << 32) | 
		(bitTable[(key >> 32) & 0xffff] << 16) |
		(bitTable[(key >> 48) & 0xffff]);  

	      reversedKey=reversedKey>>(64-kmerSize*2);

	      
	      //dynamic_bitset<> bitRepRev(64, reversedKey);
	      //cout<<bitRepRev<<endl;
	      //cout<<line.substr(i, kmerSize)<<"  "<<reversedKey<<" "<<bitRepRev<<endl;
	      //cout<<bit2String(reversedKey, kmerSize)<<endl;
	      //exit(1);

	      //return;

	      //foo=bit2String(key, kmerSize);

	      //cout<<foo<<endl;
	      

	      //	      cout<<foo<<endl;


	      //if(foo=="TTCTCCCTCCTTGTCGAGTTTCCAATAAAAG")
	      //{
	      //cerr<<"found string "<<endl;
		  
	      //}
	     



	      /*
	      revCtrHash=controlKmers.count(reversedKey);
	      ctrHash=controlKmers.count(key);
	      */
	      // ctrHash=1;
	      //revCtrHash=1;
	      
	     
	      //if(key==66316990208)
	      //{
	      //  cerr<<"key count "<<endl;
	      //}


	  //#####hash version#####
	      //if either the word or its reverse complement is found in the control don't keep the read
	      if((controlKmers.count(key)>=1) || (controlKmers.count(reversedKey))>=1)
	    {    
	      continue;
	    }else{

		if(qualityReadIn==false)
		  {
		//checking to see if base has poor quality if so reset the string 
		    getline(seqFile, temp);
		    getline(seqFile, qualityScores);
		    quality=int(qualityScores[i])-33;
		    ctr=1;

		    qualityReadIn=true;
		  }else
		  {
		 //assuming sanger style quality scores
		    quality=int(qualityScores[i])-33;

		  }
		 
		 
		 if(quality<10)
		   {
		       key=0;
		       positionCtr=0;
		       continue;
		   }

		 //if uniqueKmer is already present count it and move on
		 if(uniqueKmers.count(key)>1)
		   {
		    
		     uniqueKmers[key]=uniqueKmers[key]+1;
		     continue;
		   }
		 
		
		 //if unique kmer is not already present check to ensure both sides are novel before moving on. 
     
		 secondTime=true; //comment code to check both sides of variant make sure it is novel

	      //if one side of the variant is novel check the other side to see if it is novel
	      if(secondTime==false)
		{
		  firstKey=key;

		  //cout<<"key in count unique kmers are "<<key<<endl;

		  //check to see if read has enough space to check the other side of the variant
		  //if not do not record the variant move on to the next read
		  if((i>=kmerSize))
		    {
		      positionCtr=0;
		      key=0;
		      flag=true;
		      //revKey=0;

		    }else
		    {
		      break;
		    }

		    //need to do this because need to stay on the variant just check the other side of it
		    i=i+1;

		  secondTime=true;
		}else
		{
		
		  //cout<<"key in count unique kmers that go in to hash are "<<key<<endl;

		  //dynamic_bitset<> bitRep(64, firstKey);
		  //cout<<bitRep<<endl;
		  //cout<<line.substr(i, kmerSize)<<"  "<<reversedKey<<" "<<bitRepRev<<endl;
		  //cout<<bit2String(firstKey, kmerSize)<<endl;


		  //dynamic_bitset<> bitRepRev(64, controlCtr1[key/64]);
		  //cout<<bitRepRev<<endl;
		  //cout<<line.substr(i, kmerSize)<<"  "<<reversedKey<<" "<<bitRepRev<<endl;
	
		//cout<<"key in count unique kmers  is "<<key<<" "<<bit2String(key, kmerSize)<<" bit to check is  "<<key%64<<" "<<bitRepRev<<endl;
		//cout<<(controlCtr1[key/64] & (1ull<<(key%64)))<<endl;

		 
		  //uniqueKmers[firstKey]=uniqueKmers[firstKey]+1;
		  uniqueKmers[key]=uniqueKmers[key]+1;
		  uniqueKmers[firstKey]=uniqueKmers[firstKey]+1;
		  //uniqueKmers[reversedKey]=uniqueKmers[reversedKey]+1;
		  
		  //allKmers[firstKey]=allKmers[firstKey]+1;
		  //allKmers[key]=allKmers[key]+1;

		  //secondTime=false;

		}
	     
	      }
      
	    //}
	    
	  
	       
	}

	     	      



    }


  seqFile.close();
}

void readUniqueKmers(dense_hash_map<uint_fast64_t, int, customHash>  &uniqueKmers,  char nextLineFlag, string inputFile, int kmerSize, int ctrCutoff)
{
  ifstream seqFile;
  string line, qualityScores, temp, foo, word, kmer;
  char continueFlag, start;
  int ctr, i, revControlValue, controlValue, quality, kmerCtr;
  long totalCtr;
  bool fastq, flag, secondTime, usedRead, qualityReadIn;
  uint_fast64_t index, reversedKey, firstKey, secondKey;

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  uint_fast64_t key, positionCtr, AbitString, CbitString, GbitString, TbitString; 

  
  vector<string> splitLine;

  //sparse_hash_map<uint_fast64_t, int, customHash> allKmers;
  //dense_hash_map<uint_fast64_t, int, customHash> denseKmers;

  //denseKmers.set_empty_key(-1);
  //denseKmers.resize(3000000000);

  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));
  
  //places into right side of bit string
 GbitString=0*pow(2, (2*kmerSize))+1*pow(2, ((2*kmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

 //place into right side of bitstring
 CbitString=1*pow(2, (2*kmerSize)-2)+0*pow(2, ((2*kmerSize)));


 continueFlag=nextLineFlag;

 if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      exit(EXIT_FAILURE);
    }


 // stringstream ss(inputFile);

 
 seqFile.open(inputFile.c_str(), ifstream::in);
   

  if(!seqFile.is_open())
    {
      cerr<<"could not open file check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }
 

  //building the look up table to reverse nucleotide bit strings
   int bitWord=16;
   std::vector<uint_fast32_t> bitTable(pow(2,bitWord));

   createBitRevTable(bitWord, bitTable);


   //revKey=0;
  key=0;
  positionCtr=0;
  totalCtr=0;
  flag=true;
  ctr=3;
  //arrayIndex=0;
  flag=false;
  firstKey=0;
  secondKey=0;
  while(seqFile.good())  
  {
      totalCtr++;
      ctr++;
      //read in every line of the sequence file
      //getline(seqFile, line);

      //splitLine.clear();
      //reading in kmer counts from a file
      //while(std::getline(seqFile, word, ' ')) {
      //splitLine.push_back(word);
      
      //cout<<"first token is "<<splitLine[0]<<" second token is \t"<<splitLine[1]<<" length of split line is "<<splitLine.size()<<endl;
      
      //}
  
      
      getline(seqFile, line);
   
      
      //object will convert string to a stream so can use things like getline on it
      istringstream iss(line);

      //can parse string this way because words seperated by a space
      iss>>kmer>>word;

      //convert from string into and integer
      kmerCtr=stoi(word);

      //cout<<"kmer is "<<kmer<<" count is  "<<kmerCtr<<endl;

      
      


      if(totalCtr % 10000000==0)
	{
	  
	  cerr<<"lines processed from input file "<<totalCtr<<" size of hash table is "<<uniqueKmers.size()<<endl;
	  	
	  //return;

	  //return(0);
	 
	}
  
      //checking to make sure kmer is counted enough times
      if(kmerCtr<ctrCutoff || kmerCtr>=255)
	{
	  continue;
	}
	
       //ensuring the line is large enough to fit a kmer in it
      if(!(kmer.length()==kmerSize))
       	{
	  cerr<<"kmer processed from input file is wrong size it is  "<<kmer.length()<<" should be "<<kmerSize<<endl;
	  exit(1);
	}
      
      
      //new way of comparision looking to make sure line starts with a valid
      //sequence char
      if (validChar.find(kmer[0]) == std::string::npos) { 
	continue;
      }

  

      positionCtr=0;
      //starting a new sequence
      for(i=kmer.length()-1; i>=0; i--)
	{

	   //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(kmer[i]) == std::string::npos) {
		
	    //revKey=0;
		  key=0;
		  positionCtr=0;
		  //flag=true;
		  continue;
	      }else
		{
		  //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
		  key=key>>2;
	      
		  switch(kmer[i])
		    {
		    case 'A' :
		      key=key|0;
		      break;
		    
		    case 'a' :
		      key=key|0;
		      break;
		   
		    case 'T' :
		      key=key|TbitString;
		      break;
		    
		    case 't' :
		      key=key|TbitString;
		      break;
 		    
		    case 'G' :
		      key=key|GbitString;
		      break;
		    
		    case 'g' :
		      key=key|GbitString;
		      break;

		    case 'C' :
		      key=key|CbitString;
		      break;

		    case 'c' :
		      key=key|CbitString;
		      break;
		      
		    default:
		      cerr<<kmer<<endl;
		      cerr<<"Base is not a nucleotide exiting "<<endl;
		      exit(EXIT_FAILURE);

		    }

		  positionCtr++;		  
		  
		}


	  
	  

	   //not yet reached a full kmer in size
	      if((positionCtr<kmerSize))
		{
		  continue;  
		}
	      
	      /*
	      //16 bit reverse complement 
	      reversedKey = (bitTable[key & 0xffff] << 48) | 
		(bitTable[(key >> 16) & 0xffff] << 32) | 
		(bitTable[(key >> 32) & 0xffff] << 16) |
		(bitTable[(key >> 48) & 0xffff]);  

	      reversedKey=reversedKey>>(64-kmerSize*2);
	      */

	      //cout<<"key is "<<key<<"\t"<<bit2String(key, kmerSize)<<endl;

	      uniqueKmers[key]=kmerCtr;
	 

		 
		  //uniqueKmers[firstKey]=uniqueKmers[firstKey]+1;
	      //		  uniqueKmers[key]=uniqueKmers[key]+1;
	      //  uniqueKmers[firstKey]=uniqueKmers[firstKey]+1;
		  //uniqueKmers[reversedKey]=uniqueKmers[reversedKey]+1;
	

	}
	     
  }
      	  
	 


  seqFile.close();
}

//debugging function for a given cluster ID will print all kmers that point to that cluster
void printSingleCluster( dense_hash_map<uint_fast64_t, long, customHash> &clusterKmers, uint_fast64_t clusterID, int kmerSize)
{

  auto iter=clusterKmers.begin();
  //auto revIter=uniqueKmers.find(reversedKey);

  for(iter; iter!=clusterKmers.end(); iter++)
    {

      if(iter->second==clusterID)
	{
	  cout<<clusterID<<"\t"<<bit2String(iter->first, kmerSize)<<endl;
	}


             //cout<<myHashIterator->first<<"\n\n";

       //cerr<<"right before printing kmer "<<endl;
       //cout<<bit2String(myHashIterator->first, kmerSize)<<"\n\n";
       
       //cerr<<"right before printing reads "<<endl;
       //myHashIterator->second->printReads();
       //continue;
   


    }


}

	     	      
//function to print out the clusters of reads and all unique reads
void printClusters(vector< vector<string> > &clusterBuffer, dense_hash_map<uint_fast64_t, int, customHash> &clusterFiles, vector< string > &uniqueReadsBuffer, ofstream &uniqueOut, vector<std::shared_ptr<ofstream> > &files, int tid)
{
  long j;



  //print out every read and the cluster it belongs to
  for(j=0; j<clusterBuffer.size(); j++)
    {




      //clusterFiles[stoull(clusterBuffer[j][1])]
			  
      //somewhat confusing did this for performance reasons clusterBuffer[j][1] contains the clusterID
      //clusterFiles takes a cluster ID and returns the index into the vector containing the file that cluster belongs to
      //stoull converts from a string to an unsigned long long
      //the vector files contains a set of ofstream pointers to different file output objects
#pragma omp critical(PRINT_ALL_UNIQUE_READS)
      {
	//*(files[clusterFiles[stoull(clusterBuffer[j][1])]])<<clusterBuffer[j][0]<<"\t"<<clusterBuffer[j][1]<<"\n";
      }
		 
    }

  //cerr<<"size of vector is "<<uniqueReadsBuffer.size()<<"\n";
  //cerr<<uniqueReadsBuffer[uniqueReadsBuffer.size()-1]<<"\n";
  //uniqueOut<<uniqueReadsBuffer[1]<<"\n";

  



  		      
  for(j=0; j<(uniqueReadsBuffer.size()-1); j+=2)
    {


#pragma omp critical(PRINT_CLUSTERS)
      {

	uniqueOut<<"@"<<j<<"\n"; //fastq read header
	uniqueOut<<uniqueReadsBuffer[j]<<"\n"; //sequence
	uniqueOut<<"+"<<"\n"; //spacer line
	uniqueOut<<uniqueReadsBuffer[j+1]<<"\n"; //quality scores
      
      }
    
    }
  		      
		      //ensure all output is flushed to the disk only need this if debugging
		      //uniqueOut.flush();
}





void assignClusters(node * workNodePtr, //a pointer pointing to a chunk of work 
		    dense_hash_map<uint_fast64_t, long, customHash> &clusterKmers, //hash table key is kmer value is the cluster that kmer belongs to
		    dense_hash_map<uint_fast64_t, int, customHash> &clusterFiles,  //hash table key is a cluster id value is the index into the vector of output file names
		    vector<std::shared_ptr<ofstream> > &files, //set of opened file handlers that clusters will be written to
		    dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers, //hash table key is a kmer that is unqiue to the exp read library value is the number of times kmer occurs in experiment read library
		    ofstream &uniqueOut, //ofstream object file handler contains location of output file for all unique reads
		    int kmerSize, //kmer size
		    int numFiles, //number of files to print clusters to
		    uint_fast64_t &clusterCtr, //count of how many clusters there currently are will be shared among threads
		    int tid //thread id number passing this for debugging purposes
)
{


  string line, header, tempLine;
  char continueFlag, start;
  long ctr, i, revControlValue, controlValue, kmerCounter, j, NCtr;
  long totalCtr, startFirstKmer, startSecondKmer, lineCtr;
  bool fastq, flag, secondTime, usedRead, presentCluster, foundCluster, isUnique;
  uint_fast64_t index, reversedKey, firstKey, secondKey, tempKey, posKey, negKey, first, middle, end, currentUniqueKey;
  uint_fast64_t lastKey, clusterIndex, randomFileIndex, tempValue;
  int posMaxInter, posSecondBestInter, negMaxInter, negSecondBestInter;
 
  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  uint_fast64_t key, positionCtr, AbitString, CbitString, GbitString, TbitString; 


  //creating a vector of vectors to hold the set of reads that belong to a cluster
  //first element of vector is the read second element is the cluster the read belongs to
  vector< vector<string> > clusterBuffer;
  clusterBuffer.reserve(85000);


  //creating a vector to act as a buffer to hold reads that contain a unique word
  vector< string > uniqueReadsBuffer;
  clusterBuffer.reserve(100000);

 
  //building the look up table to reverse nucleotide bit strings
   int bitWord=16;
   std::vector<uint_fast32_t> bitTable(pow(2,bitWord));

   createBitRevTable(bitWord, bitTable);




  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));
  
  //places into right side of bit string
 GbitString=0*pow(2, (2*kmerSize))+1*pow(2, ((2*kmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

 //place into right side of bitstring
 CbitString=1*pow(2, (2*kmerSize)-2)+0*pow(2, ((2*kmerSize)));





   //#pragma omp parallel num_threads(thread_count) default(shared)	\
   //private(tid, i, workNode, sum)
  
  
 continueFlag='0'; //hard coded value a relic of older code that had functionality written into it such that 2 different lines could be concatenated together. 


   //revKey=0;
  key=0;
  positionCtr=0;
  totalCtr=0;
  ctr=3;
  firstKey=0;
  fastq=true;



  for(j=0; j<workNodePtr->chunk.size(); j++)
    {



      //if sequence does not continue on the next line reset all values to zero 
      //since sequence on next line is unrelated to current sequence
      if(continueFlag=='0')
	{
	  //revKey=0;
	  key=0;
	  positionCtr=0;
	  flag=true;
	}
      
      secondTime=false;
      usedRead=false;
      firstKey=0; 
      secondKey=0;
      kmerCounter=0;
      NCtr=0;
    


      //cerr<<"line is "<<line<<endl;

//starting a new sequence
      //confuseing doing this for performance reasons
      //workNodePtr points to a vector of sequences
      //I am working with pointers to avoid unnecassary copying of sequence
      //chunk[j] is a string chunk[j][i] is the ith character in the string
      for(i=workNodePtr->chunk[j].length()-1; i>=0; i--)
	{

	  //secondKey=reversedKey;
	  //startSecondKmer=i+1;

	   //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(workNodePtr->chunk[j][i]) == std::string::npos) {
		
	    //revKey=0;
	    key=0;
	    positionCtr=0;
	    //NCtr+=1; //Counts number of Ns in line
		  //flag=true;
	    continue;
	  }else
	    {
	      //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
	      key=key>>2;
	      
	      switch(workNodePtr->chunk[j][i])
		{
		case 'A' :
		  key=key|0;
		  break;
		  
		case 'a' :
		  key=key|0;
		  break;
		  
		case 'T' :
		  key=key|TbitString;
		  break;
		      
		case 't' :
		  key=key|TbitString;
		  break;
 		  
		case 'G' :
		  key=key|GbitString;
		  break;
		  
		case 'g' :
		  key=key|GbitString;
		  break;

		case 'C' :
		  key=key|CbitString;
		  break;

		case 'c' :
		  key=key|CbitString;
		  break;
		      
		default:
		  
		  cerr<<"line is "<<workNodePtr->chunk[j]<<endl;
		  cerr<<"Base is not a nucleotide exiting "<<endl;
		  exit(EXIT_FAILURE);

		}

	      positionCtr++;		  
		  
	    }

	   //not yet reached a full kmer in size
	      if((positionCtr<kmerSize))
		{
		  continue;  
		}
	      

	      //16 bit reverse complement 
	      reversedKey = (bitTable[key & 0xffff] << 48) | 
		(bitTable[(key >> 16) & 0xffff] << 32) | 
		(bitTable[(key >> 32) & 0xffff] << 16) |
		(bitTable[(key >> 48) & 0xffff]);  

	      reversedKey=reversedKey>>(64-kmerSize*2);


	      isUnique=false;

	      //	      continue;

	      
	      //lookup key if key not found will return iterator that points to end of hash table
	      auto iter=uniqueKmers.find(key);
	      auto revIter=uniqueKmers.find(reversedKey);

	      	
	      if(iter!=uniqueKmers.end() || revIter!=uniqueKmers.end() )
	      {


		//if count of unique kmers is less than the cutoff
		if(iter->second < 5 && revIter->second < 5)
		  {
		    continue;
		  }


		  //if any unique kmer is larger than the cutoff do not use that read
		  if(iter->second > 254 || revIter->second > 254)
		  {
		    break;
		  }

		  
		  
                 
		 
		  
		  //if clusterBuffer reaches the maximum size then output each read sequence to the correct 
		  //output file the cluster the read belongs to is stored in
		  if(clusterBuffer.size()>10000)
		    {

		     
		      if(uniqueReadsBuffer.size() > 0 && clusterBuffer.size() > 0)
			{
			  printClusters(clusterBuffer, clusterFiles, uniqueReadsBuffer, uniqueOut, files, tid);
		      	}

	      
		      clusterBuffer.clear();
		      uniqueReadsBuffer.clear();
		    }


		  isUnique=true;
		  currentUniqueKey=key;
		 
		  
		  //dynamic_bitset<> bitRep(64, key);
		  //cout<<bitRep<<endl;
		  //cout<<i<<" "<<workNodePtr->chunk[j].substr(i, kmerSize)<<"  before reversal "<<key<<" "<<bitRep<<endl;
		  //cout<<workNodePtr->chunk[j]<<endl;
		  //cout<<bit2String(key, kmerSize)<<endl;
		  


		  //cout<<"key in get REads is "<<uniqueKmers[key]<<" "<<uniqueKmers[reversedKey]<<endl;
		  
		  //checking to make sure unique kmer is represented on the other strand  
		  //if((uniqueKmers[key]>=1) && (uniqueKmers[reversedKey]>=1))
		  //		  if((uniqueKmers[key]>=1) || (uniqueKmers[reversedKey]>=1))  
		  //if((uniqueKmers[key]>=1) && (uniqueKmers[reversedKey]>=1))
		  //{  
		      
		  if(kmerCounter==0){

		    firstKey=key;
		   
		    //maybe keep this in
		    //if(i==0)
		    //{
		    //	break;
		    //}
    
		    
		    if(fastq)
		      {
			 
			if(kmerCounter==0)
			  {
			  	      
			    
			    uniqueReadsBuffer.push_back(workNodePtr->chunk[j]); //store the read sequence
			    uniqueReadsBuffer.push_back(workNodePtr->qualityScores[j]); //store the quality scores

			    //cout<<tempLine<<endl;
			      //cout<<">?>:BCC<ADCCBA;AABA7@ABBAA@CCBA>ABCAAAAACCAA@><B<A??BCA?>>?AAA=>?A:;5:6AB=6/5;;37:6@=A><;>;@:>:B>7###"<<endl;

			    ctr=1;
			    
	
			  }

		      }else
		      {

			cerr<<"Need to implement output unique reads for fasta format in getReads function "<<"\n";
			exit(EXIT_FAILURE);

			if(kmerCounter==1)
			  {
			    cout<<">"<<key<<endl;
			    cout<<workNodePtr->chunk[j]<<endl;
			  }
		      }
		    
		      
		    
		
			//cout<<"first word pos is "<<key<<"\t"<<bit2String(key, kmerSize)<<endl;
			//cout<<"first word neg is "<<reversedKey<<"\t"<<bit2String(reversedKey, kmerSize)<<endl;
		  
			//cerr<<"first key ihigher up is "<<firstKey<<endl;

		  }

		  
		            
		      
		      //dynamic_bitset<> bitRep(64, key);
		      //cout<<bitRep<<endl;
		      //cout<<i<<" "<<workNodePtr->chunk[j].substr(i, kmerSize)<<"  before reversal "<<key<<" "<<bitRep<<endl;
		      //cout<<workNodePtr->chunk[j]<<endl;
		      //cout<<bit2String(key, kmerSize)<<endl;
		      //exit(1);
		          
		      
		  kmerCounter+=1;	
		      
	      }




	      //this means already have found at least one unique kmer that fits the criteria earlier in the read
	      //add one to each kmer
	      if(kmerCounter>=1)
		{


		  //continue;

		  if(isUnique)
		    {
		  //check to see if read can be assigned to an existing cluster
		      if(clusterKmers.count(key)>0)
			{


			  clusterIndex=clusterKmers[key];

			  
			  //dynamic_bitset<> bitRep(64, key);
			  //cout<<bitRep<<endl;
			  //cout<<workNodePtr->chunk[j]<<endl;
			  //cout<<i<<" "<<workNodePtr->chunk[j].substr(i, kmerSize)<<"  before reversal "<<key<<" cluster is "<<clusterIndex<<endl;
			  //cout<<"position ctr is "<<positionCtr<<endl;
			  //cout<<bit2String(key, kmerSize)<<endl;

			  
			  /*
			  if(workNodePtr->chunk[j].compare("TCTAACCTACTGGAGGATAAGGGGTGAGCCCAAGAGCCTCAAGGCTCCCATCAACAGCCAGTCCTGTGAGTGCGGCCATCTTGGACCTGCCAGCTCAGTAA")==0)
			    {
			      cerr<<i<<" "<<workNodePtr->chunk[j].substr(i, kmerSize)<<"  before reversal "<<key<<" cluster is "<<clusterIndex<<endl;
			      cerr<<"position ctr is "<<positionCtr<<endl;
			      cerr<<bit2String(key, kmerSize)<<endl;

			      exit(0);
			    }
			  */

		
			  //adding a new line to the existing buffer that holds cluster reads
			  clusterBuffer.push_back(vector<string>{workNodePtr->chunk[j], to_string(clusterIndex)} );
		     
		      //if first word not already in cluster then add it
			  //do not have to worry about last word since last word is  the current word 
			  if(kmerCounter>1)
			    {
			      if(clusterKmers.count(firstKey)==0)
				{
                                 #pragma omp critical(ADD_KEY_EXISTING_CLUSTER)
				  {
				  
				    clusterKmers[firstKey]=clusterIndex;
				
				  }
				}
			    }

			  break;
			}else if(clusterKmers.count(reversedKey)>0)
			{
			      //if read matches existing cluster on reverse strand reverse complement and add to cluster
			  clusterIndex=clusterKmers[reversedKey];
			  revComplement(workNodePtr->chunk[j]);

			  //adding a new line to the existing cluster
			  clusterBuffer.push_back(vector<string>{workNodePtr->chunk[j], to_string(clusterIndex)} );
		      

			  

			  /*
			  
			  if(workNodePtr->chunk[j].compare("TCTAACCTACTGGAGGATAAGGGGTGAGCCCAAGAGCCTCAAGGCTCCCATCAACAGCCAGTCCTGTGAGTGCGGCCATCTTGGACCTGCCAGCTCAGTAA")==0)
			    {
			      //cout<<i<<" "<<workNodePtr->chunk[j].substr(i, kmerSize)<<"  before reversal "<<key<<" cluster is "<<clusterIndex<<endl;
			      cout<<"position ctr in reversed if condition is "<<positionCtr<<endl;
			      cout<<bit2String(reversedKey, kmerSize)<<endl;

			      
			      cout<<"cluster it is joining is "<<endl;
			      printSingleCluster(clusterKmers, clusterIndex, kmerSize);


			      exit(0);
			    }
			  */







			  //reverse complement the first key in the new read and add it to the cluster assuming it does not already exist
			  reversedKey = (bitTable[firstKey & 0xffff] << 48) | 
			    (bitTable[(firstKey >> 16) & 0xffff] << 32) | 
			    (bitTable[(firstKey >> 32) & 0xffff] << 16) |
			    (bitTable[(firstKey >> 48) & 0xffff]);  

			  reversedKey=reversedKey>>(64-kmerSize*2);
			  
			  if(kmerCounter>1)
			    {
			      if(clusterKmers.count(reversedKey)==0)
				{
				 #pragma omp critical(ADD_REV_KEY_EXISTING_CLUSTER)
				  {
				 
				    clusterKmers[reversedKey]=clusterIndex;

				  }
				}
			    }

			  break;  
			}
		    }
		

		  //if have obtained the set of all possible unique kmers or reached the end of the read or reached the end of unique kmers
		  //then generate a new cluster if no assignment can be made
		  ///		  if(i==0)
		  //if(kmerCounter==kmerSize || i==0 || i==workNodePtr->chunk[j].length()-1 || !(isUnique))
		  if(i==0)  
		  {

		      //cerr<<"inside this if statement "<<endl;
		      //exit(0);



                     #pragma omp critical(CREATE_CLUSTER)
		      {
		

			/*
 			if(workNodePtr->chunk[j].compare("ATCAGAGTTAACTTTGCAGTGAGAGCGGCCTTGCTGCGGCCAAAGAACATGGAAAAGCATGANTGGGGTGATGTGCCTTAAAGCGTCAGACACTTGGGCCT")==0)
			  {
			    cerr<<"first key in make new cluster is "<<bit2String(firstKey, kmerSize)<<endl;
			    cerr<<"second key is "<<bit2String(currentUniqueKey, kmerSize)<<endl;
			  

			    printClusters(clusterBuffer, clusterFiles, uniqueReadsBuffer, uniqueOut, files, tid);

			    clusterBuffer.clear();
			    uniqueReadsBuffer.clear();

			    //exit(0);
			  }
			*/
		      


	
			clusterCtr++;

			clusterKmers[firstKey]=clusterCtr;
			clusterKmers[currentUniqueKey]=clusterCtr;
		
			randomFileIndex=rand()%numFiles;
		     
		      //pick a random file to output cluster to
			clusterFiles[clusterCtr]=randomFileIndex;

			//assign add a new row to the cluster matrix holding the reads
			//only possible using C++11
			clusterBuffer.push_back(vector<string>{workNodePtr->chunk[j], to_string(clusterCtr)} );
                
		      }

		      break;
					

		      
		    
 
		      
		    }
		    

		}

 		   
	 
	      
     
		       
	}

    }





  if(uniqueReadsBuffer.size() > 0 && clusterBuffer.size() > 0)
    {
  //flush whatever is left in the buffers
      printClusters(clusterBuffer, clusterFiles, uniqueReadsBuffer, uniqueOut, files, tid);
    }
		      

  //free the memory allocated to the work node
  delete workNodePtr;
		      

  //sparse_hash_map<uint_fast64_t, long, customHash> clusterKmers; //hash table key is kmer value is the cluster that kmer belongs to

  /*
  sparse_hash_map<uint_fast64_t, long, customHash>::iterator myHashIteratorTemp;
  
   //unordered_map<uint_fast64_t, char>::iterator myHashIterator;
  
   for(myHashIteratorTemp=clusterKmers.begin(); myHashIteratorTemp!=clusterKmers.end(); myHashIteratorTemp++)
     {
       
	   cout<<myHashIteratorTemp->first<<"\t"<<myHashIteratorTemp->second<<endl;
       
     }
  */


  




}

//dense_hash_map<uint_fast64_t, int, customHash>

//Given a hash table containing unique kmers and their counts get reads that have unique kmers above a certain threshold

//void getReads(sparse_hash_map<uint_fast64_t, int, customHash> &uniqueKmers, sparse_hash_map<uint_fast64_t, ReadCluster *, customHash> &readClusters, char nextLineFlag, string inputFile, int kmerSize, dense_hash_map<uint_fast64_t, uint_fast64_t, customHash> &masterKey )

vector<string> getReads(dense_hash_map<uint_fast64_t, int, customHash> &uniqueKmers, int numFiles, char nextLineFlag, string inputFile, int kmerSize)
{
  ifstream seqFile;
  string line, header, tempLine, qualityScore;
  char continueFlag, start;
  long ctr, i, revControlValue, controlValue, kmerCounter, j, NCtr;
  long totalCtr, startFirstKmer, startSecondKmer, lineCtr;
  bool fastq, flag, secondTime, usedRead, presentCluster, foundCluster, isUnique;
  uint_fast64_t index, reversedKey, firstKey, secondKey, tempKey, posKey, negKey, first, middle, end, currentUniqueKey;
  uint_fast64_t lastKey, clusterIndex, clusterCtr, randomFileIndex, tempValue;
  int posMaxInter, posSecondBestInter, negMaxInter, negSecondBestInter;
 
  int chunkSize=40000, loadedFileFlag, thread_count, tid, chunkCtr, maxChunksWork;
  node *currentNodePtr, *workNodePtr;

  list<node*> workList;
 

  vector<string> fileNames; //vector containing files names clustes are placed into

  //initilize random number generator
  srand (time(NULL));

  dense_hash_map<uint_fast64_t, int, customHash> refPositions; //hash table contains the positions of every unique word in the reference read
  refPositions.set_empty_key(-1);

 

  dense_hash_map<uint_fast64_t, int, customHash> clusterFiles; //hash table key is a cluster id value is the index into the vector of output file names
  clusterFiles.set_empty_key(-1);
  clusterFiles.resize(10000000);

//hash table to hold the kmers assocated with each cluster
  dense_hash_map<uint_fast64_t, long, customHash> clusterKmers(10000000); //hash table key is kmer value is the cluster that kmer belongs to
  clusterKmers.resize(100000000);
//dense_hash_map<uint_fast64_t, long, customHash> clusterKmers; //hash table key is kmer value is the cluster that kmer belongs to
  clusterKmers.set_empty_key(-1);



  
  //sparse_hash_map<uint_fast64_t, vector<string>, customHash> clusterBuffer(100000); //hash table to temporarily store a set of sequences and the cluster they belong to
 

  //creating a vector of vectors to hold the set of reads that belong to a cluster
  //first element of vector is the read second element is the cluster the read belongs to
  vector< vector<string> > clusterBuffer;
  clusterBuffer.reserve(85000);


  //creating a vector to act as a buffer to hold reads that contain a unique word
  vector< string > uniqueReadsBuffer;
  clusterBuffer.reserve(100000);



//will periodically be flushed to a set of output files key is cluster id value is sequence
 //clusterBuffer.resize(1000000);

  vector<uint_fast64_t> posWords; //words generate from sequence variant on pos strand
  vector<uint_fast64_t> negWords; //words generated from sequence variant on neg strand
  vector<uint_fast64_t> intersection; //vector containing the intersection of words
  vector<uint_fast64_t>::iterator it; //iterator points to the end of the set intersection comparing the words of each read to the words of each read cluster 
  vector<uint_fast64_t> posMatchingWords; //list of words from a read that match words on the master read
  vector<uint_fast64_t> negMatchingWords;


  //vector of pointers to ofstream objects 
  vector<std::shared_ptr<ofstream> > files;

  //creating and opening a set of output files to place clusters into
  for(i=0; i<numFiles; i++)
    {
      line="/home/massa/Temp/file"+to_string(i)+".dat";
      
      //line="/data1/HEM0013-131-LYMPH.bam.aspera-env et al/Temp/"+to_string(i)+".dat";
      fileNames.push_back(line);
      //line="/data1/HEM0013-131-LYMPH.bam.aspera-env et al/SingleCore/"+to_string(i)+".dat";
      files.push_back(make_shared<ofstream>(line, std::ios::out ));
    }
  

  //ofstream uniqueOut("/data7/SeqDiffResults/Results/uniqueReads.fastq");

  //ofstream uniqueOut("/data1/HEM0013-131-LYMPH.bam.aspera-env et al/uniqueReads.fastq");
    ofstream uniqueOut("/data/Sequencing/kmerAnalysis/uniqueReads.fastq");



  //*(files[0])<<"testing"<<endl;
  //*(files[0])<<"testing123"<<endl;

 
//files[0]->close();


  //exit(0);

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";

  uint_fast64_t key, positionCtr, AbitString, CbitString, GbitString, TbitString; 

  sparse_hash_map<uint_fast64_t, ReadCluster *, customHash> ::iterator myHashIterator; //iterater for read cluster
  
  dense_hash_map<uint_fast64_t, int, customHash> readPositions;
  readPositions.set_empty_key(-1);


  //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*kmerSize))+1*pow(2, 64-((2*kmerSize)-1));
  
  //places into right side of bit string
 GbitString=0*pow(2, (2*kmerSize))+1*pow(2, ((2*kmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*kmerSize))+0*pow(2, 64-((2*kmerSize)-1));

 //place into right side of bitstring
 CbitString=1*pow(2, (2*kmerSize)-2)+0*pow(2, ((2*kmerSize)));


 continueFlag=nextLineFlag;

 if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      exit(EXIT_FAILURE);
    }


  seqFile.open(inputFile.c_str());
   

  if(!seqFile.is_open())
    {
      cerr<<"could not open file check to see if exists"<<endl;
      exit(EXIT_FAILURE);
    }
  
  
   getline(seqFile, header);
   start=toupper(header[0]);

   //header=line;
   
   fastq=false;
   if(start!='>')
     {
       fastq=true; 
     
       getline(seqFile, line);
       getline(seqFile, tempLine);
       getline(seqFile, qualityScore);

    
     }
     
   

   //################################33

   key=0;
   positionCtr=0;
   totalCtr=0;
   ctr=3;
   firstKey=0;
   clusterCtr=0;
   

  loadedFileFlag=0;
  thread_count=20;
  lineCtr=1;
  chunkCtr=0;
  maxChunksWork=10;

  //adding the first node
  currentNodePtr=new node;

#pragma omp parallel num_threads(thread_count) default(shared)\
  private(tid, i, workNodePtr)
  {

    tid=omp_get_thread_num();
  
    //#pragma omp critical
    //{
    //cerr<<"tid is "<<tid<<endl;
    //}

    while(loadedFileFlag==0)
      {

	//only access input file if you are the master thread
	#pragma omp master
	{
      //if list is empty get some work
	  if(workList.empty())
	    {
	      chunkCtr=0;

	      //while there is still data and the file and have not filled up the work list 
	      //add chunks of work in this case add no more than maxChunksWork or work 
	      while(seqFile.good() && chunkCtr < maxChunksWork)
		{
		 
		  if((lineCtr % chunkSize)==0)
		    {
		      lineCtr=1;
		      chunkCtr++;

		      //if(lineCtr>0)
		      //{
			  
                         #pragma omp critical(ADD_WORK)
			  {
			  
			    //cerr<<"pushing back  a chunk of work  chunkCtr is "<<chunkCtr<<endl;
			    workList.push_back(currentNodePtr);
			    //cerr<<"number of chunks of work is "<<workList.size()<<" tid is "<<tid<<endl;
			  }
			  
			  //}

		      currentNodePtr=new node;
		    }



		  
		  totalCtr++;
		  //ctr++;

		  if(fastq)
		    {
		      
		      //read in an entire fastq record

		      getline(seqFile, header);
		      getline(seqFile, line);
		      getline(seqFile, tempLine);
		      getline(seqFile, qualityScore);

		    }
     		  
       
      //cout<<totalCtr<<endl;
      if(totalCtr % 10000000==0)
	{	  
	  //	  cerr<<"lines processed from input file "<<totalCtr<< " number of clusters is "<<readClusters.size()<<endl;
	
         #pragma omp critical(PRINT_UPDATE)
	  {

	    cerr<<"lines processed from input file "<<totalCtr<<"   number of cluster kmers is "<<clusterKmers.size()<<"\n";
	    cerr<<"number of clusters is "<<clusterFiles.size()<<"\n";
	  
	  }
//cerr<<"number of chunks of work is "<<workList.size()<<endl;
  //cerr<<"number of reads assigned to clusters are "<<clusterBuffer.size()<<"\n";
	  //cerr<<"number of clusters is "<<clusterFiles.size()<<"\n";

	  // cerr<<"number of unique Kmers is "<<uniqueKmers.size()<<"\n";
	  //cerr<<"bucket count of unique kmers is "<<uniqueKmers.bucket_count()<<"\n";

	  /*
	  for(int k=0; k<10000000000; k++)
	    {
	      int foo=1;
	    }
	  */
	 
	  //return;

	  //return(0);
	 
	}
   
     
      //continue;
      /*	
      if(fastq)
	{
      
	  if(ctr==1)
	    {
	      qualityScore=line;
	    }

	  //if(ctr==2)
	  // {
	  //  header=line;
	  //}

	  //very very odd if ctr is set to 0 
	  //code slows down dramatically I have no 
	  //idea why
	  if(ctr<3)
	    {	      
	      continue;
	    }else
	    {
	     

	      ctr=-1;
	    }
	  
	  
	}
      */
      //continue;
	
       //ensuring the line is large enough to fit a kmer in it
      if(!(line.length()>= kmerSize))
	{
	  continue;
	}
      
      
      //new way of comparision looking to make sure line starts and ends with a valid
      //sequence char
      if ((validChar.find(line[(line.size()-1)]) == std::string::npos) || (validChar.find(line[0]) == std::string::npos) ) { 
	continue;
      }

   
     



       //if sequence does not continue on the next line reset all values to zero 
      //since sequence on next line is unrelated to current sequence
      // if(continueFlag=='0')
      //{
	  //revKey=0;
      //  key=0;
      //  positionCtr=0;
      //  flag=true;
	  
      //}





      //		  getline(seqFile, line);

      ///REMEMBER TO DELETE ALLOCATED WORK NODE WHEN YOU POP IT OFF////
		  lineCtr++;
		  //cerr<<"line ctr is  "<<lineCtr<<" tid is "<<tid<<"\t"<<endl;
		  currentNodePtr->chunk.push_back(line);
		  currentNodePtr->qualityScores.push_back(qualityScore);
		  //cerr<<"right after push back operation "<<"tid is "<<tid<<endl;

		  //cerr<<"line is "<<line.length()<<endl;
		  //cerr<<"quality is "<<qualityScore.length()<<endl;
	      
		  //exit(0);

		}
	    }
	    
	    
	  if(!seqFile.good())
	    {
	
	      //cerr<<"last Node size is "<<currentNodePtr->chunk.size()<<"\n";
	  //add the very last bit of file
            #pragma omp critical(ADD_LAST_WORK)
	      {
		workList.push_back(currentNodePtr);
	      }
	      
	      loadedFileFlag=1;
	
	    }


	}

	workNodePtr=NULL;
     	
           #pragma omp critical(GET_WORK) 
	    {
	      if(!workList.empty())
		{
		  
		  
		  //cerr<<"current size of workList is "<<workList.size()<<endl;
		  workNodePtr=workList.back(); //get a chunk of work from the list
		
		  workList.pop_back(); //remove the chunk of work that was just copied so another thread doesn't do the same chunk of work
		  //working with pointers to nodes 
	
		}
	    } 
		  /*
		    cout<<"thread is ********** "<<tid<<"\n";
		    for(i=0; i<workNodePtr->chunk.size(); i++)
		      {	  
		    
			cout<<workNodePtr->chunk[i]<<"\n";
				
		      } 
		    
		  */

	    //have some work 
	    if(workNodePtr)
	      {
	
		assignClusters(workNodePtr, clusterKmers, clusterFiles, files, uniqueKmers, uniqueOut, kmerSize, numFiles, clusterCtr, tid);
	      
	      }
    
	
      }

  }

  //flush whatever is left in the buffers
  //printClusters(clusterBuffer, clusterFiles, uniqueReadsBuffer, uniqueOut, files);		      

  seqFile.close();


  //sparse_hash_map<uint_fast64_t, long, customHash> clusterKmers; //hash table key is kmer value is the cluster that kmer belongs to

  /*
  sparse_hash_map<uint_fast64_t, long, customHash>::iterator myHashIteratorTemp;
  
   //unordered_map<uint_fast64_t, char>::iterator myHashIterator;
  
   for(myHashIteratorTemp=clusterKmers.begin(); myHashIteratorTemp!=clusterKmers.end(); myHashIteratorTemp++)
     {
       
	   cout<<myHashIteratorTemp->first<<"\t"<<myHashIteratorTemp->second<<endl;
       
     }
  */


  return(fileNames);

}



//function to read in Clusters do some filtering and pass the filtered Clusters to be assembled
void readInClusters(vector<string> &fileNames, int cutoffClusterSize, int clusterKmerSize)
{
  
  ifstream clusterFile;
  string line, read, clusterIdentifier;
  uint_fast64_t clusterID, i;
  // ReadCluster cluster;


  dense_hash_map<uint_fast64_t, vector<string> *, customHash> clusters; //hash table key is cluster ID value is a pointer to the  set of reads that belong to that cluster
  clusters.set_empty_key(-1);

  for(i=0; i<fileNames.size(); i++)
    {
  
      clusterFile.open(fileNames[0].c_str(), ifstream::in);
      if(!clusterFile.is_open())
	{
	  cerr<<"could not open file "<<fileNames[0]<< " check to see if exists"<<endl;
	  exit(EXIT_FAILURE);
	}

      while(clusterFile.good())  
	{
      
	  getline(clusterFile, line);
        
      //object will convert string to a stream so can use things like getline on it
	  istringstream iss(line);
      
      //can parse string this way because words seperated by a space
	  iss>>read>>clusterIdentifier;

      //convert from string into and integer
	  clusterID=stoi(clusterIdentifier);

      //check to see if cluster already exists in hash table
      //if cluster exists then add read to existing cluster
      //otherwise add an entry into hash table
	  if(clusters.count(clusterID)>0)
	    {
	      clusters[clusterID]->push_back(read);
	    }else
	    {
	      clusters[clusterID]=new vector<string>; //##########NEED to Delete Memory MEMORY LEEK REMEMBER TO HANDLE###################
	      clusters[clusterID]->push_back(read);//###############IMPORTANT LOREN SEE ABOVE##################
	    }


	}
    
      clusterFile.close();
    


      auto iter=clusters.begin();
      
      //int temp=0;

      //go through the set of clusters removing duplicates doing some filtering based on number of remaining reads 
      //then assemble reads using seqAn library functions
      for(iter; iter!=clusters.end(); iter++)
	{

	  //cout<<"size of cluster "<<iter->first<<" is "<<iter->second->size()<<endl;

	  /*
	  auto vectorIter=iter->second->begin();
	  for(vectorIter; vectorIter!=iter->second->end(); vectorIter++)
	    {
	      cout<<*vectorIter<<endl;
	    }
	  */

	  //deduplicating reads
	  sort( iter->second->begin(), iter->second->end() );
	  //iter->second->erase( unique( iter->second->begin(), iter->second->end() ), iter->second->end() );
	  auto last=unique(iter->second->begin(), iter->second->end());
	  iter->second->erase(last, iter->second->end());


	  //number of reads in cluster cutoff
	  if(iter->second->size()> cutoffClusterSize)
	    {
	     

	      ReadCluster cluster(iter->second->size()); //input paramater is number of reads in the cluster

	    
	      //adding the sequences to the cluster
	      cluster.addSequences(*(iter->second));
	      cluster.setKmerSize(clusterKmerSize);
	      cluster.getKmers();
	      
	      std::pair<uint_fast64_t, long> pair=cluster.getPair(0);
	      
	      //cout<<cluster.getNumReads()<<"\t"<<pair.second<<endl;

	      //if(pair.second==cluster.getNumReads())
	      //{
		  //cerr<<pair.second<<"\t"<<pair.first<<"\t number of Reads is "<<cluster.getNumReads()<<endl;
		  
		  cluster.printReads();
		  cluster.setStartPositions(pair.first);
		  cluster.printStartPositions();
	
		  cluster.mergeReads(clusterKmerSize, 3);
		  //cout<<bit2String(pair.first, clusterKmerSize)<<endl;
		  
		  //temp++;

		  //if(temp>5)
		  //{
		  //  exit(0);
		  //}
		  //cout<<"################\n\n\n\n";
		  //}


	      //cluster.printKmerPositions();

	      //cluster.printReads();

	      //cout<<"Cluster now is #######################"<<"\n\n\n\n\n\n";

	      //cout<<"Cluster ID is "<<iter->first<<" Number of reads are "<<iter->second->size()<<endl;

	    }

	  //cout<<"###############################################################\n\n";

	  /*
	  cout<<"size now is "<<iter->second->size()<<endl;
	  vectorIter=iter->second->begin();
	  for(vectorIter; vectorIter!=iter->second->end(); vectorIter++)
	    {
	      cout<<*vectorIter<<endl;
	    }
	    
	  cout<<"###################\n\n";

	  */

             //cout<<myHashIterator->first<<"\n\n";

       //cerr<<"right before printing kmer "<<endl;
       //cout<<bit2String(myHashIterator->first, kmerSize)<<"\n\n";
       
       //cerr<<"right before printing reads "<<endl;
       //myHashIterator->second->printReads();
       //continue;
   


	}



    }
 

  

  

}
