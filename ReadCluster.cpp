#include "ReadCluster.h"

using google::dense_hash_map;
using google::sparse_hash_map;
using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::ofstream;
using std::ios;
using std::cerr;
using std::pair;

ReadCluster::ReadCluster(){

  contig="0";
  //readSeq.push_back("test2");

}

int ReadCluster::getNumReads()
{

  return readSeqs.size();
}


vector<string> ReadCluster::getSeqs(){
  
  return readSeqs;

}

string ReadCluster::getContig()
{
  return contig;
}


void ReadCluster::addSeq(string &read)
{
  readSeqs.push_back(read);
}


void ReadCluster::addSequences(vector<string>& allReads)
{
  readSeqs=allReads;
  
  }

void ReadCluster::printReads()
{

  int i;

  //cout<<readSeqs[i]<<endl;
  for(i=0; i<readSeqs.size(); i++)
    {
      //readSeq[i].insert(startPositions[i], "  ");
      //cout<<readSeq[i]<<"\t"<<startPositions[i]<<endl;

      cout<<readSeqs[i]<<endl;

    }


  /*
  string::iterator it;

  for(it=readSeq.begin(); it<readSeq.end(); it++)
    {
      cout<<it*<<"\n";
    }
  */
}

//kmer size and a cutoff representing the minimum number of times the max nucleotide must occur to be counted
string ReadCluster::mergeReads(int kmerSize, int cutoffMinNuc)
{

  long int numCol;
  int numReads=readSeqs.size();
  int i, j, z, maxCtr, indexMax, posN, readSize;
  vector<int> nuc(4, 0); //A T G C
  string combinedNuc;
  string allNucs="ATGC";
  bool flag=false;
  double Nctr, nucCtr, badColCtr;
 


  //finding the maximum read size
  readSize=0;
  for(i=0; i<numReads; i++)
    {
      if(readSeqs[i].length() > readSize)
	{
	  readSize=readSeqs[i].length();
	}
      
    }

  numCol=(readSize*2);
  

			   
  //matrix to hold reads
  vector<vector<char>>matrix(numReads, vector<char>(numCol, 'N'));

  //cout<<"size of matrix is "<<matrix.size()<<" number of cols is "<<matrix[0].size()<<endl;


  


  /*
  //creating a matrix  sequences are lined up according to matching kmers if a particular sequence has an element in the matrix it will be filled in
  //otherwise it will be left as an N
  for(i=0; i<numReads; i++)
    {
      for(j=0; j<readSeq[i].length(); j++)
	{
	  matrix[i][((numCol/2)-startPositions[i])+j]=readSeq[i][j];
	}
    }
  */

  /*
//creating a matrix  sequences are lined up according to matching kmers if a particular sequence has an element in the matrix it will be filled in
  //otherwise it will be left as an N
  for(i=0; i<numReads; i++)
    {
      if(i>0)
	{
	  cout<<"printing start positions for reads "<<i<<endl;
	  cout<<"number of reads is "<<numReads<<endl;
	  for(int z=0; z<startPositions.size(); z++)
	    {
	      cout<<startPositions[z]<<endl;
	    }
	}
    

      cout<<"finished printing start positions "<<endl;
    }

  cout<<"test"<<endl;
  */

  //file to put debuggin info on contig assembly
  ofstream myfile;
  myfile.open ("debuggingMatrix.dat", ios::app);


  
  /*
  for(i=0; i<numReads; i++)
    {

      //cout<<"on Read "<<i<<endl;

      for(j=0; j<readSeqs[i].length(); j++)
	{
	  if(i>0)
	    {
	      
	      matrix[i][startPositions[i-1]+readSize+j]=readSeqs[i][j];	      //Need to comment out to complile need to change
	    }else
	    {
	       matrix[i][readSize+j]=readSeqs[i][j];
	      
	    }
	}
    }
  */

  
  for(i=0; i<numReads; i++)
    {

      //cerr<<"on Read "<<i<<endl;
      
      int start=(numCol/2)-startPositions[i];
      int ctr=0;
      //cerr<<"start is "<<start<<endl;
      for(j=0; j<readSeqs[i].length(); j++)
	{
	      
	      matrix[i][j+start]=readSeqs[i][ctr];
	      //cerr<<readSeqs[i][ctr];
	      ctr++;
	}
      //cerr<<endl;
    }
    


  
  //printing the matrix
  //raw matrix before filtering


  /*
  for(i=0; i<numReads; i++)
    {
      for(j=0; j<numCol; j++)
	{
	  //matrix[i][((numCol/2)-startPositions[i])+j]=readSeq[i][j];
	
	  myfile<<matrix[i][j];
	}
      myfile<<endl;
    }
  
    
  myfile<<combinedNuc<<endl;


  myfile<<"finished printing out the matrix "<<endl;
 
  */


string validChar = "ACGTacgt";
  
flag=false;

 Nctr=0;

 badColCtr=0;

  for(i=0; i<numCol; i++)
    {
      
      nuc[0]=0;
      nuc[1]=0;
      nuc[2]=0;
      nuc[3]=0;

      nucCtr=0;
      
      for(j=0; j<numReads; j++)
	{
	  /*
	  if (validChar.find(matrix[j][i]) == std::string::npos) { 
	    continue;
	  }else
	    {
	      cout<<"nuc is "<<matrix[j][i]<<endl;
	    }
	  */

	   switch(matrix[j][i])
		    {
		    case 'A' :
		      nuc[0]=nuc[0]+1;
		      nucCtr++;
		      break;
		    
		    case 'a' :
		      nuc[0]=nuc[0]+1;
		      nucCtr++;
		      
		      break;
		   
		    case 'T' :
		      nuc[1]=nuc[1]+1;
		      nucCtr++;
		      break;
		    
		    case 't' :
		      nuc[1]=nuc[1]+1;
		      nucCtr++;
		      break;
		    
		    case 'G' :
		      nuc[2]=nuc[2]+1;
		      nucCtr++;
		      break;
		    
		    case 'g' :
		      nuc[2]=nuc[2]+1;
		      nucCtr++;
		      break;

		    case 'C' :
		      nuc[3]=nuc[3]+1;;
		      nucCtr++;
		      break;

		    case 'c' :
		      nuc[3]=nuc[3]+1;
		      nucCtr++;
		      break;
		      
		  
		    }


	}

      

      indexMax=-1;
      indexMax=distance(nuc.begin(), std::max_element(nuc.begin(), nuc.end()));

      if(nuc[indexMax] > cutoffMinNuc && ((nuc[indexMax]/nucCtr) > 0.75))
	{
	  flag=true;
	  
	  //cout<<nuc[indexMax]<<endl;
	  combinedNuc+=allNucs[indexMax];

	  /*	  
	  for(z=0; z<nuc.size(); z++)
	    {
	      cout<<nuc[z]<<" ";

	    }
	  cout<<endl;

	  cout<<"column is "<<i<<" adding  "<<allNucs[indexMax]<<" count is "<<nuc[indexMax]<<" index max is "<<indexMax<<endl;
	  */
	}else
	{

	  //if have already started the contig but not confident in next base call then place an N. 
	  if(flag)
	    {
	      

	      combinedNuc += 'N';
	      Nctr++;
	    }

	}


      if(((nuc[indexMax]/nucCtr) < 0.75))
	{
	  badColCtr++;
	}

    
    }
  

  

  posN=combinedNuc.find_last_not_of( 'N' ) +1;
  
  Nctr=Nctr-(combinedNuc.length()-posN);


//remove trailing Ns
  combinedNuc.erase( posN ); 


  /*

  if(combinedNuc.length() < 51)
    {
      return "0";
    }

  if((badColCtr/combinedNuc.length())>0.1)
    {
      return "0";
    }

  //if number of not confident nucleotide calls is greater than 10% of the contig throw it out
  if((Nctr/combinedNuc.length())>0.1)
    {
      combinedNuc="0";
    }
    
  */


  //cerr<<"reached this point in merge cluster "<<combinedNuc<<endl;
  
    contig=combinedNuc;

      
    myfile<<"####################Start of matrix###########################"<<endl;
             
    if(combinedNuc.length()>=10)
      { 
	for(i=0; i<numReads; i++)
	  {
	    for(j=0; j<numCol; j++)
	      {
	  //matrix[i][((numCol/2)-startPositions[i])+j]=readSeq[i][j];
	
		myfile<<matrix[i][j];
	      }
	    myfile<<endl;
	  }
    
	myfile<<endl;
	myfile<<endl;

	myfile<<"contig is "<<combinedNuc<<endl;
      }

    myfile<<endl;
    myfile<<"percent of bad columns is "<<(badColCtr/combinedNuc.length())<<endl;
    myfile<<endl;

    myfile.close();

  //return "test";
  return contig;

}

//give a vector of reads obtain the position of every kmer in every read. 
//Do not need to worry about reverse complement all reads should be on the 
//same strand enforced this earlier in code
void ReadCluster::getKmers()
{
  uint_fast64_t kmer, positionCtr, AbitString, CbitString, GbitString, TbitString; 
  long  i, j, k;
  

  string validChar = "ACGTacgtN";
  string DNAchar="ACGTacgt";


 //code to place binary represenation into right side of bit string. 
  //00 bitstring represents A
  AbitString=0*pow(2, 64-(2*clusterKmerSize))+0*pow(2, 64-((2*clusterKmerSize)-1));

    //11 bitstring represents T place into left side of bit string
  //TbitString=1*pow(2, 64-(2*clusterKmerSize))+1*pow(2, 64-((2*clusterKmerSize)-1));

  //place into right side of bitstring
  TbitString=1*pow(2, (2*clusterKmerSize)-2)+1*pow(2, ((2*clusterKmerSize)-1));
 
    //10 bitstring represents G places into left side of bit string
  //GbitString=0*pow(2, 64-(2*clusterKmerSize))+1*pow(2, 64-((2*clusterKmerSize)-1));
  
  //places into right side of bit string
  GbitString=0*pow(2, (2*clusterKmerSize))+1*pow(2, ((2*clusterKmerSize)-1));


  //01 bitstring represents C places into left side of bitstring
  //CbitString=1*pow(2, 64-(2*clusterKmerSize))+0*pow(2, 64-((2*clusterKmerSize)-1));

 //place into right side of bitstring
  CbitString=1*pow(2, (2*clusterKmerSize)-2)+0*pow(2, ((2*clusterKmerSize)));



 if(clusterKmerSize==0)
    {
      cerr<<"Should never get a cluster kmer  size of zero"<<endl;
      exit(EXIT_FAILURE);
    }


 
 dense_hash_map<uint_fast64_t, long, customHash> kmerInReadCounts; //hash table contains a count of how many reads contain a given kmer key is the kmer value is the count 
 kmerInReadCounts.set_empty_key(-1);
 kmerInReadCounts.resize(1000000);



   //revKey=0;
  kmer=0;
  positionCtr=0;

  //iterator through ever read in the cluster getting a list of the kmers 
  for(j=0; j<readSeqs.size(); j++)  
  {

    //kmerInReadCounts[j]=1;
    
    //    if(j%1000==0)
    //{
    //cerr<<"on Read "<<j<<endl;
	//}


    dense_hash_map<uint_fast64_t, long, customHash> readKmers; //hash table contains the kmers and their position in the read 
    readKmers.set_empty_key(-1);

    dense_hash_map<uint_fast64_t, long, customHash> alreadyPresent; //hash table keeps track of kmers that have already been found earlier in the read key is the kmer the value is a flag set to 1 indicating the 
    //kmer is present in the read

    alreadyPresent.set_empty_key(-1);


    //     std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> > kmerPositions;

   	
      positionCtr=0;
      //starting a new sequence
      for(i=readSeqs[j].length()-1; i>=0; i--)
	{
	  
	   //if character is not a nucleotide reset the bit string and start over
	  if(DNAchar.find(readSeqs[j][i]) == std::string::npos) {
		
	    //revKey=0;
		  kmer=0;
		  positionCtr=0;
		  //flag=true;
		  continue;
	      }else
		{
		  //right shift by 2 (C++ function) since added a single nucleotide and each nucleotide is represented by two bits.  
		  kmer=kmer>>2;
	      
		  switch(readSeqs[j][i])
		    {
		    case 'A' :
		      kmer=kmer|0;
		      break;
		    
		    case 'a' :
		      kmer=kmer|0;
		      break;
		   
		    case 'T' :
		      kmer=kmer|TbitString;
		      break;
		    
		    case 't' :
		      kmer=kmer|TbitString;
		      break;
 		    
		    case 'G' :
		      kmer=kmer|GbitString;
		      break;
		    
		    case 'g' :
		      kmer=kmer|GbitString;
		      break;

		    case 'C' :
		      kmer=kmer|CbitString;
		      break;

		    case 'c' :
		      kmer=kmer|CbitString;
		      break;
		      
		    default:
		      cerr<<kmer<<endl;
		      cerr<<"Base is not a nucleotide exiting "<<endl;
		      exit(EXIT_FAILURE);

		    }

		  positionCtr++;		  
		  
		}


	  
	  
	  //cerr<<"reached this point "<<positionCtr<<"\t"<<clusterKmerSize<<endl;
	   //not yet reached a full kmer in size
	      if((positionCtr<clusterKmerSize))
		{
		  continue;  
		}

	      
	      readKmers[kmer]=i;
	      
	      //lookup key if key not found will return iterator that points to end of hash table
	      auto iter=alreadyPresent.find(kmer);
	      
	      
	      if(iter==alreadyPresent.end()) //means kmer is not present earlier in the read this is important because only want to count a kmer once per read
		{
		  kmerInReadCounts[kmer]=kmerInReadCounts[kmer]+1;
		 
		  alreadyPresent[kmer]=1; //setting value to 1 indicating kmer has been found in this read
		}
	      
	}

      //cout<<"size of hash table in this function is "<<readKmers.size()<<endl;
      kmerPositions.push_back(readKmers);


  }

  
  //going through hash table which counts the number of times a kmer is found in the population of reads storing the hash table in a vector of pairs
  //doing this because can easily sort the vector 
  auto iterHash=kmerInReadCounts.begin();
  for(iterHash; iterHash!=kmerInReadCounts.end(); iterHash++)
    {
      pair<uint_fast64_t, long> keyPair(iterHash->first, iterHash->second);

      kmerCounts.push_back(keyPair);
    }


  sort(kmerCounts.begin(), kmerCounts.end(), comparePairs);

  /*
  cout<<"paired kmers are "<<endl;

  for(k=0;k<kmerCounts.size();k++)
        cout << kmerCounts.at(k).second << " (" << kmerCounts.at(k).first << "%)"<< endl; 

  cout<<"#############################\n\n\n\n\n";
  */
}


//very simple compare function to pass to sort in order to sort pairs by the second element
bool comparePairs(const std::pair<uint_fast64_t, long>&i, const std::pair<uint_fast64_t, long>&j)
{
    return i.second > j.second;
}


std::pair<uint_fast64_t, long> ReadCluster::getPair(long index) //given an index return the pair kmer count pair corresponding to that index
{

  return(kmerCounts[index]);

}

void ReadCluster::setKmerSize(int size)
{

  clusterKmerSize=size;

}

void ReadCluster::setStartPositions(uint_fast64_t kmer)
{
  

  auto iter=kmerPositions.begin();

    for(iter; iter!=kmerPositions.end(); iter++)
      {
	
	startPositions.push_back((*iter)[kmer]);

      }


}

void ReadCluster::printStartPositions()
{

  auto iter=startPositions.begin();

    for(iter; iter!=startPositions.end(); iter++)
      {
	cout<<(*iter)<<endl;


      }

}



//function to print out the kmerPositions member variable
//which is a vector of hash tables each hash table contains the position of each kmer in a read the 
//key is the kmer compressed as a 64 bit integer the value is the position on the read of
//the start of the kmer 

void ReadCluster::printKmerPositions()
{

  //std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> > kmerPositions; //a vector of hash tables each hash table contains the position of each kmer in a read the 

  

  for(auto iter=kmerPositions.begin(); iter!=kmerPositions.end(); iter++)
    {

      dense_hash_map<uint_fast64_t, long, customHash> readKmers;

      readKmers=(*iter);

      //cout<<"size of hash table is "<<readKmers.size()<<"\n";

      for(auto iterHash=readKmers.begin(); iterHash!=readKmers.end(); iterHash++)
	{
	  
	  cout<<"kmer is "<<iterHash->first<<"\t position is "<<iterHash->second<<"\n";
	  
	}
      
      cout<<"\n\n\n\n\n";


    } 


}
