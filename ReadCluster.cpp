#include "ReadCluster.h"

using google::dense_hash_map;
using google::sparse_hash_map;
using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::ofstream;
using std::ios;


ReadCluster::ReadCluster(){

  contig="0";
  //readSeq.push_back("test2");
  key2=-1;
  mainKey=false;
  refPositions.set_empty_key(-1);

}

int ReadCluster::getNumReads()
{

  return readSeq.size();
}


bool ReadCluster::allReadyPresent(string &read)
{
  int numReads=readSeq.size();
  int i;

  for(i=0; i<numReads; i++)
    {
      if(read.compare(readSeq[i])==0)
	{
	  return true;
	}

    }

  return false;
}



void ReadCluster::addPosWords(vector<uint_fast64_t> &newWords){
  posWords=newWords;
}

void ReadCluster::addNegWords(vector<uint_fast64_t> &newWords){
  negWords=newWords;
}





void ReadCluster::setStrand(char s){
  
  strand=s;

}



char ReadCluster::getStrand()
{

  return strand;

}

vector<string> ReadCluster::getSeqs(){
  
  return readSeq;

}

string ReadCluster::getContig()
{
  return contig;
}


bool ReadCluster::getMainKey()
{
  return mainKey;
}

uint_fast64_t ReadCluster::getKey2()
{
  return key2;

}

void ReadCluster::setMainKey(){

  mainKey=true;
}


void ReadCluster::setKey2(uint_fast64_t key)
{
  key2=key;
}

void ReadCluster::addSeq(string &read)
{
  readSeq.push_back(read);

}

void ReadCluster::printReads()
{

  int i;

  cout<<readSeq[i]<<endl;
  for(i=1; i<readSeq.size(); i++)
    {
      //readSeq[i].insert(startPositions[i], "  ");
      //cout<<readSeq[i]<<"\t"<<startPositions[i]<<endl;

      cout<<readSeq[i]<<endl;

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
  int numReads=readSeq.size();
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
      if(readSeq[i].length() > readSize)
	{
	  readSize=readSeq[i].length();
	}
      
    }

  numCol=(readSize*3);
  

			   
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


  

  for(i=0; i<numReads; i++)
    {

      //cout<<"on Read "<<i<<endl;

      for(j=0; j<readSeq[i].length(); j++)
	{
	  if(i>0)
	    {
	      
	      matrix[i][startPositions[i-1]+readSize+j]=readSeq[i][j];	      
	    }else
	    {
	      matrix[i][readSize+j]=readSeq[i][j];
	      
	    }
	}
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
      
	myfile<<combinedNuc<<endl;
      }
    
    myfile<<"percent of bad columns is "<<(badColCtr/combinedNuc.length())<<endl;


    myfile.close();

  //return "test";
  return contig;

}
