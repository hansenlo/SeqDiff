#include "ReadCluster.h"

#include <bitset>



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
using std::list;

/*
ReadCluster::ReadCluster(long numReads){

  cerr<<"reached this point number of Reads is "<<numReads<<endl;

  contig="0";
  //startPositions.assign(-1, numReads); //set intital start positions to -1
  
  startPositions(-1, 5); //set intital start positions to -1
  

  //readSeq.push_back("test2");

}
*/
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


void ReadCluster::addSequences(std::unordered_map<string, string> *allReads)
{

  auto hashIter=allReads->begin();

  for(hashIter; hashIter!=allReads->end(); hashIter++)
    {
      readSeqs.push_back(hashIter->first);
      qualityStrings.push_back(hashIter->second);
   }

 


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
      cout<<qualityStrings[i]<<endl;
      cout<<endl;
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
void ReadCluster::mergeReads(std::vector<std::string> &contigs, int kmerSize, int cutoffMinNuc, std::ofstream &debugging, vector<string> &clusterID)
{

  double Nctr, nucCtr, badColCtr;
  long int numCol;
  int numReads=readSeqs.size();
  int i, j, z, maxCtr, indexMax, posN, readSize;

  double startContig;

  string combinedNuc;
  string allNucs="ATGC";
  bool flag=false;

  vector<int> nuc(4, 0); //A T G C
 
  
  int qualityCutoff=20; //This is the base score quality cutoff if the base quality is less than this number than ignore that base

  //finding the maximum read size
  readSize=0;
  for(i=0; i<numReads; i++)
    {
      if(readSeqs[i].length() > readSize)
	{
	  readSize=readSeqs[i].length();
	}
      
    }

  numCol=(readSize*3);
  
  vector<int> startMatrix(numReads, -1); //starting index in matrix of every assembled read

			   
  //matrix to hold reads
  vector<vector<char>>matrix(numReads, vector<char>(numCol, 'N'));

  //matrix to hold qualityScores
  //vector<vector<long>>qualityMatrix(numReads, vector<long>(numCol, -1));


  //cerr<<"size of matrix is "<<matrix.size()<<" number of cols is "<<matrix[0].size()<<endl;

 
 
  //file to put debugging info on contig assembly
  //ofstream myfile;
  //myfile.open ("debuggingMatrix.dat", ios::app);

  //assembling reads that all share the most common kmer should take care of the majority of reads
  for(i=0; i<numReads; i++)
    {
     
      if(startPositions[i]==-1) //only assemble those reads that share a kmer
      {
        continue;
      }

      int start=(numCol/2)-startPositions[i];
      int ctr=0;
      //cerr<<"start is "<<start<<endl;
      startMatrix[i]=start; //storing the start position in the matrix of the assembled read


      for(j=0; j<readSeqs[i].length(); j++)
	{
	  int quality=int(qualityStrings[i][ctr])-33;
	  //cerr<<"quality character is "<<qualityStrings[i][ctr]<<endl;

	  /*
	  if(clusterID==373437)
	    {
	      cerr<<quality<<" "<<qualityStrings[i][ctr];
	    }
	  */
	  

	  if(quality>=0 && quality<=qualityCutoff) 
	    {
	      ctr++;
	      //cerr<<"inside this condition "<<endl;
	      continue; //value in matrix will be left as an N i.e. that base will be ignored

	    }else
	    {
	      matrix[i][j+start]=readSeqs[i][ctr]; //if a good quality base than use that bases nucleotide value	 
	    }

	      //cerr<<readSeqs[i][ctr];
	      ctr++;
	}


      /*
      if(clusterID==373437)
	{
	  cerr<<endl;
	  cerr<<"end of processing the read "<<endl;

	}
      */


//cerr<<endl;
    }
    


double percentBadCol, percentNs;
 
/*
  
 std::vector<int> test(matrix.size());
 std::iota(test.begin(), test.end(), 0); //for the initial assembly use all reads
  

 combinedNuc=assembleContig(matrix, test, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig);

 checkContig(combinedNuc, readSize, badColCtr, Nctr, percentBadCol, percentNs);


printMatrix(matrix, test, combinedNuc, percentBadCol, clusterID[0], debugging, percentNs); //debugging only print alignment Matrix

 return;
*/

  char temp;
  int k;
  uint_fast64_t kmer;  
  int ctrUnassembled=1;
  bool foundFlag=false;
   
  //Now adding in reads that were not assembled in the previous step
  //algorithim: For a read that was not assembled look for a kmer in common with a read that was assembled 
  //line up the unassembled read with the read that is assembled based on the kmer in common

  //while(ctrUnassembled!=0)
  //{

      ctrUnassembled=0;
      for(i=0; i!=startPositions.size(); ++i) //go through every read looking for unassembled reads
	{
	  foundFlag=false;

      //int temp=0;

	  if(startPositions[i]==-1) //read was not assembled
	    {
	      ctrUnassembled++;

	  //go through every kmer in the unassembled read
	      auto HashIter=kmerPositions[i].begin();
	      for(HashIter; HashIter!=kmerPositions[i].end(); HashIter++)
		{
		  if(foundFlag) //if have already found a matching kmer break out of loop
		    {
		      break;
		    }
	      
		  kmer=HashIter->first;
	      

		  //check to see if the kmer was present more than once in any read if so do not use that kmer
		  if(presentMultipleTimes.count(kmer)>0)
		    {
		      continue;
		    }

		  /*
		  if(clusterID[0]=="25411")
		    {
		      cerr<<"all multiple kmers are "<<endl;
		      
		      for(auto iter=presentMultipleTimes.begin(); iter!=presentMultipleTimes.end(); iter++)
			{
			  cerr<<int2String(iter->first, 25)<<endl;
			}
		      

		      cerr<<"kmer is "<<kmer<<endl;
		      cerr<<"string representation is "<<int2String(kmer, 25)<<"\n\n\n";
		    }
		  */

	      //cerr<<"for read "<<i<<" on kmer "<<kmer<<" kmer is "<<temp<<endl;
     
	      //temp++;

		  foundFlag=false;
	      
		  for(j=0; j!=kmerPositions.size(); ++j) //go through every read looking for ones that are assembled 
		    {
		      if(startPositions[j] > -1) //if a read has been previously assembled check to see if it has a kmer in common with the unassembled read
			{
			  auto findIter=kmerPositions[j].find(kmer); //look for kmer in hash table
		      
			  if(findIter!=kmerPositions[j].end()) //kmer is present in an assembled read
			    {

			      //get a read that is assembled to compare against
			     
			      int index=0;

			      while(startPositions[index]<0 || index == j)
				{
				  index++;
				}

			      int numMismatches=0;
			      int sizeAligned=0;
			      
			      //if only one read has been assembled than do not check that one read against a neighbor because it will not have any neighbors
			      if(index<kmerPositions.size())
				{

			      //check to make sure assembled read matches well with a neighboring read if it has to many mismatches do not use it
				  for(int y=0; y<readSeqs[j].size(); y++)
				    {
				      if(matrix[j][startMatrix[j]+y]!='N' && matrix[index][startMatrix[j]+y]!='N') 
					{
					  sizeAligned++;
					  if(matrix[j][startMatrix[j]+y]!=matrix[index][startMatrix[j]+y])
					    {
					      numMismatches++;
					    }

					}

				    }
				}else
			      {
				  numMismatches=0;
			      }
			      
			      /*
			       if(clusterID==228237)
				{
				  cerr<<"index of comparing read is "<<index<<"  index of assembled read is "<<j<<endl;
				  cerr<<"number of mismatches is "<<numMismatches<<"  size Aligned is "<<sizeAligned<<endl;
				  exit(0);
				}
			      */

			      if((double(numMismatches)/double(sizeAligned))>0.1) //if more than 10% of read mismatches do not use that read as a template to aligne other reads
				{
				  continue;
				}
			      

			      int start=(startMatrix[j]+kmerPositions[j][kmer])-kmerPositions[i][kmer]; //line up the start positions of the kmer in common between the assembled and unassembled read 
			      int ctr=0;
	
			      /*
			      if(clusterID==228237)
				{
				  cerr<<" assembled read is "<<j<<" index of unassembled read is "<<i<<endl;
				  cerr<<"kmer in assembled Read is "<<bit2String(kmer, kmerSize)<<" string extracted from read is "<<readSeqs[j].substr(kmerPositions[j][kmer], kmerSize)<<endl;
				  cerr<<"kmer in unassembled Read is "<<bit2String(kmer, kmerSize)<<" string extracted from unassembled read is "<<readSeqs[i].substr(kmerPositions[i][kmer], kmerSize)<<endl;

				}
			      */
		  
			      /*
			      if(readSeqs[i].compare("CACCAATATGGCACATGTATACATATGTAACAAACCTGCACGTTGTGCACATGTACCCTAGAACTTAAAGTATAATGAAAAAAAAAAGCAATATAGATCGG")==0)
				{
				  //				  cout<<"start position is "<<start<<" startMatrix is "<<startMatrix[j]<<" kmer position in assembled read is "<<kmerPositions[j][kmer]<<" kmer position in unassembled read is "<<kmerPositions[i][kmer]<<endl;
				  //cout<<" assembled read is "<<j<<" index of unassembled read is "<<i<<endl;
				  
				  cout<<"kmer in assembled Read is "<<bit2String(kmer, kmerSize)<<" string extracted from read is "<<readSeqs[j].substr(kmerPositions[j][kmer], kmerSize)<<endl;
				  cout<<"kmer in unassembled Read is "<<bit2String(kmer, kmerSize)<<" string extracted from unassembled read is "<<readSeqs[i].substr(kmerPositions[i][kmer], kmerSize)<<endl;
				  
				  //exit(0);

				}
			      */


			  //cerr<<"for read "<<i<<" on kmer "<<kmer<<" matching kmer position is "<<kmerPositions[j][kmer]<<" start is "<<start<<endl;
     
			  //cerr<<"for read "<<i<<"start is "<<start<<endl;
     
			  

			      if(start<0) //if read extends off the beginning of the matrix dynamically insert new columns at the beginning of the matrix very inefficient operation :(
				{
				  //cerr<<"reallocating at the beginning "<<endl;
				 
				 
				  start=start+(readSize*2);

				  int z;
				  for(z=0; z<matrix.size(); z++)
				    {
				      matrix[z].insert(matrix[z].begin(), readSize*2, 'N'); //adding Ns at the beginning of the matrix adding 2X a read number of Ns
				      //qualityMatrix[z].insert(qualityMatrix[z].begin(), readSize*2, -1); //adding Ns at the beginning of the quality matrix same as read matrix

				      startMatrix[z]=startMatrix[z]+(2*readSize);
				    }


				}
			      
			      if((start+readSize)>matrix[i].size()) //if read extends off the end of the matrix dynamically insert new columns at the end of the matrix 
				{
				  //cerr<<"reallocating at the end "<<endl;
				 

				  for(z=0; z<matrix.size(); z++)
				    {
				      matrix[z].insert(matrix[z].end(), readSize*2, 'N'); //adding Ns at the end of the matrix adding 2X a read number of Ns
				      //qualityMatrix[z].insert(qualityMatrix[z].end(), readSize*2, -1); //adding Ns at the end of the quality matrix adding 2X a read number of Ns

				    }
				
				}				

			  
			  //add the unassembled read to the assembled reads
			      for(k=0; k<readSeqs[i].length(); k++)
				{

		

			      //matrix[i][k+start]='A';
				  int quality=int(qualityStrings[i][ctr])-33;
				  if(quality>=0 && quality<=qualityCutoff) 
				    {
				      ctr++;
				      continue; //value in matrix will be left as an N i.e. that base will be ignored
				      
				    }else
				    {
				      matrix[i][k+start]=readSeqs[i][ctr]; //if a good quality base than use that bases nucleotide value	 
				    }
				
			      
				  ctr++;
				}
			  			  

			      startMatrix[i]=start;
			      startPositions[i]=kmerPositions[i][kmer]; //adding the new assembled read to the set of reads already assembled
			      foundFlag=true;
			      break;
			    }
		      

			}
		    }
	      
	      
		}	  
	  
	    }

	}
    
      //}



      double longestDistN=0;
      //Assembling the contig from the matrix of aligned reads
      std::vector<int> rowsToAssemble(matrix.size());
      std::iota(rowsToAssemble.begin(), rowsToAssemble.end(), 0); //for the initial assembly use all reads

      int subClusterCtr=0;

      vector<int> problematicCols;

      combinedNuc=assembleContig(matrix, rowsToAssemble, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr, problematicCols);
  
      int SubClusterFlag=subClusterCtr;

      /*
#pragma omp critical(DEBUGGING_CLUSTER)
  {  
  //printing the matrix
  //raw matrix before filtering
  
    numCol=matrix[0].size(); //all rows are the same size  
    for(i=0; i<numReads; i++)
      {
	for(j=0; j<numCol; j++)
	  {
	  //matrix[i][((numCol/2)-startPositions[i])+j]=readSeq[i][j];
	
	    debugging<<matrix[i][j];
	  }
	debugging<<endl;
      }

    debugging<<"cluster being assembled is  "<<clusterID<<endl;
    debugging<<"contig is "<<combinedNuc<<endl;
    debugging<<"percentage of bad columns is "<<badColCtr/combinedNuc.length()<<endl;
    debugging<<"finished printing out the matrix "<<endl;
    debugging<<"##############################################\n\n\n\n"<<endl;
  

  }
      */

      
      bool goodContig=false;
      
       goodContig=checkContig(combinedNuc, readSize, badColCtr, Nctr, percentBadCol, percentNs);


       string printableClusterID=clusterID[0];


      if(goodContig==true)
	{
	  contigs.push_back(combinedNuc);

	  if(subClusterCtr>2)
	    {
	      bool subCluster=false;
	      subCluster=assembleSubClusters(matrix, rowsToAssemble, combinedNuc, contigs, startContig, cutoffMinNuc, readSize, problematicCols);
	      
	      
	      if(subCluster==true)
		{
		  clusterID.push_back(clusterID[0]+"_"+"subCluster");
		}
	      
	    }

	}else
	{
	  contigs.push_back("0");

	  
	  printMatrix(matrix, rowsToAssemble, combinedNuc, percentBadCol, printableClusterID, debugging, percentNs);
	}


      int oldSize;

      if(percentNs>0 && percentNs<0.03 && goodContig==true)
	{	  
	  oldSize=contigs.size();
	      
	  extractHetro(matrix, rowsToAssemble, combinedNuc, contigs, startContig, cutoffMinNuc);
	  
	  if(contigs.size()==oldSize+1)
	    {
	      clusterID.push_back(clusterID[0]+"_"+"0_hetro");
	    }
	  
	  if(contigs.size()==oldSize+2)
	    {
	
	      clusterID.push_back(clusterID[0]+"_"+"0_hetro");
	      clusterID.push_back(clusterID[0]+"_"+"1_hetro");

	    }
	}

    
      //printMatrix(matrix, rowsToAssemble, combinedNuc, percentBadCol, printableClusterID, debugging, percentNs); //debugging only print alignment Matrix


      


      //if the contig is to noisy look to break down the cluster into two different variants
      //if(goodContig==false && ((longestDistN/combinedNuc.size())>0.3)  && percentNs>0.03 && combinedNuc.size()>=readSize)
	  if(goodContig==false  && percentNs>0.03 && combinedNuc.size()>=readSize)
	    {
	      oldSize=contigs.size();
	  
	      //cout<<"Nctr is "<<Nctr<<" "<<combinedNuc<<" percentNs is "<<percentNs<<endl;


	      extractVariants(matrix, rowsToAssemble, combinedNuc, contigs, startContig, cutoffMinNuc, readSize);
	  	 

	  if(contigs.size()>oldSize)
	    { 

	      
	      //if(clusterID[0]=="24459")
	      //{
	      //      cerr<<"successfully extracted variants"<<endl;
	      //}



	      clusterID.push_back(clusterID[0]+"_"+"0_splitVariant");
	      clusterID.push_back(clusterID[0]+"_"+"1_splitVariant");

	    }

	}

  if(goodContig==false) //if the contig is noisy it may be the cluster is impure try and find subclusters that assembly better
    {
      if(numReads>(2*cutoffMinNuc)) //making sure can break up the cluster into 2 or more clusters with the mininum number of reads to pass the cutoff
	{
	  //cerr<<"debugging contig is "<<combinedNuc<<endl;	  

	  //generating the matrix of distances between every 2 aligned sequences in the alignmentMatrix
	  //The distance will be stored only if the number of aligned bases meet a certain standard otherwise it will be set to 0
	  //vector<vector<int>> numMatches(numReads, vector<int>(numReads, -1)); //number of basepair matches between 2 sets of aligned sequences

	  //vector<vector<int>>numMatches(numReads, vector<int>(numReads, -1));

	  Graph graphMatches(numReads);

	  getDistanceGraph(matrix, startMatrix, graphMatches, readSize);
	   
	  vector<bool> visited(numReads, false);
	  vector<int> cluster;
	  
	  int clusterCtr=0;
	  for(int d=0; d<numReads; d++)
	    {
	      if(visited[d])
		{
		  continue;
		}else
		{

		  cluster.clear();
		   //vector to hold the given cluster
		  //cerr<<"getting the subcluster starting with read "<<d<<endl;
		  graphMatches.DFS(d, cluster, visited); //print the subcluster that includes read 0
	
		  /*
		  int z;
		  for( z=0; z<cluster.size(); z++)
		    {
		      cerr<<cluster[z]<<" ";
    
		    }
		  cerr<<endl;
		  */

		}

	      if(cluster.size() > cutoffMinNuc)
		{

		  clusterCtr++;

		  combinedNuc=assembleContig(matrix, cluster, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr, problematicCols);
		  
		  
		  
		  //cerr<<"contig of subcluster is "<<combinedNuc<<endl;

		  goodContig=false;
		  goodContig=checkContig(combinedNuc, readSize, badColCtr, Nctr, percentBadCol, percentNs);
	      
		  string subClusterID=clusterID[0]+"_"+std::to_string(clusterCtr);

		  //printMatrix(matrix, cluster, combinedNuc, percentBadCol, subClusterID, debugging, percentNs); //debugging only print alignment Matrix
	  
		  if(goodContig==true)
		    {
		  //adding the new contig to the set of contigs
		      contigs.push_back(combinedNuc); 
		      clusterID.push_back(subClusterID);
		    }else
		    {
		      contigs.push_back("0"); 
		      clusterID.push_back(subClusterID);

		      printMatrix(matrix, cluster, combinedNuc, percentBadCol, subClusterID, debugging, percentNs); //debugging only print alignment Matrix
	  

		    }
		}

	    }
	  
	  /*
	  cerr<<"got distance matrix starting sub clusters"<<endl;

	  vector<vector<int>> clusters; //matrix to hold all the subclusters

	  getSubClusters(numMatches, clusters); //extracting subclusters from the data
	  
	  cerr<<"number of total reads is "<<numReads<<endl;

	  cerr<<"number of clusters is "<<clusters.size()<<endl;

	  for(i=0; i<clusters.size(); i++)
	    {
	      combinedNuc=assembleContig(matrix, clusters[i], nucCtr, badColCtr, cutoffMinNuc, Nctr);

	      //cerr<<"contig of subcluster "<<i<<" is "<<combinedNuc<<endl;

	      checkContig(combinedNuc, readSize, badColCtr, Nctr, percentBadCol, percentNs);
	      
	      string subClusterID=std::to_string(clusterID)+std::to_string(i);

	      printMatrix(matrix, clusters[i], combinedNuc, percentBadCol, subClusterID, debugging); //debugging only print alignment Matrix

	
	    }

	  */
	 


	  //exit(0);
	}
    
    }




//Given the alignment matrix of reads and the row index and start positition of the 2 reads calculate percent
//difference of the 2 reads for the aligned regions row1 and row2 are the column indexes for the two reads, start1 and start2 are the 
//the start columns of the reads in the alignment matrix sizeAligned is the number of aligned bases between the 2 reads the aligned bases may not necassarily match
//function returns the number of matches
//uint_fast64_t ReadCluster::numberDiff(vector<vector<char>> &alignmentMatrix, int row1, int row2, int start1, int start2, int &sizeAligned)

 





  //cerr<<"reached this point in merge cluster "<<combinedNuc<<endl;
  
  //  contig=combinedNuc;

    /*
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
    */

  //return "test";
  //return contig;

}


 //checks to see if columns that are Ns are approximately 50% to different bases if so will return contigs for both cases if successfull function will return a 1 and a vector of new contigs
bool ReadCluster::extractHetro(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &contigs, int startContig, int cutoffMinNuc)
  {
    
    
    bool flag;
    int numCol, i, j, k, sizeReadsToUse, indexMax, posN;
    vector<double> nuc(4, 0); //A T G C
    string allNucs="ATGC";

    double nucCtr;

    string validChar = "ACGTacgt";
  
    vector<int> positionsN; // holds all the positions that sub occurs within str

    //find allNs in a contig
    double pos = combinedNuc.find("N", 0);
    while(pos != string::npos)
      {
	positionsN.push_back(pos);
	pos = combinedNuc.find("N",pos+1);
      }


    int currentRow=0; //index into the current rows to assemble vector


    numCol=alignmentMatrix[0].size(); //all rows are the same size
    sizeReadsToUse=rowsToAssemble.size();

   

    flag=false;
    for(i=0; i<positionsN.size(); i++)
      {
      
	//only look at Ns that sit well inside the contig and are not close to either end
	if( (((positionsN[i])/double(combinedNuc.size())) < 0.3) || (((positionsN[i])/double(combinedNuc.size())) > 0.7))
	  {
	    continue;
	  }

	
	nuc[0]=0;
	nuc[1]=0;
	nuc[2]=0;
	nuc[3]=0;

	nucCtr=0;
      
      //for every column run through all the rows that will be used in this assembly
	for(k=0; k<sizeReadsToUse; k++)
	{
	  j=rowsToAssemble[k]; //j is row index to use in assembly
  
	  switch(alignmentMatrix[j][positionsN[i]+startContig])
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

      


	//indexMax=-1;
      //indexMax=distance(nuc.begin(), std::max_element(nuc.begin(), nuc.end()));

      //sort in descending order
      //std::sort(nuc.rbegin(), nuc.rend()); 
      
	//sort the nuc vector but retain the index of the sorted vector
      vector<size_t> idx(nuc.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
	   [&nuc](size_t i1, size_t i2) {return nuc[i1] < nuc[i2];});


      /*
      for(int z=0; z<nuc.size(); z++)
	{
	  cout<<nuc[z]<<endl;
	}

      cout<<"end of nuc count \n"<<endl;

      for(int z=0; z<idx.size(); z++)
	{
	  cout<<idx[z]<<endl;
	}	

      cout<<"end if index count \n\n\n\n";
      */



      //check to make sure Ns are highly covered and that the nucleotides with the second and first highest counts are very close to fifty percent
      if(nuc[idx[idx.size()-1]]>(2*cutoffMinNuc))
	{

	  if(((nuc[idx[idx.size()-1]]/nucCtr)>0.4) && ((nuc[idx[idx.size()-1]]/nucCtr)<0.6) )
	    {
	     
	      if(((nuc[idx[idx.size()-2]]/nucCtr)>0.4) && ((nuc[idx[idx.size()-2]]/nucCtr)<0.6) )
		{

		  /*
		   for(int z=0; z<nuc.size(); z++)
		     {
		       cout<<nuc[z]<<endl;
		     }

		   cout<<"end of nuc count \n"<<endl;
		   
		   for(int z=0; z<idx.size(); z++)
		     {
		       cout<<idx[z]<<endl;
		     }	
		   
		   cout<<"end of index count \n\n\n\n";

		  */

		  if(flag==false)
		    {

		      contigs.push_back(combinedNuc);
		      contigs.push_back(combinedNuc);
		  //co
		  //cout<<"N position is "<<newContigs[newContigs.size()-1][positionsN[i]]<<endl;
		      
		      flag=true;
		    }

		  
		  if(flag==true)
		    {
		      contigs[contigs.size()-1][positionsN[i]]=allNucs[idx[idx.size()-1]];
		      contigs[contigs.size()-2][positionsN[i]]=allNucs[idx[idx.size()-2]];

		      //		      cerr<<"replacing it with "<<allNucs[idx[idx.size()-1]]<<" "<<allNucs[idx[idx.size()-2]]<<endl;
		    }

		  
		  
	      //newContigs.push_back(combinedNuc);
	      //  newContigs[positionsN[i]]=allNucs[idx[idx.size()-1]];

						    
		}

	    }

	}
      




      }
  
    


  }


bool ReadCluster::extractVariants(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &contigs, double &startContig, int cutoffMinNuc, int readSize)
{
    
    
    bool flag;
    int numCol, i, j, k, sizeReadsToUse, indexMax, posN;
    vector<double> nuc(4, 0); //A T G C
    string allNucs="ATGC";

    double nucCtr;

    string validChar = "ACGTacgt";
  
    vector<int> positionsN; // holds all the positions that sub occurs within str

    //find allNs in a contig
    double pos = combinedNuc.find("N", 0);
    while(pos != string::npos)
      {
	positionsN.push_back(pos);
	pos = combinedNuc.find("N",pos+1);
      }


    double smallestDist=10000000;
    int indexSmallestDist=-1;
    //go through the Ns looking for the one closest to the middle of the contig
    for(i=0; i<positionsN.size(); i++)
      {
	if(abs(static_cast<double>(positionsN[i]-(combinedNuc.size()/2))) < smallestDist )
	  {
	    smallestDist=abs(static_cast<double>(positionsN[i]-(combinedNuc.size()/2)));
	    indexSmallestDist=i;
	  }
      }

    if(indexSmallestDist==-1)
      {
	return(false);
      }


    numCol=alignmentMatrix[0].size(); //all rows are the same size
    sizeReadsToUse=rowsToAssemble.size();

   

    flag=false;
      

	
	nuc[0]=0;
	nuc[1]=0;
	nuc[2]=0;
	nuc[3]=0;

	nucCtr=0;
      
      //run through the 
	for(k=0; k<sizeReadsToUse; k++)
	{
	  j=rowsToAssemble[k]; //j is row index to use in assembly
  
	  switch(alignmentMatrix[j][positionsN[indexSmallestDist]+startContig])
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
       

      



      
	//sort the nuc vector but retain the index of the sorted vector
      vector<size_t> idx(nuc.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
	   [&nuc](size_t i1, size_t i2) {return nuc[i1] < nuc[i2];});


      /*
      for(int z=0; z<nuc.size(); z++)
	{
	  cout<<nuc[z]<<endl;
	}

      cout<<"end of nuc count \n"<<endl;

      for(int z=0; z<idx.size(); z++)
	{
	  cout<<idx[z]<<endl;
	}	

      cout<<"end if index count \n\n\n\n";
      */



      //check to make sure the second highest nucleotide in the column is greater than the min necassary to make a cluster
      //if so run through the column and seperate rows into the first and second highest matching nucleotides then try and assemble contigs based on the two seperate 
      //groupings
      if(nuc[idx[idx.size()-1]]>(cutoffMinNuc))
	{

	  double Nctr, nucCtr, badColCtr, longestDistN, percentNs, percentBadCol;
	  int subClusterCtr; //this variable is not used in this call to assembleContig

	  vector<int> problematicCols; //this variable is not used just using it so can call the assembleContig function

	  std::vector<int> rowsCluster1, rowsCluster2;
	  bool goodContig=false;

	  for(k=0; k<sizeReadsToUse; k++)
	    {
	      j=rowsToAssemble[k]; //j is row index to use in assembly
  
	      if(alignmentMatrix[j][positionsN[indexSmallestDist]+startContig]==allNucs[idx[idx.size()-1]])
		{
		  rowsCluster1.push_back(j);
		}

	      if(alignmentMatrix[j][positionsN[indexSmallestDist]+startContig]==allNucs[idx[idx.size()-2]])
		{
		  rowsCluster2.push_back(j);
		}

	    }  

	  string contig1=assembleContig(alignmentMatrix, rowsCluster1, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr, problematicCols);
  
	  goodContig=false;
	  
	  goodContig=checkContig(contig1, readSize, badColCtr, Nctr, percentBadCol, percentNs);


      
	       if(goodContig==true)
		 {
		   contigs.push_back(contig1);
		   
		 }
       
	       if(nuc[idx[idx.size()-2]]>(cutoffMinNuc))
		 {
	       
		   string contig2=assembleContig(alignmentMatrix, rowsCluster2, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr, problematicCols);
	       
		   goodContig=false;

		   goodContig=checkContig(contig2, readSize, badColCtr, Nctr, percentBadCol, percentNs);
	       

      
		   if(goodContig==true)
		     {
		       contigs.push_back(contig2);
		   
		     }
		 }


	    

	}
	  /*
		  if(flag==false)
		    {

		      contigs.push_back(combinedNuc);
		      contigs.push_back(combinedNuc);
		  //co
		  //cout<<"N position is "<<newContigs[newContigs.size()-1][positionsN[i]]<<endl;
		      
		      flag=true;
		    }

		  
		  if(flag==true)
		    {
		      contigs[contigs.size()-1][positionsN[i]]=allNucs[idx[idx.size()-1]];
		      contigs[contigs.size()-2][positionsN[i]]=allNucs[idx[idx.size()-2]];

		      //		      cerr<<"replacing it with "<<allNucs[idx[idx.size()-1]]<<" "<<allNucs[idx[idx.size()-2]]<<endl;
		    }
	  */

      
  
    


}


bool ReadCluster::assembleSubClusters(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, std::string &combinedNuc, std::vector<std::string> &contigs, double &startContig, int cutoffMinNuc, int readSize, vector<int> &problematicCols)
{
    
    
    bool flag;
    int numCol, i, j, k, sizeReadsToUse, indexMax, posN;
    vector<double> nuc(4, 0); //A T G C
    string allNucs="ATGC";

    double nucCtr;

    string validChar = "ACGTacgt";
  
 
    vector<vector<int>> potentialSubClusters; //a vector of vectors each row is the set of reads that potentially belong to the same cluster every number in a row is a row index into the alignmentMatrix

    vector<int> intersection;

    numCol=alignmentMatrix[0].size(); //all rows are the same size
    sizeReadsToUse=rowsToAssemble.size();

    vector<int> dominantCluster;

    flag=false;
      

    for(i=0; i<problematicCols.size(); i++)
      {
	
	nuc[0]=0;
	nuc[1]=0;
	nuc[2]=0;
	nuc[3]=0;
	
	nucCtr=0;
      
	potentialSubClusters.push_back(vector<int>());

      //run through the 
	for(k=0; k<sizeReadsToUse; k++)
	{
	  j=rowsToAssemble[k]; //j is row index to use in assembly
  
	  switch(alignmentMatrix[j][problematicCols[i]])
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
       

	

	//sort the nuc vector but retain the index of the sorted vector
	vector<size_t> idx(nuc.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
	     [&nuc](size_t i1, size_t i2) {return nuc[i1] < nuc[i2];});


	//run through all the rows looking to see if a given row mismatches the consensus
	for(k=0; k<sizeReadsToUse; k++)
	  {
	    j=rowsToAssemble[k]; //j is row index to use in assembly
  
	   
		if(alignmentMatrix[j][problematicCols[i]]==allNucs[idx[idx.size()-2]])
		  {
		    potentialSubClusters[i].push_back(j);
		  }	
	    
	  }  



	//if have stored at least two problematic columns than check to see if they have similar rows if so then assemble two seperate clusters
	if(i==1)
	  {
	    sort(potentialSubClusters[0].begin(), potentialSubClusters[0].end());  
	    sort(potentialSubClusters[1].begin(), potentialSubClusters[1].end());
	    
	    set_intersection(potentialSubClusters[0].begin(), potentialSubClusters[0].end(), potentialSubClusters[1].begin(), potentialSubClusters[1].end(), back_inserter(intersection)); 

	    //if the rows that mismatch the consensus all agree with each other than evidence that two variants have been placed in the same cluster
	    if(intersection.size()>(2*cutoffMinNuc))
	      {
		
		 double Nctr, nucCtr, badColCtr, longestDistN, percentNs, percentBadCol;
		 int subClusterCtr; //this variable is not used in this call to assembleContig
		
		 bool goodContig=false;

		 vector<int> problematicCols;
		
		 string contig1=assembleContig(alignmentMatrix, potentialSubClusters[0], nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr, problematicCols);
  
		 goodContig=checkContig(contig1, readSize, badColCtr, Nctr, percentBadCol, percentNs);
      
		 if(goodContig==true)
		  {
		    contigs.push_back(contig1);
		     
		    return(true);

		  }else
		   {
		     return(false);
		   }
		
	      }
	    
	    

	  }

      }



    return(false);



      //      potentialSubClusters

      /*
      for(int z=0; z<nuc.size(); z++)
	{
	  cout<<nuc[z]<<endl;
	}

      cout<<"end of nuc count \n"<<endl;

      for(int z=0; z<idx.size(); z++)
	{
	  cout<<idx[z]<<endl;
	}	

      cout<<"end if index count \n\n\n\n";
      */

      /*

      //check to make sure the second highest nucleotide in the column is greater than the min necassary to make a cluster
      //if so run through the column and seperate rows into the first and second highest matching nucleotides then try and assemble contigs based on the two seperate 
      //groupings
      if(nuc[idx[idx.size()-1]]>(cutoffMinNuc))
	{

	  double Nctr, nucCtr, badColCtr, longestDistN, percentNs, percentBadCol;
	  int subClusterCtr; //this variable is not used in this call to assembleContig

	  std::vector<int> rowsCluster1, rowsCluster2;
	  bool goodContig;

	  for(k=0; k<sizeReadsToUse; k++)
	    {
	      j=rowsToAssemble[k]; //j is row index to use in assembly
  
	      if(alignmentMatrix[j][positionsN[indexSmallestDist]+startContig]==allNucs[idx[idx.size()-1]])
		{
		  rowsCluster1.push_back(j);
		}

	      if(alignmentMatrix[j][positionsN[indexSmallestDist]+startContig]==allNucs[idx[idx.size()-2]])
		{
		  rowsCluster2.push_back(j);
		}

	    }  

	  string contig1=assembleContig(alignmentMatrix, rowsCluster1, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr);
  
	       goodContig=checkContig(contig1, readSize, badColCtr, Nctr, percentBadCol, percentNs);


      
	       if(goodContig==true)
		 {
		   contigs.push_back(contig1);
		   
		 }
       
	       if(nuc[idx[idx.size()-2]]>(cutoffMinNuc))
		 {
	       
		   string contig2=assembleContig(alignmentMatrix, rowsCluster2, nucCtr, badColCtr, cutoffMinNuc, Nctr, startContig, longestDistN, subClusterCtr);
	       
		   checkContig(contig2, readSize, badColCtr, Nctr, percentBadCol, percentNs);
	       

      
		   if(goodContig==true)
		     {
		       contigs.push_back(contig2);
		   
		     }
		 }


	    

	}

      */
	  /*
		  if(flag==false)
		    {

		      contigs.push_back(combinedNuc);
		      contigs.push_back(combinedNuc);
		  //co
		  //cout<<"N position is "<<newContigs[newContigs.size()-1][positionsN[i]]<<endl;
		      
		      flag=true;
		    }

		  
		  if(flag==true)
		    {
		      contigs[contigs.size()-1][positionsN[i]]=allNucs[idx[idx.size()-1]];
		      contigs[contigs.size()-2][positionsN[i]]=allNucs[idx[idx.size()-2]];

		      //		      cerr<<"replacing it with "<<allNucs[idx[idx.size()-1]]<<" "<<allNucs[idx[idx.size()-2]]<<endl;
		    }
	  */

      
  
    


}






 //Given the alignment matrix of reads and the row index and start positition of the 2 reads calculate percent
//difference of the 2 reads for the aligned regions row1 and row2 are the column indexes for the two reads, start1 and start2 are the 
//the start columns of the reads in the alignment matrix sizeAligned is the number of aligned bases between the 2 reads the aligned bases may not necassarily match
//function returns the number of matches
uint_fast64_t ReadCluster::numberDiff(vector<vector<char>> &alignmentMatrix, int row1, int row2, int start1, int start2, int &sizeAligned, int readSize)
{
  int startShared, i, j, numCols;
  int numMatches;

  sizeAligned=0;

  if(start1>start2)
    {
      startShared=start1;

    }else
    {
      startShared=start2;
    }

  int endAlignmentIndex;

  if((startShared+readSize) > alignmentMatrix[0].size())
    {

      endAlignmentIndex=alignmentMatrix[0].size();

    }else
    {
      endAlignmentIndex=startShared+readSize;
    }


  /*
  cerr<<"row1 and row2 are "<<row1<<" "<<row2<<endl;
  cerr<<"start1 and start2 are "<<start1<<" "<<start2<<" startShared is "<<startShared<<endl;
  cerr<<"number of columns in matris is "<<alignmentMatrix[0].size()<<endl;
  */

  numMatches=0;

  for(i=startShared; i<endAlignmentIndex; i++)
    {
      
      if(alignmentMatrix[row1][i]!='N' && alignmentMatrix[row2][i]!='N') 
	{
	  sizeAligned++;
	  if(alignmentMatrix[row1][i]==alignmentMatrix[row2][i])
	    {
	      numMatches++;
	    }

	}

    }

  
  return numMatches;

}



string ReadCluster::assembleContig(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, double &nucCtr, double &badColCtr, int cutoffMinNuc, double &Nctr, double &startContig, double &longestDistN,
				   int &subClusterCtr, vector<int> &problematicCols)
{
  
  bool flag;
  int numCol, i, j, k, sizeReadsToUse, indexMax, posN;
  vector<int> nuc(4, 0); //A T G C
  string combinedNuc="";
  string allNucs="ATGC";
  int Ndistance=0;

  //int subClusterCtr=0; //this is a count of how often the second highest nucleotide count is higher than the mininum necassary to make a cluster

  vector<int> potentialColsSubCluster; //this is a vector of vectors this records rows that possibly are in the same sub cluster 

  string validChar = "ACGTacgt";
  
  flag=false;

  Nctr=0; //counter to keep track of number of Ns in the contig getting assembled

  badColCtr=0;

  numCol=alignmentMatrix[0].size(); //all rows are the same size
  sizeReadsToUse=rowsToAssemble.size();


  spp::sparse_hash_map<uint_fast64_t, double> mismatchedCount;

  spp::sparse_hash_map<uint_fast64_t, double> nucleotideCtr; //two hash tables to keep track fo which rows consitently mismatch the consensus
  // if a row mismatches the consensus at a N then you start counting 
  //the number of mismatches to the consesus for that row this count is stored in mismatchedCount the total nubmer of nucleotides 
  //for that row is stored in nucleotideCtr if the number of mismatches every gets higher than 25% of the total number of nucleotides
  //than start ignoring that row in assembly 

 
  
  double cleanDistance=0;
  double maxCleanDistance=0;

  int currentRow=0; //index into the current rows to assemble vector

  for(i=0; i<numCol; i++)
    {
      
      nuc[0]=0;
      nuc[1]=0;
      nuc[2]=0;
      nuc[3]=0;

      nucCtr=0;
      
      //for every column run through all the rows that will be used in this assembly
      for(k=0; k<sizeReadsToUse; k++)
	{
	  j=rowsToAssemble[k]; //j is row index to use in assembly
	   
	  
	  if(nucleotideCtr.size()>0)
	    {
	  //if the row has to many mismatches relative to the consensus than do not consider this read any longer
	      if(nucleotideCtr.count(j)>0)
		{
		  if(nucleotideCtr[j]==-1)
		    {
		      continue;
		    }
		}
	    }


	  switch(alignmentMatrix[j][i])
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

      


     

      //indexMax=-1;
      //indexMax=distance(nuc.begin(), std::max_element(nuc.begin(), nuc.end()));
       
      
      //sort the nuc vector but retain the index of the sorted vector
      vector<size_t> idx(nuc.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
	   [&nuc](size_t i1, size_t i2) {return nuc[i1] < nuc[i2];});

      indexMax=idx[idx.size()-1];


      

      
      if(Nctr>0)
	{

	  //check to see if those rows which are suspect are mismatching the consensus
	  for(auto iter=nucleotideCtr.begin(); iter!=nucleotideCtr.end(); iter++)
	    {
	      //if read has already been considered bad than do not consider it again just move on
	      if(iter->second==-1)
		{
		  continue;
		}


	      if(alignmentMatrix[iter->first][i]!=allNucs[indexMax] && alignmentMatrix[iter->first][i]!='N')
		{
		  mismatchedCount[iter->first]=mismatchedCount[iter->first]+1;
		  nucleotideCtr[iter->first]=nucleotideCtr[iter->first]+1;
	
		  
		  //if the suspect rows ever get a mismatch relative to the consensus of greater than 25% than do not use those rows in the future
		  if((nucleotideCtr[iter->first]>5) &&  ((mismatchedCount[iter->first]/nucleotideCtr[iter->first]) > 0.25) )
		    {
		      
		      iter->second=-1;
		      
		    }
		  
		}

	    }
      

	 
	      
	 
	}



      if(nuc[indexMax] >= cutoffMinNuc && ((nuc[indexMax]/nucCtr) >= 0.75)) //filtering clusters
	{
	  //if the second highest nucleotide counter is higher than the minimum necassary to call a base
	  //then increment the subClustercounter
	  if(nuc[idx[idx.size()-2]] >= 2*cutoffMinNuc)
	    {
	      subClusterCtr++;
	      problematicCols.push_back(i);
	    }



	  if(flag==false)
	    {
	      startContig=i;
	    }

	  flag=true;	  
	  //cout<<nuc[indexMax]<<endl;	  

	  combinedNuc+=allNucs[indexMax];

	  cleanDistance++;

	  //keeping track of how many bases in a row in the contig are called cleanly 
	  if(cleanDistance>maxCleanDistance)
	    {
	      maxCleanDistance=cleanDistance;
	    }

	}else
	{

	  //if have already started the contig but not confident in next base call then place an N. 
	  if(flag)
	    {
	      //if you have a long column of nucleotides that don't agree than increment the N counter 
	      if(nucCtr>=cutoffMinNuc)
		{
	      
		  Nctr++;
		  		  
		  cleanDistance=0;
		  

		  //check to see which rows mismatched the consensus if they consistentally mismatch the consensus ignore those rows in the future
		  for(k=0; k<sizeReadsToUse; k++)
		    {
		      j=rowsToAssemble[k]; //j is row index to use in assembly
	   
		      //if a read has already been considered and discarded than move on 
		      if(nucleotideCtr.count(j)>0)
			{
			  if(nucleotideCtr[j]==-1)
			    {
			      continue;
			    }
			}

		      //does not match the consesus
		      if(alignmentMatrix[j][i]!=allNucs[indexMax] && alignmentMatrix[j][i]!='N')
			{
			  mismatchedCount[j]=1;
			  nucleotideCtr[j]=1;
			    
			}

		    }
		  


		}


	     
	      //cout<<"combinedNuc is "<<combinedNuc<<endl;
	      combinedNuc+='N';
		 
	    }

	}


      

      /*
      //debugging code REMEMBER to comment out!!
      if(flag==false)
	{
	  combinedNuc+='N';
	}
      
      */
 

/*
      if(nuc[indexMax]>cutoffMinNuc)
	{	
	  combinedNuc+=allNucs[indexMax];	  
	}

*/



      if(((nuc[indexMax]/nucCtr) < 0.75))
	{
	  badColCtr++;
	}
      
    


    }
  
  
  /*uncomment when running the code for real 
   */


  posN=combinedNuc.find_last_not_of( 'N' ) +1;
  
  // Nctr=Nctr-(combinedNuc.length()-posN);


//remove trailing Ns
  combinedNuc.erase( posN ); 
  
  longestDistN=maxCleanDistance;

  return(combinedNuc);

}


void ReadCluster::getDistanceGraph(std::vector<std::vector<char>> &alignmentMatrix, vector<int> &startMatrix,  Graph &numMatches, int readSize)
{
  int numReads=alignmentMatrix.size();

  int z, y;


  //cerr<<"printing the distances "<<endl;
	  for(z=0; z<numReads; z++)
	    {

	      
	      for(y=(z+1); y<numReads; y++)
		{
		  
		  int sizeAligned=0;

		  

		  int numMatch=numberDiff(alignmentMatrix, z, y, startMatrix[z], startMatrix[y], sizeAligned, readSize);

	
		  if((double(numMatch)/double(sizeAligned)) > 0.95 && sizeAligned > (readSize*0.9) ) //at least 95 percent of the aligned Reads must match with an alignment length of at least 20
		    {
		      //cerr<<numMatch<<" ";

		      //numMatches[z][y]=numMatch;
		      
		      //add new edges to the graph linking the two reads
		      numMatches.addEdge(z, y);
		      numMatches.addEdge(y, z);


		    }else if((double(numMatch)/double(sizeAligned)) < 0.5 && sizeAligned>20) //if 2 reads align but they do not match then set the number of matches between them to a negative number
		    {
		      //numMatches[z][y]=-1*numMatch;
		    }else
		    {
		      //numMatches[z][y]=0;
		    }
		    
		  
		}

	      //cerr<<endl;
	    }

	  
}

void ReadCluster::getSubClusters(vector<vector<int>> &numMatches, vector<vector<int>> &clusters)
{
  
  int numReads=numMatches.size();
  int z, y, clusterCtr;
  bool addedARead; //a variable to record if having gone through every read in the aligment matrix at least one of them was added to a cluster if I never add a read to a cluster then stop the algorithm


  int notInClusterFlag;

  dense_hash_map<int_fast64_t, int_fast64_t, customHash> inCluster; //hash table the key is a row index the value is which cluster the read corresponding to that row is found in. The value is an index into the clusters vector
  inCluster.set_empty_key(-1);

  dense_hash_map<int_fast64_t, int_fast64_t, customHash> needsCluster; //hash table the key is a row index the value does not matter just need to record which rows are not in a cluster but also do not match any existing clusters will set the value to 1
  needsCluster.set_empty_key(-1);
  needsCluster.set_deleted_key(-2);



  clusterCtr=0;

  for(z=0; z<numReads; z++)
    {
      /*
      if(inCluster.count(z)>0) //if read is already assigned to a cluster skip it
	{
	  continue;
	}
      */


      for(y=(z+1); y<numReads; y++)
	{

	  //see getDistanceFunction I set a cutoff there if percent of matches for the alignment is 0.95 or greater and the size of the Aligned region is 20 then number of Matches is positive otherwise
	  //it is set to 0 or a negative number
	  if(numMatches[z][y]>0) 
	    {
	      if(clusters.size()==0) //add a new cluster if there are no clusters
		{
		  vector<int> temp;

		  temp.push_back(z);
		  temp.push_back(y);

		  clusters.push_back(temp);
		  
		  inCluster[z]=clusterCtr;
		  inCluster[y]=clusterCtr;
		  
		  clusterCtr++;
		  continue;
		}

	      if(needsCluster.count(z)>0 || needsCluster.count(y)>0) //if read 1 or read 2 needs a cluster
		{
		  vector<int> temp;

		  temp.push_back(z);
		  temp.push_back(y);

		  clusters.push_back(temp);
		  
		  inCluster[z]=clusterCtr;
		  inCluster[y]=clusterCtr;
 
		  clusterCtr++;
		  
		  needsCluster.erase(z);
		  needsCluster.erase(y);

		  if(needsCluster.count(z)>0)
		    {
		      cerr<<"Mistake!! element should not exist "<<endl;
		    }

		  continue;
		} 

	      

	      if(inCluster.count(z)==0 && inCluster.count(y)>0) //if second read is in cluster but first read isn't add the first read to the second reads cluster
		{
		  clusters[inCluster[y]].push_back(z);  //value in hash table is a row into the clusters vector
		  inCluster[z]=inCluster[y]; //adding the first read to the hash table recording which cluster a read belongs to 

		}

	      if(inCluster.count(y)==0 && inCluster.count(z)>0) //if first read is in cluster but second read isn't add the second read to the first reads cluster
		{
		  clusters[inCluster[z]].push_back(y);  //value in hash table is a row into the clusters vector
		  inCluster[y]=inCluster[z]; //adding the 2nd read to the hash table recording which cluster a read belongs to 

		}

	    }


	  if(numMatches[z][y] < -1) //if number of matches is less than -1 than means reads overlap well but do not align well with many mismatches
	    {
	      if(inCluster.count(z)==0 && inCluster.count(y)>0) //first read is not in a cluster second read is in a cluster hence first read does not belong in the same cluster as the second read 
		{
		  needsCluster[z]=1;
		}

	      if(inCluster.count(y)==0 && inCluster.count(z)>0) //first read is not in a cluster second read is in a cluster hence first read does not belong in the same cluster as the second read 
		{
		  needsCluster[y]=1;
		}


	    }


	}
    }


}

void ReadCluster::printMatrix(std::vector<std::vector<char>> &alignmentMatrix, std::vector<int> &rowsToAssemble, string &combinedNuc, double &percentBadCol, string &clusterID, std::ofstream &debugging, double &percentNs)
{
  int numCol, i, j, k;

  int sizeReadsToUse=rowsToAssemble.size();

#pragma omp critical(DEBUGGING_CLUSTER)
  {  
  //printing the matrix
  //raw matrix before filtering
  
    numCol=alignmentMatrix[0].size(); //all rows are the same size  

    for(k=0; k<sizeReadsToUse; k++)
      {
	i=rowsToAssemble[k];

	for(j=0; j<numCol; j++)
	  {
	  //matrix[i][((numCol/2)-startPositions[i])+j]=readSeq[i][j];
	
	    debugging<<alignmentMatrix[i][j];
	  
	  }
	debugging<<endl;
      }

    
    debugging<<"cluster being assembled is  "<<clusterID<<endl;
    debugging<<"contig is "<<combinedNuc<<endl;
    debugging<<"percentage of bad columns is "<<percentBadCol<<endl;
    debugging<<"percentage of Ns is "<<percentNs<<endl;
    debugging<<"finished printing out the matrix "<<endl;
    debugging<<"##############################################\n\n\n\n"<<endl;
    

  }




}



//give a vector of reads obtain the position of every kmer in every read. 
//Do not need to worry about reverse complement all reads should be on the 
//same strand enforced this earlier in code
uint_fast64_t ReadCluster::getKmers()
{

  uint_fast64_t kmer, positionCtr, AbitString, CbitString, GbitString, TbitString; 
  long  i, j, k;

  kmerPositions.clear();
  kmerInReadCounts.clear();
  presentMultipleTimes.clear();

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




 //return(1);

 
 //dense_hash_map<uint_fast64_t, long, customHash> kmerInReadCounts; //hash table contains a count of how many reads contain a given kmer key is the kmer value is the count 
 //kmerInReadCounts.set_empty_key(-20);
 //kmerInReadCounts.resize(1000000);

 //return(1);

   //revKey=0;
  kmer=0;
  positionCtr=0;

  
  //than is not an appropriate kmer to use in assembly 
  //presentMultipleTimes.set_empty_key(false);


  //iterator through ever read in the cluster getting a list of the kmers 
  for(j=0; j<readSeqs.size(); j++)  
  {

    //kmerInReadCounts[j]=1;
    
    //    if(j%1000==0)
    //{
    //cerr<<"on Read "<<j<<endl;
	//}

    //continue;


    dense_hash_map<uint_fast64_t, long, customHash> readKmers; //hash table contains the kmers and their position in the read 
    readKmers.set_empty_key(-20);

    //dense_hash_map<uint_fast64_t, long, customHash> alreadyPresent; //hash table keeps track of kmers that have already been found earlier in the read key is the kmer the value is a flag set to 1 indicating the 
    //kmer is present in the read
    //alreadyPresent.set_empty_key(-20);


    //     std::vector< google::dense_hash_map<uint_fast64_t, long, customHash> > kmerPositions;

   

    kmer=0;
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

	 
	      /*
#pragma omp critical(DEBUGGING_READKMERS)
	      {
	      
		cout<<"kmer is  "<<kmer<<"\t kmer translated from bit string is "<<bit2String(kmer, 32)<<" bit representation is "<<std::bitset<64>(kmer)<<endl;
	      
		cout<<"hashed value for key is "<<SpookyHash::Hash64(&kmer, 8, 20)<<endl;
	      }
	      */
	      

	      //lookup key if key not found will return iterator that points to end of hash table
	      auto iter=readKmers.find(kmer);
	      
	      if(iter!=readKmers.end())
		{
		  presentMultipleTimes[kmer]=true;

		}

	      
	      if(iter==readKmers.end()) //means kmer is not present earlier in the read this is important because only want to count a kmer once per read
		{
		  kmerInReadCounts[kmer]=kmerInReadCounts[kmer]+1;		  
		  readKmers[kmer]=i;	 
		  
		  //alreadyPresent[kmer]=1; //setting value to 1 indicating kmer has been found in this read
		}



	      
	      //continue;
	     
	      
     
	   
	      
	}

      //cout<<"size of hash table in this function is "<<readKmers.size()<<endl;
      kmerPositions.push_back(readKmers);


  }


  
  signed long max=-1;
  uint_fast64_t maxKmer=0;

  

  //going through hash table which counts the number of times a kmer is found in the population of reads obtaining the kmer that occurs most often
  auto iterHash=kmerInReadCounts.begin();
  for(iterHash; iterHash!=kmerInReadCounts.end(); iterHash++)
    {

      if(iterHash->second > max)
	{

	  
	  
	   //checking to see if a potential max kmer is present multiple times in at least one read if it is then do not use that kmer as the max kmer
	  if(presentMultipleTimes.count(iterHash->first)>0)
	    {
	      //cerr<<"read "<<j<<" "<<maxKmer<<endl;
	      
	      continue;
	    }
	  
	  

	  maxKmer=iterHash->first;
	  max=iterHash->second;
	}
      //pair<uint_fast64_t, long> keyPair(iterHash->first, iterHash->second);

      //kmerCounts.push_back(keyPair);
    }


  
  //sort(kmerCounts.begin(), kmerCounts.end(), comparePairs);

 

  /*
  cout<<"paired kmers are "<<endl;

  for(k=0;k<kmerCounts.size();k++)
        cout << kmerCounts.at(k).second << " (" << kmerCounts.at(k).first << "%)"<< endl; 

  cout<<"#############################\n\n\n\n\n";
  
  */

  return(maxKmer);

}

//before calling this function make sure getKmers and setStartPositions has already been called
void ReadCluster::revCompCluster(uint_fast64_t maxKmer)
{

  
  uint_fast64_t key, revKey;
 
  
  //building the look up table to reverse nucleotide bit strings
   int bitWord=16;
   std::vector<uint_fast32_t> bitTable(pow(2,bitWord));

   //createBitRevTableMachineWord(bitWord, bitTable);

   createBitRevTableMachineWord(bitWord, bitTable);


   revKey=0;
   revComplementMachineWord(maxKmer, revKey, clusterKmerSize, bitTable);
    
   vector<bool> flipped(sameStrand.size(), false); //this is to keep track of those reads that have been reversed complemented if a read has been reverse complemented and another read matches it 
   //then you know the other read should also be reverse complemented
   


   //iterating through the set of reads asking if the reverse complement of the max kmer is present in any reads if it is then reverse complement the read and its quality
  auto iter=kmerPositions.begin(); //vectorIterator
  long ctr=0;
    for(iter; iter!=kmerPositions.end(); iter++)
      {
	
	auto Hashiter=(*iter).find(revKey); //look for the kmer in the hash table that contains every kmer for the given read
	if(Hashiter!=(*iter).end()) //kmer present in read
	  {
	    if(sameStrand[ctr]==false)
	      {
		sameStrand[ctr]=true; //if have not already assembled this read into the cluster
		revComplement(readSeqs[ctr]); //reverse complement and replace the orginal sequence
		std::reverse(qualityStrings[ctr].begin(), qualityStrings[ctr].end()); //reverse the quality scores
		flipped[ctr]=true;

	      }
	  }

	ctr++;
     }


    /*    
    int x;
    for(x=0; x!=sameStrand.size(); x++)
      {

	cerr<<sameStrand[x]<<endl;
       
      }
    

    */

    /*   
   int y;
   for(y=0; y!=readSeqs.size(); y++)
     {

       cerr<<readSeqs[y]<<endl;

     }
    */



    ctr=0;
    bool foundFlag=false;
    int i, j;
    uint_fast64_t unknownKmer, unknownRevKmer;
    uint_fast32_t ctrUnassigned=1;
    uint_fast32_t ctrUnassignedPrev=1000;


    //if all the reads 
    while(ctrUnassigned!=0 && (ctrUnassignedPrev!=ctrUnassigned))
      {   
	ctrUnassigned=ctrUnassignedPrev;

	for(i=0; i!=sameStrand.size(); ++i)
	  {
	    ctrUnassigned=0;

	    if(sameStrand[i]==false)
	      {
		foundFlag=false;

	    //look at every kmer on the read for which the strand is not known try and find a match for reads on which the strand is known
		auto HashIter=kmerPositions[i].begin();
		for(HashIter; HashIter!=kmerPositions[i].end(); HashIter++)
		  {
		    if(foundFlag) //if have already found a matching kmer break out of loop
		      {
			break;
		      }
	      
		    unknownKmer=HashIter->first;


		    for(j=0; j!=kmerPositions.size(); ++j) //go through every read looking for ones that the strand is known 
		      {
			if(sameStrand[j]==true) //if a read has been assigned a strand check to see if it has a kmer in common with the unknown strand read
			  {
			    auto findIter=kmerPositions[j].find(unknownKmer); //look for kmer in hash table
		      
			    if(findIter!=kmerPositions[j].end()) //kmer is present on the same strand
			      {
				sameStrand[i]=true;
				foundFlag=true;
				//cerr<<"on same strand for unkonwn read "<<i<<" it matches read "<<j<<endl;

			      //if the read that there is a match to has itself been reverse complemented than rev complement the unknown read to match 
			      if(flipped[j]==true)
				{
				  //				  cerr<<"on same strand for unkonwn read  rev complemtned "<<i<<" it matches read "<<j<<endl;
				      

				  revComplement(readSeqs[i]); //reverse complement and replace the orginal sequence
				  std::reverse(qualityStrings[i].begin(), qualityStrings[i].end()); //reverse the quality scores
				  
				  flipped[i]=true;

				}

				break;
			      }

			    unknownRevKmer=0;
			    revComplementMachineWord(unknownKmer, unknownRevKmer, clusterKmerSize, bitTable);
			    auto findIterRev=kmerPositions[j].find(unknownRevKmer); //look for kmer in hash table
		      
			    if(findIterRev!=kmerPositions[j].end()) //kmer is present on the opposite strand
			      {
				sameStrand[i]=true;
			
				if(flipped[j]==false)
				  {
				    revComplement(readSeqs[i]); //reverse complement and replace the orginal sequence
				    std::reverse(qualityStrings[i].begin(), qualityStrings[i].end()); //reverse the quality scores
			
				    //  cerr<<"on opposite strand for unkonwn read rev complemented "<<i<<" it matches read "<<j<<endl;
			
	    
				    flipped[i]=true;
				  }
				
				foundFlag=true;
				//cerr<<"on opposite strand for unkonwn read "<<i<<" it matches read "<<j<<endl;
			      
				break;
			      
			      }


			  }

		      }

		  }

		if(foundFlag==false)
		  {
		    ctrUnassigned;
		  }

	      }

	
	    
	  }

      }

    /*
    int z;
   for(z=0; z!=sameStrand.size(); z++)
     {

       cerr<<sameStrand[z]<<endl;

     }
    */

    /*          
   int z;
   for(z=0; z!=readSeqs.size(); z++)
     {

       cerr<<readSeqs[z]<<endl;

     }
    
    */


}




//very simple compare function to pass to sort in order to sort pairs by the second element
bool comparePairs(const std::pair<uint_fast64_t, long>&i, const std::pair<uint_fast64_t, long>&j)
{
    return i.second > j.second;
}



/*
std::pair<uint_fast64_t, long> ReadCluster::getPair(long index) //given an index return the pair kmer count pair corresponding to that index
{

  return(kmerCounts[index]);

}
*/


void ReadCluster::setKmerSize(int size)
{

  clusterKmerSize=size;

}

void ReadCluster::setStartPositions(uint_fast64_t kmer)
{
  sameStrand.clear();
  std::fill(startPositions.begin(), startPositions.end(), -1);
  
  auto iter=kmerPositions.begin(); //vectorIterator
  long ctr=0;
    for(iter; iter!=kmerPositions.end(); iter++)
      {
	
	auto Hashiter=(*iter).find(kmer); //look for the kmer in the hash table that contains every kmer for the given read
	if(Hashiter!=(*iter).end()) //kmer present in read
	  {
	    
	    if(startPositions[ctr]==-1) //if have not already assembled this read into the cluster
	      {
		startPositions[ctr]=(*iter)[kmer];
		sameStrand.push_back(true);
		//		cerr<<"inside set StartPositions function "<<endl;
	      }
	  }else
	  {
	    sameStrand.push_back(false);
	  }

	ctr++;
      }


    
    /*
    cerr<<"starting to print same Staand stuff size is"<<sameStrand.size()<<endl;
    auto iterStrand=sameStrand.begin();

    for(iterStrand; iterStrand!=sameStrand.end(); iterStrand++)
      {
	cerr<<(*iterStrand)<<endl;


      }
    */

    //cerr<<"sameSTrand size is "<<sameStrand.size()<<endl;
}

void ReadCluster::printStartPositions()
{

  
  auto iter=startPositions.begin();

    for(iter; iter!=startPositions.end(); iter++)
      {
	cerr<<(*iter)<<endl;


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

 
void Graph::addEdge(int v, int w)
{
    graph[v].push_back(w); // Add w to v’s list.
}
 
void Graph::DFSUtil(int v, std::vector<bool> &visited, std::vector<int> &cluster)
{
    // Mark the current node as visited and store it
  visited[v] = true;
    //cerr << v << " ";
 
  cluster.push_back(v);

    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for (i = graph[v].begin(); i != graph[v].end(); ++i)
        if (!visited[*i])
	  DFSUtil(*i, visited, cluster);
}
 
// DFS traversal of the vertices reachable from v. It uses recursive DFSUtil()
void Graph::DFS(int v, std::vector<int> &cluster, std::vector<bool> &visited)
{
    // Mark all the vertices as not visited
    //bool *visited = new bool[V];
    
  // vector<bool> visited(graph.size());

  //for (int i = 0; i < graph.size(); i++)
  //visited[i] = false;
 
  //  cerr<<"started the recursive helper function "<<endl;
    // Call the recursive helper function to print DFS traversal
  DFSUtil(v, visited, cluster);
}
 



/*
void ReadCluster::setUsedReadsFalse()
{
  if(clusterKmerSize<1)
    {
      cerr<<"Error number of reads in cluster is less than 1 "<<endl;
      exit(0);

    }else
    {

      usedReads.assign(clusterKmerSize, false);
      
    }



}
*/
