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
//example

//variantFinder <exp kmer count file> <exp sequencing lib fastq or fasta> <kmerSize>

int main(int argc, char *argv[] )
{

 if(argc!=4)
    {
      cerr<<"expecting five command line arguments exiting"<<endl;
      return(0);
    }

  
  int kmerSize, cutoff, readLength, timesRun;
  char continueFlag;
  string controlSeqLib, outputFile, expSeqLib, uniqueExpKmerCountFile;
  ifstream seqFile;
  string line;
  char start;
  bool fastq;
  vector<string> fileNames; //names of the files clusters are stored in

  //number of characters in kmer
  kmerSize=atoi(argv[3]);
  
  //Number of variants allowed in the control 
  //cutoff=atoi(argv[5]);

  //file that contains the counts of unique kmers in the experiment sequencing  library
  uniqueExpKmerCountFile=argv[1];

  //experiment sequencing library
  expSeqLib=argv[2];

  //name of outputFile that control kmers could be printed to
  //outputFile=argv[6];



  //dense_hash_map<uint_fast64_t, ReadCluster *, customHash> allClusters;

  sparse_hash_map<uint_fast64_t, ReadCluster *, customHash> allClusters;
  allClusters.resize(30000);
  
  dense_hash_map<uint_fast64_t, uint_fast64_t, customHash> masterKey; //hash table contains the first middle and end kmer of the master read that defines a cluster
  masterKey.set_empty_key(-1);

  /*
  vector<string> test;

  test.push_back("/data/Temp/0.dat");

  //test is the file containing the reads in a contig
  //4 is the cutoff number of reads 
  //20 is the kmer size used to assemble the reads into a contig
    readInCluster(test[0], 4, 20);


  //assembleContigs();

  return(0);
  */

   //sparse_hash_map<uint_fast64_t, int, customHash> uniqueKmers;

   dense_hash_map<uint_fast64_t, int, customHash> uniqueKmers;

   sparse_hash_map<uint_fast64_t, int, customHash> controlKmers;

   //dense_hash_map<uint_fast64_t, int, customHash> uniqueKmers;
  //dense_hash_map<uint_fast64_t, ReadClusters, customHash> clusters;  

  //clusters.set_empty_key(-1);
  uniqueKmers.set_empty_key(-1);
  uniqueKmers.set_deleted_key(-2);
  uniqueKmers.resize(100000000);
  //  controlKmers.resize(10000000000);
  

  //whether or not sequence on continous lines should be connected
  // a 1 indicates should be connected a 0 indicates should not be connected
  //continueFlag=*argv[4];
  continueFlag='0';


  if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      return(1);
    }

  //4 is the kmer count cutoff
    readUniqueKmers(uniqueKmers, continueFlag, uniqueExpKmerCountFile, kmerSize, 5);


  /*   
   //fout<<"size of hash table is "<<sizeof(allKmers)<<endl;

   sparse_hash_map<uint_fast64_t, int, customHash>::iterator myHashIteratorTemp;
  
   //unordered_map<uint_fast64_t, char>::iterator myHashIterator;
  
   for(myHashIteratorTemp=controlKmers.begin(); myHashIteratorTemp!=controlKmers.end(); myHashIteratorTemp++)
     {

       if(myHashIteratorTemp->second > 1)
       {

	   //foutTemp<<bit2String(myHashIteratorTemp->first, kmerSize)<<"\t"<<myHashIteratorTemp->second<<endl;
       
	   cout<<bit2String(myHashIteratorTemp->first, kmerSize)<<"\t"<<myHashIteratorTemp->second<<endl;
       
	   

	   }

       
     }
  */

    cerr<<"finished getting control kmers"<<endl;

      //256 is number of files to split clusters into
    fileNames=getReads(uniqueKmers, 256, continueFlag, expSeqLib, kmerSize);


    int currentClusterFilePublic=0; //index of the current cluster file that has not yet been assembled ctr will be shared by all threads
    int currentClusterFilePrivate=0; //index of the cluster file the thread is currently working on will be private

    //test.push_back("/data/Temp/0.dat");

    cerr<<"Starting to assemble clusters "<<endl;

    int thread_count=20;
    int tid;


    
    ofstream contigOut("/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/contigs.fa"); //location of where to put file containing contigs
    long clusterID=0;

    ofstream debuggingMatrix("/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/debuggingMatrix.dat"); //location of where to put file containing the matrix of aligned reads used to assemble contigs

    
    //parallizing the cluster assembly 
#pragma omp parallel num_threads(thread_count) default(shared)	\
  private(tid, currentClusterFilePrivate)
  {

    tid=omp_get_thread_num();

    while(currentClusterFilePublic<fileNames.size())
      {

    //get the current cluster file that needs to be assembled 
#pragma omp critical(GET_CLUSTER)
	{
	  currentClusterFilePrivate=currentClusterFilePublic;
	  currentClusterFilePublic++;
	
	}

	//#pragma omp critical(TEST)
	//{
	  //cerr<<"currentClusterFilePrivate index is "<<currentClusterFilePrivate<<endl;


	//4 is the cutoff number of reads 
	//20 is the kmer size used to assemble the reads into a contig
	
	readInCluster(fileNames[currentClusterFilePrivate], 4, 20, tid, contigOut, clusterID, debuggingMatrix);

	  // readInCluster(fileNames[1], 4, 20, tid);
	 
	   //currentClusterFilePrivate++;
  
	   //cerr<<"after readINCluster function call "<<endl;

	  //}
    
      }



  }
    

    /*
    int i;

    for(i=0; i<fileNames.size(); i++)
      {

	readInCluster(fileNames[i], 4, 20, tid);
      }
    */
  	  
  contigOut.close();
  debuggingMatrix.close();

    return(0);
   

}
