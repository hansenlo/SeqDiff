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

#include <unordered_map>

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
  masterKey.set_empty_key(-20);




  
  
   //sparse_hash_map<uint_fast64_t, int, customHash> uniqueKmers;

  //dense_hash_map<uint_fast64_t, int, customHash> uniqueKmers;

  bitset<bitSetSize> test;

   //test=1*pow(2, (2*kmerSize)-2)+1*pow(2, ((2*kmerSize)-1));

  test.set(test.size()-1);

  //sparse_hash_map<bitset<192>, int, stdHash> temp;


  //std::unordered_map<uint_fast64_t, int, customHash> temp;


  //std::unordered_map<bitset<bitSetSize>, int, stdHash> temp;

  //  sparse_hash_map<bitset<bitSetSize>, int, stdHash> temp;

  //dense_hash_map<bitset<bitSetSize>, int, stdHash> uniqueKmers;

  spp::sparse_hash_map<bitset<bitSetSize>, int, stdHash> uniqueKmers;
  
  

  //size_t hashvalue=boost::hash_value(test.m_bits);

   
  //hash<bitset<30>> foo;

  //size_t hashvalue=foo(test);

							  


   //temp.set_empty_key(NULL);

   
   //SpookyHash::Hash64(&value, 8, 20)

   


  bitset<30> foo, temp, foo2;
									 
   foo.set(foo.size()-1);


   //cout<<"size of bitset object is "<<test.size()<<" number of bytes is "<<sizeof(test)<<" bit representation is "<<test<<"equality test is "<<(test==foo)<<endl;
   
   //   return(0);





   //dense_hash_map<bit, int, customHash> uniqueKmers;



															  //return(0);


   sparse_hash_map<uint_fast64_t, int, customHash> controlKmers;

   //dense_hash_map<uint_fast64_t, int, customHash> uniqueKmers;
  //dense_hash_map<uint_fast64_t, ReadClusters, customHash> clusters;  

   //cerr<<"reached this point before set empty key"<<endl;

   bitset<bitSetSize> emptyKey;

   //################NEED TO ADD CODE MAKING SURE kmer size is small enough this never happesn##############
   emptyKey.set(bitSetSize-1);
   
  //clusters.set_empty_key(-20);
  //uniqueKmers.set_empty_key(emptyKey);
  //uniqueKmers.set_deleted_key(-2);
  //uniqueKmers.resize(100000000);
  
   uniqueKmers.reserve(450000000);
//  controlKmers.resize(10000000000);
  

  //cerr<<"reached this point after set empty key"<<endl;

  

  //whether or not sequence on continous lines should be connected
  // a 1 indicates should be connected a 0 indicates should not be connected
  //continueFlag=*argv[4];
  continueFlag='0';


  if(kmerSize==0)
    {
      cerr<<"Should never get a kmer size of zero"<<endl<<"conversion from string to int must not have worked"<<endl;
      return(1);
    }

  cerr<<"starting to read in unique kmers "<<endl;

  //8 is the kmer count cutoff
  readUniqueKmers(uniqueKmers, continueFlag, uniqueExpKmerCountFile, kmerSize, 8); //Need to uncomment for code to work

  
  //return(0);
														

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
  fileNames=getReads(uniqueKmers, 256, continueFlag, expSeqLib, kmerSize); //NEED TO UNCOMMENT
   
    uniqueKmers.clear(); //reallocating all the memory held by the hash table of unique kmers

    //mergeClusters(fileNames); //Does not work!! 



    int currentClusterFilePublic=0; //index of the current cluster file that has not yet been assembled ctr will be shared by all threads
    int currentClusterFilePrivate=0; //index of the cluster file the thread is currently working on will be private

    //test.push_back("/data/Temp/0.dat");

    cerr<<"Starting to assemble clusters "<<endl;

    cerr<<"Number of files is "<<fileNames.size()<<endl;

    int thread_count=35;
    int tid;


    
    ofstream contigOut("/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/contigs.fa"); //location of where to put file containing contigs
    long clusterID=0;

    ofstream debuggingMatrix("/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/debuggingMatrix.dat"); //location of where to put file containing the matrix of aligned reads used to assemble contigs

    
    /*
    vector<string> testAssembly;

    testAssembly.push_back("/data2/Temp/0.dat");

  //test is the file containing the reads in a contig
  //4 is the cutoff number of reads 
  //20 is the kmer size used to assemble the reads into a contig
  //readInCluster(test[0], 4, 20);

    cerr<<"starting to test cluster assembly "<<endl;
    readInCluster(testAssembly[0], 4, 25, tid, contigOut, clusterID, debuggingMatrix);
  //assembleContigs();

  return(0);
    */




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

	  //currentClusterFilePrivate=1; //MUST REMEBER TO DELETE FOR DEBUGGING PURPOSES ONLY
	
	}

	//#pragma omp critical(TEST)
	//{
	  //cerr<<"currentClusterFilePrivate index is "<<currentClusterFilePrivate<<endl;


	//3 is the cutoff number of reads 
	//20 is the kmer size used to assemble the reads into a contig
	
	readInCluster(fileNames[currentClusterFilePrivate], 4, 25, tid, contigOut, clusterID, debuggingMatrix);

	  // readInCluster(fileNames[1], 4, 20, tid);
	 
	   //currentClusterFilePrivate++;
  
	   //cerr<<"after readINCluster function call "<<endl;

	  //}

    
      }



  }

  unordered_map<std::string, std::string> genome;


  /*
  cerr<<"calling zygozity "<<endl;

  readInFasta(genome, "test.fa");

  cerr<<"chr 1 is "<<endl;
  cerr<<genome["chr1"]<<endl;


  cerr<<"chr 2 is "<<endl;
  cerr<<genome["chr2"]<<endl;
  */


    
  	  
  contigOut.close();
  debuggingMatrix.close();

  int i;

  /*
  //deleting the cluster files
  for(i=0; i<fileNames.size(); i++)
    {

      if(remove(fileNames[i].c_str()))
	{
	  cerr<<"Temporary cluster file "<<fileNames[i]<<" Was not deleted properly please check the file/directory containing the file has proper permissions "<<endl;

	}
    }
  */

    return(0);
   

}
