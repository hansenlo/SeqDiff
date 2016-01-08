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

  vector<string> test;

  test.push_back("/data/Temp/0.dat");

  //test is the file containing the reads in a contig
  //4 is the cutoff number of reads 
  //20 is the kmer size used to assemble the reads into a contig
    readInClusters(test, 4, 20);


  //assembleContigs();

  return(0);


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
    getReads(uniqueKmers, 256, continueFlag, expSeqLib, kmerSize);


    //test.push_back("/data/Temp/0.dat");

    cerr<<"Starting to assemble clusters "<<endl;

  //test is the file containing the reads in a contig
  //4 is the cutoff number of reads 
  //20 is the kmer size used to assemble the reads into a contig
    //readInClusters(test, 4, 20);


    return(0);
   

}
