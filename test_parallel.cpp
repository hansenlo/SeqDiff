#include<omp.h>
#include<iostream>
#include<sparsehash/sparse_hash_map>
#include<sparsehash/dense_hash_map>
#include<fstream>
#include<string>
#include<stdlib.h>
#include<list>
#include <omp.h>
#include <seqan/store.h>
#include <seqan/consensus.h>

using namespace seqan;

using google::dense_hash_map;
using google::sparse_hash_map;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::list;
using namespace seqan;

struct node {
  vector<string> chunk;
};


int main(int argc, char *argv[] )
{

  int thread_count=20, tid, nthreads, i, lineCtr=0, nodeCtr=0, sum=0;
  int chunkSize=10, loadedFileFlag;
  std::ifstream seqFile;
  string file="/home/hansenlo/SeqDiff/svn/repos/Project/testData.dat";
  string line;

  node *currentNode, *workNode;

  list<node*> workList;
 

  bool test=false;

  if(!test)
    {
      cout<<"man I'm good "<<endl;
      return(0);
    }else
    {
      cout<<"Clearly I am not as good as I think I am "<<endl;

    }

  seqFile.open(file.c_str());

  
  currentNode=new node;

  currentNode->chunk.push_back("test");

  cout<<currentNode->chunk[0]<<"\n";


  uint_fast64_t temp=0;


  return(0);
  

  //first load the linked list with initial chunks of work chunkSize is number of sequences in each chunk of work
  


  loadedFileFlag=0;


#pragma omp parallel num_threads(thread_count) default(shared)\
  private(tid, i, workNode, sum)
  {

    tid=omp_get_thread_num();
  
    while(loadedFileFlag==0)
      {

	//only access input file if you are the master thread
	#pragma omp master
	{
      //if list is empty get some work
	  if(workList.empty())
	    {

	      lineCtr=0;

	      while(seqFile.good() && lineCtr < (30*chunkSize))
		{
		  if(lineCtr % chunkSize==0)
		    {
		      if(lineCtr>0)
			{
			  #pragma omp critical
			  {
			    workList.push_back(currentNode);
			  }
			}

		      currentNode=new node;
		    }


		  getline(seqFile, line);
		  lineCtr++;
		  currentNode->chunk.push_back(line);
	     
	      

		}
	    }
	    
	    
	  if(!seqFile.good())
	    {
	
	      //cout<<"last Node size is "<<currentNode->chunk.size()<<"\n";
	  //add the very last bit of file
	       #pragma omp critical
	      {
		workList.push_back(currentNode);
	      }
	      
	      loadedFileFlag=1;
	
	    }


	}

	workNode=NULL;
     
           #pragma omp critical 
	    {
	      if(!workList.empty())
		{
		  
		  workNode=workList.back(); //copy get some work
		  //workList.pop_back(); //remove the chunk of work from the queue 
		
		    workList.pop_back();

		    /*
		    cout<<"thread is ********** "<<tid<<"\n";
		    for(i=0; i<workNode->chunk.size(); i++)
		      {	  
		    
			cout<<workNode->chunk[i]<<"\n";
				
		      } 
		    */


		}
	    }  

	    
		  
		  //workList.clear();
		  
		  //cout<<"thread is ********** "<<tid<<"\n";
		  
		  /*
		  for(i=0; i<workNode->chunk.size(); i++)
		    {	  
		    
		      cout<<workNode->chunk[i]<<"\n";
		    
		    }
		  */

	    
	    if(!(workNode==NULL))
	      {
	    //cout<<"thread is ********** "<<tid<<"\n";
		sum=0;
		for(i=0; i<workNode->chunk.size(); i++)
		  {	  
		    //#pragma omp critical
		    //{
		      //cout<<workNode->chunk[i]<<"\n";
		    //}
		    sum+=atoi(workNode->chunk[i].c_str());
		    //cout<<workNode->chunk[i]<<"\n";	    
		  } 
	    

		
              #pragma omp critical
		{
	      //cout<<"thread is ********** "<<tid<<"\n";
		  cout<<"sum is "<<sum<<" for thread "<<tid<<" chunk size is "<<workNode->chunk.size()<<"\n";
		}
	    

	      }
	    
	    

		
	    
	  

	//else //do work
	     //{

	    /*
	    //check to see if queue has work in in
	    if(!workList.empty())
	      {
		//get some work from the queue
		#pragma omp critical 
		{
		  workNode=workList.back();
		  workList.pop_back();
		}

		#pragma omp critical
		{
		  cout<<"thread is **********"<<tid<<"\n";
		}

		for(i=0; i<workNode->chunk.size(); i++)
		  {	  
		    #pragma omp critical
		    {
		      cout<<workNode->chunk[i]<<"\n";
		    }
		  }
		nodeCtr++;


	      }
	    */
	//}

      
      /*
      for(list<node *>::iterator list_iter = workList.begin(); list_iter != workList.end(); list_iter++)
	{	  
	  //cout<<**list_iter.chunk[0]<<"\n";

	  cout<<"Node ***********"<<nodeCtr<<endl;
	  for(i=0; i<(*list_iter)->chunk.size(); i++)
	  {	  
	    cout<<(*list_iter)->chunk[i]<<"\n";

	  }
	  nodeCtr++;
	}      
      

      
      for(i=0; i<workList.size(); i++)
	{
	  workNode=workList.back();
	  workList.pop_back();

	}
      
      */

      }

  }
  

  /*

#  pragma omp parallel num_threads(thread_count) default(shared)\
  private(tid, i)
  {
    tid=omp_get_thread_num();

    #pragma omp critical
    {
      std::cout<<"hello world from thread "<<tid<<std::endl;
    }

    for(i=0; i<10; i++)
      {
	
	//#pragma omp atomic	
	temp++;
      }
	    
    if(tid==0)
      {
	nthreads=omp_get_num_threads();
	cout<<"number of threads from the master is "<<nthreads<<endl;


      }

  }

  cout<<"temp after many parrallel computation is "<<temp<<endl;

  */

}
