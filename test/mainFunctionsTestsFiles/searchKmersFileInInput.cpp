/*#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <string>
#include <vector>*/
#include "searchKmersFileInInput.hpp"

using namespace std;

//#define PRINT_DEBUG

/*
This program tests the number of kmers of a fasta file that are found in an other file.
parameters :
- input fasta file
- file to find the kmers
- kmer size

Originally, the file to find kmers was a file containing solid kmers found by leon kmer counting (executed during graph building).
It may also be an other file for other uses
*/


struct KmersStats{

	int nbKmers;
	int nbSolidKmers;

};


void searchKmers( vector<string>* kmers, string searchFileName, struct KmersStats *ks)
{

	//cout << "iterate" << endl;
	for (std::vector<string>::iterator it = kmers->begin() ; it != kmers->end(); ++it)
	{
    	string kmer = *it;
    	//just grep kmer 
    	string commandGrep1 = "grep " + kmer + " " + searchFileName;
    	//just grep rev comp of kmer
    	string commandGrep2 = "grep " + getReverseComplement(kmer) + " " + searchFileName;
    	//grep kmer and rev comp
    	//string commandGrep3 = "grep " + kmer + "\\|" + getReverseComplement(kmer) + " " + searchFileName;
    	string commandGrepLines1 = commandGrep1 + " | wc -l" ;
    	string commandGrepLines2 = commandGrep2 + " | wc -l" ;
    	//string commandGrepLines3 = commandGrep3 + " | wc -l" ;
    	FILE* cmd1;
    	FILE* cmd2;
    	//FILE* cmd3;
    	char result1[1024];
    	char result2[1024];

    	cmd1 = popen(commandGrepLines1.c_str(), "r");

    	if (cmd1 == NULL) {
        	perror("popen");
        	exit(EXIT_FAILURE);
   	 	}

   		cmd2 = popen(commandGrepLines2.c_str(), "r");

    	if (cmd2 == NULL) {
        	perror("popen");
        	exit(EXIT_FAILURE);
   	 	}

   	 	while (fgets(result1, sizeof(result1), cmd1) && fgets(result2, sizeof(result2), cmd2)) 
   	 	{
   	 		ks->nbKmers += 1;
   	 		#ifdef PRINT_DEBUG
        	printf("%s", result1);
        	printf("%s", result2);
        	#endif

        	if (atoi(result1) || atoi(result2))
        	{
        		#ifdef PRINT_DEBUG
        		cout << "res1 && res2" << endl;
        		#endif
        		ks->nbSolidKmers += 1;
        	}
        	#ifdef PRINT_DEBUG
        	cout << "nb kmers : " << ks->nbKmers << endl;
        	cout << "nb solid kmers : " << ks->nbSolidKmers << endl;
        	#endif
   		}
/*
   	 	while (fgets(result, sizeof(result), cmd2)) 
   	 	{
        	printf("%s", result);
   		}

   		cmd3 = popen(commandGrepLines3.c_str(), "r");
*/
    	/*if (cmd3 == NULL) {
        	perror("popen");
        	exit(EXIT_FAILURE);
   	 	}

   	 	while (fgets(result, sizeof(result), cmd3)) 
   	 	{
        	printf("%s", result);
   		}*/

        #ifdef PRINT_DEBUG
   		printf("\n\n");
   		#endif

    	pclose(cmd1);
    	pclose(cmd2);
    	//pclose(cmd3);
    }
}

int main(int argc, char const *argv[])
{

	if (argc != 4){
		cerr << "usage : arg1 fasta input; \narg2 file to search kmers in;\narg3 kmer size" << endl;
		exit(EXIT_FAILURE);
	}
 
	string input = argv[1];
	string searchFileName = argv[2];
	int kmerSize = atoi(argv[3]);

 	ifstream fastaFile;
  fastaFile.open (input);
  //ofstream searchFile;
  //searchFile.open(searchFileName);
  struct KmersStats *ks = (struct KmersStats *) malloc(sizeof(struct KmersStats));
 	ks->nbKmers = 0;
  ks->nbSolidKmers = 0;


 	string read;
 	//read two time to get the read, and skip headers
 	while(getline (fastaFile, read))
 	{
 		getline (fastaFile, read);

 		vector<string> kmers;

 		getKmers(read, &kmers, kmerSize);

 		searchKmers(&kmers, searchFileName, ks);

 	}

 	cout << "nb kmers : " << ks->nbKmers << endl;
    cout << "nb solid kmers : " << ks->nbSolidKmers << endl;

  	fastaFile.close();
    free(ks);
	
	return 0;
}