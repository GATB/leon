#include "testGraphContains.hpp"

using namespace std;

struct GraphStats {

	int nbRequestedKmers;
	int nbContainedKmers;
	int nbMissingKmers;

	//vector<string> kmersRequested;
	//vector<string> kmersContained;
	vector<string>* missingKmers;
};

void graphContains(Graph graph, vector<string>* kmers, struct GraphStats* graphStats, size_t kmerSize, KmerModel* kmerModel)
{
	kmer_type kmer, kmerMin;
	Node node;

	for (std::vector<string>::iterator it = kmers->begin() ; it != kmers->end(); ++it)
	{
    	string kmerString = *it;
		kmer = getKmerType(kmerString.c_str(), kmerModel);
		kmerMin = min(kmer, revcomp(kmer, kmerSize));
		node = Node(Node::Value(kmerMin));

		//cerr << "graphContains - kmer tested : " << kmerString << endl;

		if(graph.contains(node))
		{
			//cerr << "graphContains - kmer " << kmerString << " contained : yes" << endl;
			graphStats->nbContainedKmers += 1;
		}
		else
		{
			//cerr << "graphContains - kmer " << kmerString << " contained : no" << endl;
			graphStats->nbMissingKmers += 1;
			if (find(graphStats->missingKmers->begin(), graphStats->missingKmers->end(), kmerString) == graphStats->missingKmers->end()) {
  			// someName not in name, add it
 			 graphStats->missingKmers->push_back(kmerString);
			}

		}

		graphStats->nbRequestedKmers += 1;
	}
}

void testGraphContains(string inputFileName, Graph graph, size_t kmerSize, KmerModel* kmerModel)
{
	//cerr << "testGraphContains begin" << endl;

	ifstream fastaFile;
  	fastaFile.open (inputFileName);

  	struct GraphStats* graphStats = (struct GraphStats*) malloc (sizeof(struct GraphStats));
  	graphStats->nbRequestedKmers = 0;
	graphStats->nbContainedKmers = 0;
	graphStats->nbMissingKmers = 0;
  	graphStats->missingKmers = new vector<string>();

  	//read lines two by two to skip headers
  	string read;
  	while(getline (fastaFile, read))
 	{
 		getline (fastaFile, read);

 		vector<string> kmers;

 		getKmers(read, &kmers, kmerSize);

 		//ask graph if it contains vector's kmers
 		graphContains(graph, &kmers, graphStats, kmerSize, kmerModel);

 	}

 	cerr << "nb requested kmers : " << graphStats->nbRequestedKmers  << endl;
 	cerr << "nb contained kmers : " << graphStats->nbContainedKmers  << endl;
 	cerr << "nb missing kmers : " << graphStats->nbMissingKmers  << endl;
 	
 	cerr << "missing kmers list : " << endl;
 	std::vector<string>* missingKmers = graphStats->missingKmers;
 	for (std::vector<string>::iterator it = missingKmers->begin() ; it != missingKmers->end(); ++it)
	{
		string kmerString = *it;
		cerr << kmerString << endl;
	}

	delete(graphStats->missingKmers);
	free(graphStats);

	exit(EXIT_FAILURE);
}