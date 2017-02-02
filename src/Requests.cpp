#include "Requests.hpp"

//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
Requests::Requests(IBank* inputBank, string outputFilename, Graph graph, Kmer<>::ModelCanonical model, Partition<kmer_count> & solidCollection, size_t kmerSize){
	cout << "entering requests creator" << endl;

	_inputBank = inputBank;
	_itBank = _inputBank->iterator();
	_itBanks =  _itBank->getComposition();
	_nbBanks = _itBanks.size();

	cout << "nb banks : " << _nbBanks << endl; 

	_graph = graph;
	_model = model;

	_outputFilename = outputFilename; 
	_solidFileSize = solidCollection.getNbItems();
	_nb_kmers_infile = solidCollection.getNbItems();
	_itKmers = solidCollection.iterator();

	//_kmerSize = kmerSize;

	cout << "nb kmers : " << _nb_kmers_infile << endl;

	_signature_array =  (unsigned char  *)  malloc(_solidFileSize*sizeof(char));
    _color_array =  (unsigned char  *)  malloc(_solidFileSize*sizeof(char));


    memset(_signature_array, 0, _solidFileSize);
    memset(_color_array, 0, _solidFileSize);



    const char* signatures_file_path = outputFilename.c_str();
	char* signatures_file_ext = ".signatures_file";
	char signatures_file_path_ext[1024];
	strcpy(signatures_file_path_ext, signatures_file_path);
	strcat(signatures_file_path_ext, signatures_file_ext);
	printf("signatures_file_path_ext : %s", signatures_file_path_ext);

	const char* colors_file_path = outputFilename.c_str();
	char* colors_file_ext = ".colors_file";
	char colors_file_path_ext[1024];
	strcpy(colors_file_path_ext, colors_file_path);
	strcat(colors_file_path_ext, colors_file_ext);

    FILE* signatures_file = fopen(signatures_file_path_ext, "r");
    FILE* colors_file = fopen(colors_file_path_ext, "r");

    fread(_signature_array, 1, _solidFileSize, signatures_file); 

    fread(_color_array, 1, _solidFileSize, colors_file); 

    fclose(signatures_file);
    fclose(colors_file);

}

void Requests::printSignatures(){
	cout << "test signatures" << endl;

	cout << "signatures : \n" << endl;

	for (int i=0; i < _nb_kmers_infile; ++i){
		cout << std::bitset<8>(_signature_array[i]) << endl;
	}
}

void Requests::printColors(){
	cout << "test signatures" << endl;

	cout << "colors : \n" << endl;

	for (int i=0; i < _nb_kmers_infile; ++i){
		cout << std::bitset<8>(_color_array[i]) << endl;
	}
}




void Requests::printSequences(){

    for (_itBank->first(); !_itBank->isDone(); _itBank->next())
	{
		
		std::cout <<  _itBank->item().toString()  << endl << endl;

}
}

void Requests::printKmers(){
for (_itKmers->first(); !_itKmers->isDone(); _itKmers->next()){
			
	std::cout <<  _model.toString (_itKmers->item().getValue())  << endl;
	}

}

void Requests::printMPHFIndexes(){
for (_itKmers->first(); !_itKmers->isDone(); _itKmers->next()){
			
	Node node(Node::Value(_itKmers->item().getValue()));
		printf("_graph.nodeMPHFIndex(node) : %d\n", _graph.nodeMPHFIndex(node));
	}

}

void Requests::printTestAll(){

 for (_itKmers->first(); !_itKmers->isDone(); _itKmers->next())
	{
		
		
	Node node(Node::Value(_itKmers->item().getValue()));
	printf("_graph.nodeMPHFIndex(node) %d\n", _graph.nodeMPHFIndex(node));
		
	std::cout <<  _model.toString (_itKmers->item().getValue())  << "\t" <<
	std::bitset<8>(_color_array[_graph.nodeMPHFIndex(node)]) << "\t" <<
	std::bitset<8>(_signature_array[_graph.nodeMPHFIndex(node)]) << std::endl;

		
	} 
	}