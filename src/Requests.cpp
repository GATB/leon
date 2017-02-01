#include "Requests.hpp"

//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
Requests::Requests(IBank* inputBank, string outputFilename, Graph graph, Partition<kmer_count> & solidCollection){
	cout << "entering requests creator" << endl;

	_inputBank = inputBank;
	_itBank = _inputBank->iterator();
	_itBanks =  _itBank->getComposition();
	_nbBanks = _itBanks.size();

	_graph = graph;

	_outputFilename = outputFilename; //getInput()->getStr(STR_URI_FILE);

	_solidFileSize = solidCollection.getNbItems();
	_nb_kmers_infile = solidCollection.getNbItems();


	//int64_t nb_estimate_kmers = _inputBank->estimateNbItems();
	cout << "nb kmers : " << _nb_kmers_infile << endl;

	_signature_array =  (unsigned char  *)  malloc(_solidFileSize*sizeof(char));
    _color_array =  (unsigned char  *)  malloc(_solidFileSize*sizeof(char));



    /*Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
																solidCollection.iterator(),
																nb_kmers_infile,
																"coloriage"
																);
	LOCAL (itKmers);
	
	_signature_array =  (unsigned char  *)  malloc(solidFileSize*sizeof(unsigned char));
	_color_array = (unsigned char  *)  calloc(solidFileSize,sizeof(unsigned char));
	
	for (itKmers->first(); !itKmers->isDone(); itKmers->next())
	{
	
		uint64_t hashvalue = 	hash1(itKmers->item().getValue(),0);
		
		Node node(Node::Value(itKmers->item().getValue()));
		
		_signature_array[  _graph.nodeMPHFIndex(node) ]  = hashvalue  & 255 ;
		
	//	printf("%u \n",_signature_array[  _graph.nodeMPHFIndex(node) ]);
		
	}*/



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

	cout << "\nheyooooo\n" << endl;

    FILE* signatures_file = fopen(signatures_file_path_ext, "r");
    FILE* colors_file = fopen(colors_file_path_ext, "r");

    fread(_signature_array, 1, _solidFileSize, signatures_file); 

    fread(_color_array, 1, _solidFileSize, colors_file); 

    fclose(signatures_file);
    fclose(colors_file);

}

void Requests::printSignatures(){
	cout << "test signatures" << endl;

	//Test we have all needed data

	cout << "kmers, colors and signatures : \n" << endl;



}














