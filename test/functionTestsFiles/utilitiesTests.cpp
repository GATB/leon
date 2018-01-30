#include "utilitiesTests.hpp"

using namespace std;

kmer_type getKmerType(const char* kmer_chars, KmerModel* _kmerModel){

	kmer_type kmer = _kmerModel->codeSeed(kmer_chars, Data::ASCII).value() ;

	return kmer;
}

void getKmerChars(kmer_type kmer, char* kmer_chars, size_t _kmerSize){

	//char kmer_chars[_kmerSize+1];
	
	const char* kmerStr = kmer.toString(_kmerSize).c_str();
	strncpy(kmer_chars, kmerStr, _kmerSize);
	kmer_chars[_kmerSize] = '\0';

	//return kmer_chars;
}

string getKmerString(kmer_type kmer, size_t _kmerSize){
	return kmer.toString(_kmerSize);
}

Node getKmerNode(char* kmer_chars, KmerModel* _kmerModel){

	kmer_type kmer = getKmerType(kmer_chars, _kmerModel);
	Node node = Node(Node::Value(kmer));

	return node;
}

Node getKmerNode(string kmer_string, KmerModel* _kmerModel){

	kmer_type kmer = getKmerType(kmer_string.c_str(), _kmerModel);
	Node node = Node(Node::Value(kmer));

	return node;
}

Node getKmerNode(kmer_type kmer){

	Node node = Node(Node::Value(kmer));

	return node;
}

string getReverseComplement(string kmer)
{

	string kmerReverseComp;

	for(char& c : kmer) {

		//cout << "c : " << c << endl;

		switch(c){
			case 'A':
			kmerReverseComp.insert(0,"T");
			break;
			case 'C':
			kmerReverseComp.insert(0,"G");
			break;
			case 'G':
			kmerReverseComp.insert(0,"C");
			break;
			case 'T':
			kmerReverseComp.insert(0,"A");
			break;
			case 'a':
			kmerReverseComp.insert(0,"t");
			break;
			case 'c':
			kmerReverseComp.insert(0,"g");
			break;
			case 'g':
			kmerReverseComp.insert(0,"c");
			break;
			case 't':
			kmerReverseComp.insert(0,"a");
			break;
		}
	}

		//cerr << "kmer : " << kmer << endl;
		//cerr << "kmerReverseComp : " << kmerReverseComp << endl;

		//exit(EXIT_FAILURE);


    	//do_things_with(c);
		
	return kmerReverseComp;
}

void getKmers(string read, vector<string>* kmers, int kmerSize)
{
	for (int i=0; i<read.size()-kmerSize+1; ++i)
	{
		#ifdef PRINT_DEBUG
		cout << read.substr(i,kmerSize) << endl;
		#endif
		kmers->push_back(read.substr(i,kmerSize));
	}
}