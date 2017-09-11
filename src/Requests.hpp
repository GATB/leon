#ifndef __REQUESTS__HPP__
#define __REQUESTS__HPP__



#include <iostream>
#include <gatb/gatb_core.hpp>
#include "../thirdparty/gatb-core/gatb-core/src/gatb/gatb_core.hpp"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/Model.hpp"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/gatb_core.hpp"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/tools/misc/impl/Tool.hpp"
#include <sys/time.h>

#include "DnaCoder.hpp"
#include "OrderedBlocks.h"
#include "HeaderCoder.hpp"

#include "OrderedBlocks.h"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/tools/compression/RangeCoder.hpp"

/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;


typedef kmer::impl::Kmer<>::ModelDirect KmerModel;
typedef kmer::impl::Kmer<>::Type        kmer_type;
typedef kmer::impl::Kmer<>::Count       kmer_count;

#include <string>
#include <vector>
#include <list>
#include <ctype.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

//#include "RangeCoder.hpp"
#include "Leon.hpp"
#include <bitset>

//#include "CompressionUtils.hpp"

using namespace std;
class Leon;
class DnaDecoder;
class HeaderDecoder;
class QualDecoder;


//class HeaderEncoder;
//class HeaderDecoder;
//class DnaEncoder;

//class Requests : public misc::impl::Tool

#define NB_MAX_COLORS ((size_t) 8)
//the minimum number of required matches with a read 
#define NB_MIN_MATCHES 0
#define NB_MAX_SNPS 0
//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
class Requests
{
	public:

		/*****constructor*****/
		Requests(IBank* inputBank, string outputFilename, Graph graph, 
			Kmer<>::ModelCanonical model, 
			Partition<kmer_count> & solidCollection, size_t kmerSize, 
			Hash16<kmer_type, 
			u_int32_t >  * anchorKmers,
			Hash16<kmer_type, 
			u_int32_t >  * anchorKmersSorted,
			Leon* leon,
			DnaDecoder* dnadec/*, 
			Order0Model generalModel,
			vector<Order0Model> numericModel,
			Order0Model anchorDictModel*/);
		
		/*****utilities*****/

		/*
		**-writes the nth kmer of the sequence in kmer
		**-returns false if end of sequence
		*/
		bool getSequenceKmer(const char* seq, uint n, char* kmer);
		/*
		**-writes the next anchor of the sequence in anchor, starting from 
		**pos to the right extent
		**-writes the address of the anchor in anchorAddress
		**-returns false if end of sequence
		*/
		bool getNextAnchor(char* sequence, uint* pos, char* anchor, u_int32_t anchorAddress);
		void fillSequenceAnchorsDict(Hash16<kmer_type, list<u_int32_t>* >  * sequenceAnchorKmers,
										char* sequence);
		/*
		**delete the created lists of the hash table
		*/
		void emptySequenceAnchorDict(Hash16<kmer_type, list<u_int32_t>* >  * sequenceAnchorKmers, char* sequence);
		/*
		**-returns true if kmer is in anchor dictionnary
		*/
		bool anchorExist(char* kmer_chars, u_int32_t* anchorAddress);
		/*
		**-returns the signature of the kmer
		*/
		unsigned char getKmerSignature(kmer_type kmer);
		bitset<NB_MAX_COLORS> getReadColor(struct ReadInfos* ri); 


		//conversion functions
		kmer_type getKmerType(char* kmer_chars);
		void getKmerChars(kmer_type kmer, char* kmer_chars);
		void getSequenceChars(char* sequence_chars, string sequence);
		string getKmerString(kmer_type kmer);
		Node getKmerNode(char* kmer_chars);
		Node getKmerNode(kmer_type kmer);
		int getNodeMPHFIndex(kmer_type kmer);
		int getNodeMPHFIndex(char* kmer_chars);

		// decode functions
		//leon methods copies
		//solution :
		// - will have to call leon's method before creating requests ??
		void initializeRangeDecoder();
		void clearRangeDecoder();
		void initializeDecodingReads();
		void setupNextComponent(vector<u_int64_t> & blockSizes);
		void decodeBloom();
		void decodeAnchorDict();
		void decodeSortedAnchorDict();
		void decodeInfos();
		void headerSetUp();
		void dnaSetUp();
		void initializeDecoders();
		void clearDecoders();
		void headerDecoderSetup(int blockIndice);
		void dnaDecoderSetup(int blockIndice);
		void qualDecoderSetup(int blockIndice);	

		//TODO :
		//remove this (called in dnadecoder, initially a leon method)
		kmer_type getAnchor(ifstream* anchorDictFile, u_int32_t adress);




		/*****query*****/
		void fgetRequests();
		void fgetString(char* string, int stringLen, char* query);
		bool fgetKmer(char* kmer);
		bool fgetSequence(char* sequence);

		/*****debug*****/
		void printSignatures();
		void printColors();
		void printKmers();
		void printSequences();
		void printMPHFIndexes();

		//leon
		void printSequenceAnchors(char* sequence);
		void printIsKmerInSequenceAnchorDict(char* kmer_chars, Hash16<kmer_type, list<u_int32_t>* >* sequenceAnchorKmers);
		void printSequenceAnchorsDict(char* sequence, Hash16<kmer_type, list<u_int32_t>* >* sequenceAnchorKmers);
		void printTestAll();
		void testPrintReadsFile(bool getReads, bool getAnchors, bool getAnchorPos);
		void testPrintAllHeadersReadsFile();

		//peacock
		void testPrintPFile();

		/*****requests*****/

		void printNbBanks();
		int getNbBanks();

		void printKmerSize();
		int getKmerSize();

		// requests on kmers
		// (if kmer is in graph, it is on data)

		bool isKmerInGraph(char* kmer);
		bool isKmerInGraph(kmer_type kmer);
		bitset<NB_MAX_COLORS> getKmerColors(char* kmer);
		bitset<NB_MAX_COLORS> getKmerColors(kmer_type kmer);
		int getKmerNbColors(char* kmer);

		// First step for sequences requests (or rapid search)
		// request the graph, if not in the graph, no reason to search in data
		
		bool isSequenceInGraph(char* sequence);
		bitset<NB_MAX_COLORS> getSequenceColorsInGraph(char* sequence);
		int getSequenceNbColorsInGraph(char* sequence);

		// Second step for sequences requests
		// search in the data, extending anchor kmers

		bool isSequenceInData(char* sequence);
		bitset<NB_MAX_COLORS> getSequenceColorsInData(char* sequence);
		int getSequenceNbColorsInData(char* sequence);
		void getSequenceFileMatchesInData(char* sequence, vector<bitset<NB_MAX_COLORS>>* sequenceMatches);

		//request global variables
		bool _orderReads;
		char request[1024];
		int req_buffer_size;
		bool end_requests;
		int sequenceMaxSize;

		Hash16<kmer_type, list<u_int32_t>*>* _sequenceAnchorKmers;
		// 1 if a color can be fullfilled using the reads of the dataset of a
		// same color
		bitset<NB_MAX_COLORS> _sequenceColors;
		// 1 if a color may have been filled using the same read at least twice, 0 else
		bitset<NB_MAX_COLORS> _sequenceAmbiguousMatches;
		// 1 if a color may have been fullfilled using the same read at least twice, 0 else
		bitset<NB_MAX_COLORS> _sequenceAmbiguousColors;
		vector<bitset<NB_MAX_COLORS>>* _sequenceMatches;

		// variables from Leon

		IBank* _inputBank;
		Iterator<Sequence>* _itBank;
		std::vector<Iterator<Sequence>*> _itBanks;
		Iterator<Sequence>* _itSubBank;
		int _nbBanks;

		string _outputFilename;

		//Kmer<>::ModelCanonical::Iterator _itKmer;

		Leon* _leon;
		Graph _graph;
		size_t _kmerSize;
		Kmer<>::ModelCanonical _model;
		KmerModel* _kmerModel;
		Hash16<kmer_type, u_int32_t >  * _anchorKmers;
		Hash16<kmer_type, u_int32_t >  * _anchorKmersSorted;
		Hash16<kmer_type, u_int32_t >  * _anchorKmersSortedD; //get at decompression
		u_int32_t _anchorAdress;
		ifstream* _anchorDictFile;


		//Decode requirements
		string _decodeFilename;
		RangeDecoder _rangeDecoder;
		RangeDecoder _anchorRangeDecoder;
		//RangeDecoder _sortedAnchorRangeDecoder;
		ifstream* _descInputFile;
		ifstream* _inputFile;
		DnaDecoder* _dnadec;
		vector<u_int64_t> _headerBlockSizes;
		vector<u_int64_t> _dnaBlockSizes;
		vector<kmer_type> _vecAnchorKmers;
		u_int64_t _filePos;
		u_int64_t _blockCount;
		static const int READ_PER_BLOCK = 50000;
		bool _isFasta;
		bool _noHeader;
		HeaderDecoder* _hdecoder;
		DnaDecoder* _ddecoder;
		QualDecoder* _qdecoder;
		u_int64_t _filePosHeader = 0;
		u_int64_t _filePosDna = 0;

		// models :
		Order0Model _generalModel;
		vector<Order0Model> _numericModel;
		Order0Model _anchorDictModel;
		vector<Order0Model> _nbReadsPerAnchorModel;
		//uint64_t blablalol;
		
		//int a1;
		//int a2;
		//int tab[5];
		//int16_t b;
		//int a3;
		/*int a4;
		int a5;
		int a6;

*/		/*Iterator<kmer_count>* _itKmersAll;
		Kmer<>::ModelCanonical::Iterator _itKmersSubBank;*/
		Iterator<kmer_count>* _itKmers;

		u_int64_t _solidFileSize;
		u_int64_t _nb_kmers_infile;

		unsigned char* _signature_array;
		unsigned char* _color_array;
		

	protected:
};

#endif __REQUESTS__HPP__ 

