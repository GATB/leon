/*****************************************************************************
 *   Leon: reference free compression for NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Benoit, G.Rizk, C.Lemaitre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef __leon__Leon__
#define __leon__Leon__


#define LEON_VERSION_MAJOR 1
#define LEON_VERSION_MINOR 0
#define LEON_VERSION_PATCH 0


//#define LEON_PRINT_STAT

#include <iostream>
#include <gatb/gatb_core.hpp>
#include <sys/time.h>

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






#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>
#include <list>
#include <stdio.h>
#include <sys/wait.h>
#include <cstdlib>
#include "HeaderCoder.hpp"
#include "DnaCoder.hpp"
#include "Requests.hpp"

#include "OrderedBlocks.h"

//#include "RangeCoder.hpp"

#include <time.h> //Used to calculate time taken by decompression
#include <zlib.h> //Test bloom compression
//char char2phred(char c);
//double char2proba(char c);


#include <pthread.h>

class HeaderEncoder;
class HeaderDecoder;
class DnaEncoder;

class Leon : public misc::impl::Tool
{
	public:

		//DEBUG VARIABLES
		//TO ERASE
		int debug = 0;
		//normally private
		ifstream* _inputFile;

		//Path to binary
		std::string _binaryPath;// = "/home/tbraquel/NGS/gatb-tool-leon/";
		
		//Leon( bool compress, bool decompress);
		Leon(string binaryPath);
		~Leon();
	
		static const char* STR_COMPRESS;
		static const char* STR_DECOMPRESS;
		static const char* STR_TEST_DECOMPRESSED_FILE;
		static const char* STR_DNA_ONLY;
		static const char* STR_NOHEADER;
		static const char* STR_NOQUAL;

		//add
		static const char* STR_OUTPUT_FILE;

		//requests
		static const char* STR_REQUEST;
		static const char* STR_SEARCH_SEQUENCE;

	
		size_t          _kmerSize;
		string     _h5OutputFilename;
		static const int READ_PER_BLOCK = 50000;
		int _nb_cores;
		
		bool _compress, _decompress, _request;
		
		clock_t _time; //Used to calculate time taken by decompression
		
		//Global compression
		void writeBlock(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID);
	
		void writeBlockLena(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID);

		//Peacock file Compression (for requests)
		//void writePeacockBlock(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID);

		//Header compression
		string _firstHeader;
		u_int64_t _totalHeaderSize;
		//u_int64_t _totalHeaderCompressedSize;
		
		//Dna compression
		u_int32_t _anchorAdress;
		
		bool anchorExist(const kmer_type& kmer, u_int32_t* anchorAdress);
		int findAndInsertAnchor(const vector<kmer_type>& kmers, u_int32_t* anchorAdress);
		
		u_int64_t _totalDnaSize;
		u_int64_t _anchorDictSize;
		u_int64_t _anchorAdressSize;
		u_int64_t _anchorPosSize;
		u_int64_t _otherSize;
		u_int64_t _readSizeSize;
		u_int64_t _bifurcationSize;
		u_int64_t _noAnchorSize;
		
		//u_int64_t _totalDnaCompressedSize;
		//u_int64_t _realDnaCompressedSize;
		u_int64_t _compressedSize;
		//IBloom<kmer_type>* _bloom;
		Graph _graph;
	
		bool _isFasta;
		bool _noHeader;

		bool _orderReads;
		unsigned long _nbLinesToSort = 0;
		int _nbReducers = 4; //default 4

		//tmp tests bool
		bool _readSortedFileTest;

		ofstream unsortedReads;
		ofstream noAnchorReadsOfstream;
		ifstream noAnchorReadsIfstream;
		//list<struct ReadInfos>* anchorsSequences;
		vector< list< struct ReadInfos > > anchorsSequences;


	bool _lossless;
	//for qual compression
		u_int64_t _total_nb_quals_smoothed ;
		u_int64_t _input_qualSize;
		u_int64_t _compressed_qualSize;
	
	
		//test dna compression
		u_int64_t _MCtotal;
	
		u_int64_t _MCnoAternative;
		u_int64_t _MCuniqSolid;
		u_int64_t _MCuniqNoSolid;
		u_int64_t _MCmultipleSolid;
		//u_int64_t _MCmultipleNoSolid;
	
		int _coverage;
		u_int64_t _blockCount;
		//u_int64_t _noAnchor_full_N_kmer_count;
		//u_int64_t _noAnchor_with_N_kmer_count;
		u_int64_t _readCount;
		//double _anchorKmerCount;
		u_int64_t _readWithoutAnchorCount;
		//double _total_kmer_indexed;
		//double _uniq_mutated_kmer;
		//u_int64_t _total_kmer;
		//u_int64_t _readWithoutAnchorSize;
		//u_int64_t _readWithAnchorSize;
		//double _readWithAnchorMutationChoicesSize;
		//Hash16<kmer_type>* _kmerAbundance;
		
		int   _nb_thread_living;
		OrderedBlocks * _blockwriter;

		OrderedBlocks * _qualwriter;

		// ProgressSynchro *
		dp::IteratorListener * _progress_decode;
	
		void setBlockWriter (OrderedBlocks* blockwriter) { SP_SETATTR(blockwriter); }

		//DNA decompression
		kmer_type getAnchor(ifstream* anchorDictFile, u_int32_t adress);
		string _anchorDictFilename;
		
		
		static const int nt2binTab[128];
		static const int bin2ntTab[5];
		//static const vector<int> bin2ntTab(5;
		
		

		//Utils
		static int nt2bin(char nt){
			return nt2binTab[nt];
			
		}
		static int bin2nt(int nt){
			return bin2ntTab[nt];
		}
		
	
	void decodeSortedAnchorDict();

	private:

		unsigned char  * _signature_array;
		unsigned char  * _color_array;
	
	
	
	
	
		u_int64_t _lastAnchorValue;
		
		 struct timeval _tim;
		double _wdebut_leon, _wfin_leon;
		//static const char* STR_GZ;
		IFile* _outputFile;
		IFile* _outputFileRequests;
		ofstream* _dictAnchorFile;
		int _nks;
		
		void execute ();
		void createBloom ();
		void coloriage();
		//void createKmerAbundanceHash();
		
		//Global compression
		string _inputFilename;
		string _baseOutputname;
		string _outputFilename;
		string _outputFileRequestsName;
	

	//quals
	string _FileQualname;
	IFile* _FileQual;
	string _qualOutputFilename; //temp file
	ofstream* _qualOutputFile;

		Order0Model _generalModel;
		vector<Order0Model> _numericModel;
		vector<Order0Model> _nbReadsPerAnchorModel;
		RangeEncoder _rangeEncoder;
		vector<u_int64_t> _blockSizes;
	
		vector<u_int64_t> _qualBlockSizes;
		vector<u_int64_t> _headerBlockSizes;
		vector<u_int64_t> _dnaBlockSizes;

		IBank* _inputBank;
		void setInputBank (IBank* inputBank) { SP_SETATTR(inputBank); }

		//u_int64_t _bloomSize;
		
		void executeCompression();
		void executeDecompression();
		void executeRequest();
		void endCompression();
		//void endPeacockCompression();
		void endQualCompression();
	
		//Global decompression
		void setupNextComponent(vector<u_int64_t>   & blockSizes  );
		
		RangeDecoder _rangeDecoder;
		//ifstream* _inputFile;
		//ifstream* _inputFile;
		ifstream* _inputFileQual;

		ifstream* _descInputFile;
		u_int64_t _filePos;
	
	u_int64_t _filePosHeader;
	u_int64_t _filePosDna;

		double _headerCompRate, _dnaCompRate, _qualCompRate;
		
		//Quals
		u_int64_t _filePosQual;

	
		void startDecompressionAllStreams();

		
	
		//Header compression
		void startHeaderCompression();
		void endHeaderCompression();
	

		//DNA Compression
		void startDnaCompression();
		void endDnaCompression();
		void writeBloom();
		void writeAnchorDict();
		void writeSortedAnchorDict();
		void encodeInsertedAnchor(const kmer_type& kmer);
		void encodeInsertedSortedAnchor(const kmer_type& anchor, u_int32_t nbReads);
			
		RangeEncoder _anchorRangeEncoder;
		//RangeEncoder _sortedAnchorRangeEncoder;
		Order0Model _anchorDictModel;
		
		//map<kmer_type, u_int32_t> _anchorKmers; //uses 46 B per elem inserted
		//OAHash<kmer_type> _anchorKmers;
		Hash16<kmer_type, u_int32_t >  * _anchorKmers ; //will  use approx 20B per elem inserted
		Hash16<kmer_type, u_int32_t >  * _anchorKmersSorted;

		//tmp test decompression
		Hash16<kmer_type, u_int32_t >  * _anchorKmersSortedD;

		//Header decompression
	
		string _headerOutputFilename;
		ofstream* _headerOutputFile;
		
	   	int _auto_cutoff;
		pthread_mutex_t findAndInsert_mutex;
		pthread_mutex_t writeblock_mutex;

		//DNA Decompression
		void decodeBloom();
		void decodeAnchorDict();
		
		KmerModel* _kmerModel;
		string _dnaOutputFilename;
		ofstream* _dnaOutputFile;
		RangeDecoder _anchorRangeDecoder;
		RangeDecoder _sortedAnchorRangeDecoder;
		vector<kmer_type> _vecAnchorKmers;
		
		//Global decompression
		void endDecompression();
		
		//IFile* _outputFile;
};


#endif /* defined(__leon__Leon__) */
