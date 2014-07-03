//
//  Leon.hpp
//  leon
//
//  Created by Guillaume Rizk on 17/01/2014.
//  Copyright (c) 2014 G.Rizk, R.Uricaru. All rights reserved.
//

#ifndef __leon__Leon__
#define __leon__Leon__


#include <iostream>
#include <gatb/gatb_core.hpp>

/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;


typedef kmer::impl::Kmer<>::Model KmerModel;
typedef kmer::impl::Kmer<>::Type  kmer_type;
typedef kmer::impl::Kmer<>::Count kmer_count;


#include <string>
#include <sstream>
#include "HeaderCoder.hpp"
#include "DnaCoder.hpp"
//#include "RangeCoder.hpp"

#include <time.h> //Used to calculate time taken by decompression
#include <zlib.h> //Test bloom compression
//char char2phred(char c);
//double char2proba(char c);

#include <thread>
#include <future>

class HeaderEncoder;
class HeaderDecoder;
class DnaEncoder;

class Leon : public misc::impl::Tool
{
	public:
		
		Leon();
		
		static const char* STR_COMPRESS;
		static const char* STR_DECOMPRESS;
		
		size_t          _kmerSize;
		string     _dskOutputFilename;
		static const int READ_PER_BLOCK = 30000;
		int _nb_cores;
		
		clock_t _time; //Used to calculate time taken by decompression
		
		//pthread_mutex_t writer_mutex;
		//pthread_cond_t buffer_full_cond;
		//int _next_writer_thread_id;
		
		//Header compression
		string _firstHeader;
		void writeBlock(u_int8_t* data, u_int64_t size);
		u_int64_t _totalHeaderSize;
		u_int64_t _totalHeaderCompressedSize;
		
		//Dna compression
		u_int32_t _anchorAdress;
		
		bool anchorExist(const kmer_type& kmer, u_int32_t* anchorAdress);
		int findAndInsertAnchor(const vector<kmer_type>& kmers, u_int32_t* anchorAdress);
		
		u_int64_t _totalDnaSize;
		u_int64_t _totalDnaCompressedSize;
		u_int64_t _realDnaCompressedSize;
		Bloom<kmer_type>* _bloom;
		//test dna compression to remove			
		double _MCtotal;
		double _MCnoAternative;
		double _MCuniqSolid;
		double _MCuniqNoSolid;
		double _MCmultipleSolid;
		double _MCmultipleNoSolid;
		u_int64_t _noAnchor_full_N_kmer_count;
		u_int64_t _noAnchor_with_N_kmer_count;
		u_int64_t _readCount;
		double _anchorKmerCount;
		double _readWithoutAnchorCount;
		double _total_kmer_indexed;
		double _uniq_mutated_kmer;
		u_int64_t _total_kmer;
		u_int64_t _readWithoutAnchorSize;
		u_int64_t _readWithAnchorSize;
		double _readWithAnchorMutationChoicesSize;
		OAHash<kmer_type>* _kmerAbundance;
		
		//DNA decompression
		kmer_type getAnchor(u_int32_t adress);
		
		
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
		
		
	private:
	 
		
		//static const char* STR_GZ;
		IFile* _outputFile;
		ofstream* _dictAnchorFile;
		int _nks;
		
		void execute ();
		void createBloom ();
		void createKmerAbundanceHash();
		
		//Global compression
		bool _isFasta;
		string _inputFilename;
		string _outputFilename;
		Order0Model _generalModel;
		Order0Model _numericSizeModel;
		vector<Order0Model> _numericModel;
		RangeEncoder _rangeEncoder;
		vector<u_int64_t> _blockSizes;
		IBank* _inputBank;
		//u_int64_t _bloomSize;
		
		void executeCompression();
		void executeDecompression();
		void endCompression();
		
		//Global decompression
		RangeDecoder _rangeDecoder;
		//ifstream* _inputFile;
		ifstream* _inputFile;
		ifstream* _descInputFile;
		u_int64_t _filePos;
		double _headerCompRate, _dnaCompRate;
		
		void setupNextComponent();
		
		//Header compression
		void startHeaderCompression();
		void endHeaderCompression();
		
		//DNA Compression
		void startDnaCompression();
		void endDnaCompression();
		void writeBloom();
		void writeAnchorDict();
		void encodeInsertedAnchor(const kmer_type& kmer);
			
		RangeEncoder _anchorRangeEncoder;
		Order0Model _anchorDictModel;
		
		map<kmer_type, u_int32_t> _anchorKmers;
		//OAHash<kmer_type> _anchorKmers;
		
		//Header decompression
		void startHeaderDecompression();
		
		string _headerOutputFilename;
		ofstream* _headerOutputFile;
		
		//DNA Decompression
		void startDnaDecompression();
		void decodeBloom();
		void decodeAnchorDict();
		
		string _dnaOutputFilename;
		ofstream* _dnaOutputFile;
		RangeDecoder _anchorRangeDecoder;
		vector<kmer_type> _vecAnchorKmers;
		
		//Global decompression
		void endDecompression();
		
		//IFile* _outputFile;
};


#endif /* defined(__leon__Leon__) */
