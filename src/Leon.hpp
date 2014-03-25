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



//char char2phred(char c);
//double char2proba(char c);

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
		std::string     _solidFile;
		static const int READ_PER_BLOCK = 100000;
		int _nb_cores;
		
		//pthread_mutex_t writer_mutex;
		//pthread_cond_t buffer_full_cond;
		//int _next_writer_thread_id;
		
		//Header compression
		string _firstHeader;
		void writeBlock(u_int8_t* data, u_int64_t size);
		u_int64_t _totalHeaderSize;
		u_int64_t _totalHeaderCompressedSize;
		
		//Dna compression
		u_int64_t _totalDnaSize;
		u_int64_t _totalDnaCompressedSize;
		collections::impl::Bloom<kmer_type>* _bloom;
		//test dna compression to remove
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
	
	private:
	 
		//static const char* STR_GZ;
		IFile* _outputFile;
		bool _isFasta;

		void execute ();
		virtual collections::impl::Bloom<kmer_type>* createBloom ();
		
		//Global compression
		Order0Model _generalModel;
		RangeEncoder _rangeEncoder;
		vector<u_int64_t> _blockSizes;
		
		void executeCompression();
		void executeDecompression();
		
		//Global decompression
		RangeDecoder _rangeDecoder;
		ifstream* _inputFile;
		
		//Header compression
		Order0Model _headerNumericSizeModel;
		vector<Order0Model> _headerNumericModels;
		void startHeaderCompression(Iterator<Sequence>* itSeq);
		void endHeaderCompression();
		
		//DNA Compression
		void startDnaCompression(Iterator<Sequence>* itSeq);
		void endDnaCompression();
		
		//Header decompression
		//HeaderDecoder _headerDecoder;
		
};

/*
class CompressReads
{
	public:
	
		CompressReads(Bloom<kmer_type>* bloom, Leon * leon, u_int64_t* nb_solids_kmers, int nb_cores, int * nbliving);
		CompressReads(const CompressReads& cr);
		~CompressReads();
		
		void operator() ( Sequence& sequence);
		
	private:
		
		Bloom<kmer_type> * _bloom; // the bloom containing the solid kmers
		Leon *         _leon; // the parent Leon object
		
		int *  _nb_living;
		int _thread_id;
		
		ISynchronizer* _synchro;
		ISynchronizer* getSynchro ()  { return _synchro; }
		
		u_int64_t *  _total_nb_solid_kmers_in_reads;
		u_int64_t   _local_nb_solid_kmers_in_reads;
		
};*/

#endif /* defined(__leon__Leon__) */
