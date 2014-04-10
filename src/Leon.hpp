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
		IBank* _testBank; //debug
		Iterator<Sequence>* _testBankIt;
		kmer_type getAnchor(u_int32_t adress);
		
		//Utils
		static int nt2bin(char nt){
			/*
			int i;
			i = nt;
			i = (i>>1)&3; // that's quite clever, guillaume.
			return i;
			*/
			
			if(nt == 'A')
				return 0;
			else if(nt == 'C')
				return 1;
			else if(nt == 'T')
				return 2;
			else if(nt == 'G')
				return 3;
			else if(nt == 'N'){
				//cout << "error nt2bin N" << endl;
				return 4;
			}
				
		}
		
		static int bin2nt(int nt){
			//if(nt == 4) cout << "error bin2nt N" << endl;
			static char tab[5] = {'A', 'C', 'T', 'G', 'N'};
			return tab[nt];
			
			/*
			if(nt == 0)
				return 'A';
			else if(nt == 1)
				return 'C';
			else if(nt == 2)
				return 'G';
			else if(nt == 3)
				return 'T';
			else if(nt == 4)
				return 'N';*/
		}
		
		
	private:
	 
		//static const char* STR_GZ;
		IFile* _outputFile;
		bool _isFasta;
		int _nks;
		
		void execute ();
		void createBloom ();
		void createKmerAbundanceHash();
		
		//Global compression
		string _inputFilename;
		string _outputFilename;
		Order0Model _generalModel;
		Order0Model _numericSizeModel;
		vector<Order0Model> _numericModel;
		RangeEncoder _rangeEncoder;
		vector<u_int64_t> _blockSizes;
		IBank* _inputBank;
		u_int64_t _bloomSize;
		
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
		//HeaderDecoder _headerDecoder;
		
		//DNA Decompression
		void decodeBloom();
		void decodeAnchorDict();
		
		RangeDecoder _anchorRangeDecoder;
		vector<kmer_type> _vecAnchorKmers;
};


#endif /* defined(__leon__Leon__) */
