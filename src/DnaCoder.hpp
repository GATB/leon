
#ifndef _DNACODER_HPP_
#define _DNACODER_HPP_


#include "RangeCoder.hpp"
#include "Leon.hpp"
#include "CompressionUtils.hpp"
#include <gatb/gatb_core.hpp>

using namespace std;
class Leon;

//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
class AbstractDnaCoder
{
	public:
		AbstractDnaCoder(Leon* leon);
		
	protected:
		Leon* _leon;
		collections::impl::Bloom<kmer_type>* _bloom; // the bloom containing the solid kmers
		
};

//====================================================================================
// ** DnaEncoder
//====================================================================================
class DnaEncoder : AbstractDnaCoder
{
		
	public:
		
		DnaEncoder(Leon* leon);
		DnaEncoder(const DnaEncoder& copy);
		~DnaEncoder();

		void operator()(Sequence& sequence);
		
	private:
		
		RangeEncoder _rangeEncoder;
		
		//static void encodeFirstHeader();
		void startBlock();
		void writeBlock();
		void execute();
		
		Sequence* _sequence;
		int _readSize;
		char* _readseq;
		size_t _kmerSize;
	
		OAHash<kmer_type> _anchorKmers;
};

/*
//====================================================================================
// ** DnaDecoder
//====================================================================================
class DnaDecoder : AbstractDnaCoder
{
		
	public:
		
		HeaderDecoder(Leon* leon, ifstream* inputFile, ofstream* outputFile);
		~HeaderDecoder();
		
		//void processNextByte(u_int8_t byte);
		void setup(u_int64_t blockStartPos, u_int64_t blockSize);
	
	private:
		
		RangeDecoder _rangeDecoder;
		ifstream* _inputFile;
		//OutputFile* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
		
		void execute();
		//int _prevPos;
		void endHeader();
		//void decodeFirstHeader();
		void decodeMatch();
		void decodeAscii();
		void decodeNumeric();
		void decodeDelta();
		void decodeZero();
		
		//char _prevHeader2[1000];
		//char _currentHeader2[1000];
		//int _prevHeaderSize;
		
		
};
*/
#endif /* _DNACODER_HPP_ */

