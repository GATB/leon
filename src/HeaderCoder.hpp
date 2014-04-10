
#ifndef _HEADERCODER_HPP_
#define _HEADERCODER_HPP_

/*
#include <string>
#include <vector>
#include <ctype.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>*/

#include "RangeCoder.hpp"
#include <gatb/gatb_core.hpp>
#include "Leon.hpp"
#include "CompressionUtils.hpp"

using namespace std;
class Leon;
//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
class AbstractHeaderCoder
{
	public:
		AbstractHeaderCoder(Leon* leon);
		
	protected:
		void addFieldColumn();
	
		enum HeaderType{HEADER_END=1, FIELD_ASCII, FIELD_NUMERIC, FIELD_DELTA, FIELD_ZERO_ONLY, FIELD_ZERO_AND_NUMERIC, HEADER_TYPE_COUNT};
		//static const int MAX_FIELD_COUNT = 200;
		
		vector<Order0Model> _typeModel;
		vector<Order0Model> _fieldIndexModel;
		vector<Order0Model> _fieldColumnModel;
		vector<Order0Model> _misSizeModel;
		vector<Order0Model> _asciiModel;
		vector<Order0Model> _numericSizeModel;
		vector< vector<Order0Model> > _numericModels;
		vector<Order0Model> _zeroModel;
		
		int typeOfChar(u_int8_t c, bool* isDigit);
		void splitHeader();
		void makeField();
		void endHeader();
		
		string _prevHeader;
		string _currentHeader;
		vector<unsigned int> _prevFieldPos;
		vector<unsigned int> _currentFieldPos;
		int _currentPos;
		int _fieldStartPos;
		int _prevFieldCount;
		int _fieldIndex;
		int _misIndex;
		vector<u_int64_t> _prevFieldValues;
		vector<u_int64_t> _currentFieldValues;
		vector<u_int64_t> _prevFieldZeroValues;
		vector<u_int64_t> _currentFieldZeroValues;
		vector<HeaderType> _prevFieldTypes;
		vector<HeaderType> _currentFieldTypes;
		
		
		bool _isCurrentFieldNumeric;
		int _currentFieldCount;
		
		Leon* _leon;
		
		void startBlock();
		
};

//====================================================================================
// ** HeaderEncoder
//====================================================================================
class HeaderEncoder : AbstractHeaderCoder
{
		
	public:
		
		HeaderEncoder(Leon* leon);
		HeaderEncoder(const HeaderEncoder& copy);
		~HeaderEncoder();

		void operator()(Sequence& sequence);
		int getId();
		u_int64_t _lastSequenceIndex;
		
	private:
		
		RangeEncoder _rangeEncoder;
		
		int _fieldPos;
		//int _misPrevStartPos, _misCurrentStartPos;
		int _misCurrentStartPos;
		//int _encoderFieldIndex;
		int _prevFieldSize, _currentFieldSize;
		int _lastMatchFieldIndex;
		
		//static void encodeFirstHeader();
		void writeBlock();
		
		void processNextHeader();
		void compareHeader();
		//void encode();
		//void encodeMismatch();
		void encodeNumeric();
		void encodeAscii();
		
};

//====================================================================================
// ** HeaderDecoder
//====================================================================================
class HeaderDecoder : AbstractHeaderCoder
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

#endif /* _HEADERCODER_HPP_ */

