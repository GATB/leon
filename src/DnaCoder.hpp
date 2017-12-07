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

#ifndef _DNACODER_HPP_
#define _DNACODER_HPP_


#include <gatb/gatb_core.hpp>
//#include "RangeCoder.hpp"
#include "Leon.hpp"
 #include "Requests.hpp"
//#include "CompressionUtils.hpp"

//#define PRINT_DISTRIB

#include "OrderedBlocks.h"

#ifdef PRINT_DISTRIB
#include <unordered_set>
#include <unordered_map>
#endif

#define NB_MODELS 14

using namespace std;
class Leon;
class Requests;

struct ReadInfos{

	Sequence sequence;
	u_int8_t readType; //comp
	int readSize; //comp
	int anchorPos; //comp
	u_int64_t anchorAddress; //comp
	kmer_type anchor; //comp
	kmer_type revAnchor;
	int revcomp; //comp
	vector<u_int8_t> bifurcations; //comp
	vector<u_int8_t> binaryBifurcations; //comp
	vector<u_int8_t> bifurcationTypes; //comp
	u_int64_t NposCount;
	vector<int> Npos; //comp
	vector<int> leftErrorPos; //comp
	u_int64_t nbLeftError;
	char* cread;
	string sread;


};


struct OrderedReadsInfosBlock{

	kmer_type anchor;
	kmer_type revAnchor;
	u_int32_t nbReads;
};
/*
struct OrderedReadInfos{
	
	u_int8_t readType; //comp
	int readSize; //comp
	int anchorPos; //comp
	//u_int64_t anchorAddress; //comp
	kmer_type anchor; //comp
	int revcomp; //comp
	vector<u_int8_t> bifurcations; //comp
	vector<u_int8_t> binaryBifurcations; //comp
	vector<u_int8_t> bifurcationTypes; //comp
	vector<int> Npos; //comp
	vector<int> leftErrorPos; //comp
	kmer_type revAnchor;
	string sread;
};*/


//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
class AbstractDnaCoder
{
	public:
		AbstractDnaCoder(Leon* leon);

		//debug variables

		ofstream noRangeEncoderOfstream_readTypeModel;
		ifstream noRangeDecoderIfstream_readTypeModel;
		ofstream noRangeEncoderOfstream_anchorAddressModel;
		ifstream noRangeDecoderIfstream_anchorAddressModel;
		ofstream noRangeEncoderOfstream_anchorKmerTypeModel;
		ifstream noRangeDecoderIfstream_anchorKmerTypeModel;
		ofstream noRangeEncoderOfstream_anchorPosModel;
		ifstream noRangeDecoderIfstream_anchorPosModel;
		ofstream noRangeEncoderOfstream_numericDnaModel;
		ifstream noRangeDecoderIfstream_numericDnaModel;
		ofstream noRangeEncoderOfstream_leftErrorModel;
		ifstream noRangeDecoderIfstream_leftErrorModel;
		ofstream noRangeEncoderOfstream_rightErrorModel;
		ifstream noRangeDecoderIfstream_rightErrorModel;
		ofstream noRangeEncoderOfstream_NposModel;
		ifstream noRangeDecoderIfstream_NposModel;
		ofstream noRangeEncoderOfstream_leftErrorPosModel;
		ifstream noRangeDecoderIfstream_leftErrorPosModel;
		ofstream noRangeEncoderOfstream_rightErrorPosModel;
		ifstream noRangeDecoderIfstream_rightErrorPosModel;
		ofstream noRangeEncoderOfstream_readSizeValueModel;
		ifstream noRangeDecoderIfstream_readSizeValueModel;
		ofstream noRangeEncoderOfstream_readAnchorRevcompModel;
		ifstream noRangeDecoderIfstream_readAnchorRevcompModel;
		ofstream noRangeEncoderOfstream_bifurcationModel;
		ifstream noRangeDecoderIfstream_bifurcationModel;
		ofstream noRangeEncoderOfstream_bifurcationBinaryModel;
		ifstream noRangeDecoderIfstream_bifurcationBinaryModel;
		ofstream noRangeEncoderOfstream_noAnchorReadModel;
		ifstream noRangeDecoderIfstream_noAnchorReadModel;
		ofstream noRangeEncoderOfstream_noAnchorReadSizeValueModel;
		ifstream noRangeDecoderIfstream_noAnchorReadSizeValueModel;

		ofstream ofstreams[NB_MODELS];

		enum TypeModel{READ_TYPE_MODEL, 
			ANCHOR_ADDRESS_MODEL,
			ANCHOR_KMER_TYPE_MODEL,
			ANCHOR_POS_MODEL,
			NUMERIC_DNA_MODEL,
			LEFT_ERROR_MODEL,
			N_POS_MODEL,
			LEFT_ERROR_POS_MODEL,
			READ_SIZE_VALUE_MODEL,
			READ_ANCHOR_REVCOMP_MODEL,
			BIFURCATION_MODEL,
			BIFURCATION_BINARY_MODEL,
			NO_ANCHOR_READ_MODEL,
			NO_ANCHOR_READ_SIZE_VALUE_MODEL};



		//end debug
		
	protected:
		Leon* _leon;
		bool _orderReads;
		collections::impl::IBloom<kmer_type>* _bloom; // the bloom containing the solid kmers

		Order0Model _readSizeDeltaTypeModel;
		Order0Model _anchorPosDeltaTypeModel;
		Order0Model _anchorAddressDeltaTypeModel;
		Order0Model _NposDeltaTypeModel;
		Order0Model _errorPosDeltaTypeModel;
		
		Order0Model _readTypeModel; //only 2 value in this model: with anchor or without anchor

		//Sequence* _prevSequences;
		//Order0Model _isPrevReadAnchorableModel;
		//vector<Order0Model> _isPrevReadAnchorablePosModel;
		
		vector<Order0Model> _anchorAddressModel;
		//same use but for for sorted reads
		vector<Order0Model> _anchorKmerTypeModel;
		
		vector<Order0Model> _anchorPosModel;
		
		vector<Order0Model> _numericModel;
		vector<Order0Model> _leftErrorModel;
		vector<Order0Model> _rightErrorModel;
		
		vector<Order0Model> _NposModel;
		vector<Order0Model> _leftErrorPosModel;
		vector<Order0Model> _rightErrorPosModel;
		
		vector<Order0Model> _readSizeValueModel;
		//vector<Order0Model> _nbReadsPerAnchorModel;
		Order0Model _readAnchorRevcompModel;
		Order0Model _bifurcationModel;
		Order0Model _bifurcationBinaryModel;
		
		Order0Model _noAnchorReadModel;
		vector<Order0Model> _noAnchorReadSizeValueModel;
		
		size_t _kmerSize;
		int _readSize;
		KmerModel _kmerModel;
		
		vector<int> _leftErrorPos;
		vector<int> _rightErrorPos;
		vector<int> _Npos;
		
		void startBlock();
		void endRead();
		void codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right);
		void codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right);

		void addErrorPos(int pos, bool rightExtend);
		
		u_int64_t _seqId;
		

	
		u_int64_t _prevReadSize;
		u_int64_t _prevAnchorPos;
		u_int64_t _prevAnchorAddress;
		u_int64_t _prevNpos;
		u_int64_t _prevErrorPos;
		u_int64_t _prevNbLeftError;
	
		int _processedSequenceCount;
};

//====================================================================================
// ** DnaEncoder
//====================================================================================
class DnaEncoder : AbstractDnaCoder
{
		
	public:

		//normally private
		//TODO make it private again
		u_int64_t _readWithoutAnchorCount;
		
		DnaEncoder(Leon* leon);
		DnaEncoder(const DnaEncoder& copy);
		~DnaEncoder();

		void operator()(Sequence& sequence);


		void reset();
		void encodeSortedFileNoAnchorRead(string read);
		void encodeSortedFileAnchor(kmer_type anchor);
		void encodeSortedFileRead(kmer_type anchor, int isRevComp, int readSize, int anchorPos, /*int anchorAddress,*/ vector<int> Npos,
									vector<int> leftErrorPos, vector<u_int8_t> bifurcations, 
									vector<u_int8_t> binaryBifurcations, vector<u_int8_t> bifurcationTypes);
		void encodeSortedFileWriteBlock(int* blockID);
		//old version for ordering reads
		//void encodeReadsInfos(vector< list< struct ReadInfos > > anchorsSequences);
		
	private:
	
	
	//pour quals
	char * _qualseq;
	int * _nb_solids;
	int _smoothing_threshold;
	int _max_read_size;
	bool _trunc_mode;
	
	void storeSolidCoverageInfo();
	void smoothQuals();
	bool apply_smoothing_at_pos(int pos);

	double char2proba(char c);
	char char2phred(char c);

	char * _bufferQuals;
	int _bufferQuals_idx;
	int _bufferQuals_size;
	
#ifdef PRINT_DISTRIB
	vector<Sequence*> _sequences;
	const int maxSequences = 100;
	vector<u_int32_t> _distrib;
	u_int64_t _outDistrib;
#endif
	
	
	
	
		RangeEncoder _rangeEncoder;
		
		#ifdef LEON_PRINT_STAT
			RangeEncoder _rangeEncoder1;
			RangeEncoder _rangeEncoder2;
			RangeEncoder _rangeEncoder3;
			RangeEncoder _rangeEncoder4;
			RangeEncoder _rangeEncoder5;
			RangeEncoder _rangeEncoder6;
		#endif
		
		//static void encodeFirstHeader();
		void writeBlock();
		void execute();
		
		void buildKmers();
		bool isReadAnchorable();
		int findExistingAnchor(u_int32_t* anchorAddress);
		

		void encodeAnchorRead(int anchorPos, u_int32_t anchorAddress);

		kmer_type buildBifurcationList(int pos, kmer_type kmer, bool rightExtend);
		//int buildBifurcationList(int pos, bool rightExtend);
		int voteMutations(int pos, int depth, bool rightExtend);
		
		void encodeNoAnchorRead();
		void saveNoAnchorRead();

		int getBestPath(int pos, kmer_type& kmer, bitset<4>& initRes4, bool rightExtend);

		Sequence* _sequence;
		char* _readseq;
		vector<kmer_type> _kmers;
		KmerModel::Iterator _itKmer;
		vector<u_int8_t> _bifurcations;
		vector<u_int8_t> _binaryBifurcations;
		vector<u_int8_t> _bifurcationTypes;
		
		//bool _isPrevReadAnchorable;
		//u_int64_t _isPrevReadAnchorablePos;

		vector<int> _solidMutaChain;
		int _solidMutaChainPos;
		u_int64_t _totalDnaSize;
		u_int64_t _readCount;
		u_int64_t _MCtotal;
		//u_int64_t _readWithoutAnchorCount;
		u_int64_t _MCuniqSolid;
		u_int64_t _MCuniqNoSolid;
		u_int64_t _MCnoAternative;
		u_int64_t _MCmultipleSolid;
		//u_int64_t _MCmultipleNoSolid;
	
		int _thread_id;

		int _lala;
		int _solidMutaChainStartPos;
		int _solidMutaChainSize;
		int _solidMutaChainLockTime;
};

//====================================================================================
// ** DnaDecoder
//====================================================================================

class DnaDecoder : AbstractDnaCoder
{
		
	public:

		//debug
		int _nbReadTest = 0;
		int _nbAnchorTest = 0;
		//end debug
		
		DnaDecoder(Leon* leon, const string& inputFilename);
		DnaDecoder(Leon* leon, Requests* req, const string& inputFilename);

		~DnaDecoder();
		
		void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount);
		//return false if end of file
		//bool getNextAnchor(string* anchor);
		bool getNextReadInfos(struct ReadInfos* ri);
		bool getNextOrderedReadsInfosBLock(struct OrderedReadsInfosBlock* orib);
		bool getNextOrderedReadInfos(struct ReadInfos* ri);
		
		void execute();
	
		string _buffer;
		bool _finished;

		Requests* _requests;
		bool _decodeReq = false;

	
	private:
	

	
	
		RangeDecoder _rangeDecoder;
		ifstream* _inputFile;
		ofstream* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
		int _decodedSequenceCount;
		string _currentSeq;
		ifstream* _anchorDictFile;
		
		//sortedFile
		kmer_type _anchor;
		void decodeAnchorRead();
		void decodeSortedAnchorRead();

		kmer_type extendAnchor(kmer_type kmer, int pos, bool rightExtend);
		
		void decodeNoAnchorRead();
		void endRead();
		
		int _sequenceCount;
};

class QualDecoder
{
public:
	QualDecoder(Leon* leon, const string& inputFilename);
	~QualDecoder();
	
	void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount);
	void execute();
	
	string _buffer;
	bool _finished;
	
	
private:
	Leon* _leon;

	char * _inbuffer;
	ifstream* _inputFile;
	ofstream* _outputFile;
	u_int64_t _blockStartPos;
	u_int64_t _blockSize;
	int _decodedSequenceCount;
	string _currentSeq;
	int _sequenceCount;
	int _processedSequenceCount;

};
#endif /* _DNACODER_HPP_ */

