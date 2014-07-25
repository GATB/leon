
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
		collections::impl::IBloom<kmer_type>* _bloom; // the bloom containing the solid kmers
		
		Order0Model _readSizeDeltaTypeModel;
		Order0Model _anchorPosDeltaTypeModel;
		Order0Model _anchorAddressDeltaTypeModel;
		Order0Model _NposDeltaTypeModel;
		Order0Model _errorPosDeltaTypeModel;
		
		Order0Model _readTypeModel; //only 2 value in this model: with anchor or without anchor
		
		Order0Model _anchorAddressSizeModel;
		vector<Order0Model> _anchorAddressModel;
		
		Order0Model _anchorPosSizeModel;
		vector<Order0Model> _anchorPosModel;
		
		Order0Model _numericSizeModel;
		vector<Order0Model> _numericModel;
		
		Order0Model _NposSizeModel;
		vector<Order0Model> _NposModel;
		Order0Model _errorPosSizeModel;
		vector<Order0Model> _errorPosModel;
		
		Order0Model _readSizeModel;
		vector<Order0Model> _readSizeValueModel;
		Order0Model _readAnchorRevcompModel;
		Order0Model _mutationModel;
		
		Order0Model _noAnchorReadModel;
		Order0Model _noAnchorReadSizeModel;
		vector<Order0Model> _noAnchorReadSizeValueModel;
		
		size_t _kmerSize;
		int _readSize;
		KmerModel _kmerModel;
		
		vector<int> _errorPos;
		vector<int> _Npos;
		
		void startBlock();
		void endRead();
		void codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right);
		void codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right);
		
		u_int64_t _seqId;
		

	
		u_int64_t _prevReadSize;
		u_int64_t _prevAnchorPos;
		u_int64_t _prevAnchorAddress;
		u_int64_t _prevNpos;
		u_int64_t _prevErrorPos;
	
		int _processedSequenceCount;
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
		
		#ifdef LEON_PRINT_STAT
			RangeEncoder _rangeEncoder1;
			RangeEncoder _rangeEncoder2;
			RangeEncoder _rangeEncoder3;
			RangeEncoder _rangeEncoder4;
			RangeEncoder _rangeEncoder5;			
		#endif
		
		//static void encodeFirstHeader();
		void writeBlock();
		void execute();
		
		void buildKmers();
		int findExistingAnchor(u_int32_t* anchorAddress);
		
		void encodeAnchorRead(int anchorPos, u_int32_t anchorAddress);
		kmer_type buildBifurcationList(int pos, kmer_type kmer, bool rightExtend);
		//int buildBifurcationList(int pos, bool rightExtend);
		int voteMutations(int pos, bool rightExtend);
		
		void encodeNoAnchorRead();
		
		
		Sequence* _sequence;
		char* _readseq;
		vector<kmer_type> _kmers;
		KmerModel::Iterator _itKmer;
		vector<int> _bifurcations;
		
		vector<int> _solidMutaChain;
		int _solidMutaChainPos;
		u_int64_t _totalDnaSize;
		u_int64_t _readCount;
		u_int64_t _MCtotal;
		u_int64_t _readWithoutAnchorCount;
		u_int64_t _MCuniqSolid;
		u_int64_t _MCuniqNoSolid;
		u_int64_t _MCnoAternative;
		u_int64_t _MCmultipleSolid;
		u_int64_t _MCmultipleNoSolid;
	
		
		bool extendMutaChain(kmer_type kmer, int pos, bool rightExtend);
		bool extendMutaChainRec(vector< vector< vector<kmer_type> > >& mutaChains, bool rightExtend);
		
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
		
		DnaDecoder(Leon* leon, const string& inputFilename);
		~DnaDecoder();
		
		void setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount);
		void execute();
	
		string _buffer;
		bool _finished;
		
	private:
		
		RangeDecoder _rangeDecoder;
		ifstream* _inputFile;
		ofstream* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
		int _decodedSequenceCount;
		string _currentSeq;
		ifstream* _anchorDictFile;
		
		void decodeAnchorRead();
		kmer_type extendAnchor(kmer_type kmer, int pos, bool rightExtend);
		
		void decodeNoAnchorRead();
		void endRead();
		
		int _sequenceCount;
};

#endif /* _DNACODER_HPP_ */

