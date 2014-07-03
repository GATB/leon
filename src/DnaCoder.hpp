
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
		
		Order0Model _readTypeModel; //only 2 value in this model: with anchor or without anchor
		Order0Model _anchorAdressSizeModel;
		vector<Order0Model> _anchorAdressModel;
		Order0Model _anchorPosSizeModel;
		vector<Order0Model> _anchorPosModel;
		
		Order0Model _readSizeModel;
		Order0Model _readAnchorRevcompModel;
		vector<Order0Model> _readSizeValueModel;
		Order0Model _mutationModel;
		
		Order0Model _noAnchorReadModel;
		Order0Model _noAnchorReadSizeModel;
		vector<Order0Model> _noAnchorReadSizeValueModel;
		
		size_t _kmerSize;
		int _readSize;
		KmerModel _kmerModel;
		
		vector<int> _uniqNoOrigMutationPos;
		vector<int> _Npos;
		
		void startBlock();
		void codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right);
		void codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right);
		
		
	
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
		void writeBlock();
		void execute();
		
		void buildKmers();
		int findExistingAnchor(u_int32_t* anchorAdress);
		
		void encodeAnchorRead(int anchorPos, u_int32_t anchorAdress);
		kmer_type encodeMutations(int pos, kmer_type kmer, bool rightExtend);
		//int encodeMutations(int pos, bool rightExtend);
		int voteMutations(int pos, bool rightExtend);
		
		void encodeNoAnchorRead();
		
		
		Sequence* _sequence;
		char* _readseq;
		vector<kmer_type> _kmers;
		KmerModel::Iterator _itKmer;
		vector<int> _mutations;
		
		vector<int> _solidMutaChain;
		int _solidMutaChainPos;
		
		
		bool extendMutaChain(kmer_type kmer, int pos, bool rightExtend);
		bool extendMutaChainRec(vector< vector< vector<kmer_type> > >& mutaChains, bool rightExtend);
		
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
		
		void setup(u_int64_t blockStartPos, u_int64_t blockSize);
		void execute();
	
		string _buffer;
		bool _finished;
		
	private:
		
		RangeDecoder _rangeDecoder;
		ifstream* _inputFile;
		ofstream* _outputFile;
		u_int64_t _blockStartPos;
		u_int64_t _blockSize;
		string _currentSeq;
		
		void decodeAnchorRead();
		kmer_type decodeMutations(kmer_type kmer, int pos, bool rightExtend);
		
		void decodeNoAnchorRead();
		void endSeq();
};

#endif /* _DNACODER_HPP_ */

