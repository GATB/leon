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

#include "DnaCoder.hpp"

/*
#define PRINT_DEBUG_EXTMUTA
#define PRINT_DEBUG_ENCODER
#define PRINT_DEBUG_DECODER
*/


char bin2NTrev[4] = {'T','G','A','C'};
//char bin2NT[4] = {'A','C','T','G'};

/*
GACGCGCCGATATAACGCGCTTTCCCGGCTTTTACCACGTCGTTGAGGGCTTCCAGCGTCTCTTCGATCGGCGTGTTGTAATCCCAGCGATGAATTTG6:2308q
			Anchor pos: 8
			Anchor: GATATAACGCGCTTTCCCGGCTTTTACCACG
			§ Anchor: 8
*/

//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
AbstractDnaCoder::AbstractDnaCoder(Leon* leon) :
_kmerModel(leon->_kmerSize),
_readTypeModel(2), //only 2 value in this model: read with anchor or without anchor
//_isPrevReadAnchorableModel(2),
_noAnchorReadModel(5), _bifurcationModel(5), //5value: A, C, G, T, N
_bifurcationBinaryModel(2), //0 or 1 (lowest or highest bifurcation by alphabetical order)
_readAnchorRevcompModel(2),
_readSizeDeltaTypeModel(3),
_anchorPosDeltaTypeModel(3),
_anchorAddressDeltaTypeModel(3),
_NposDeltaTypeModel(3),
_errorPosDeltaTypeModel(3),_seqId(0)
{
	//cerr << "debug AbstractDnaCoder::AbstractDnaCoder - seg falt : not here (1) " << endl;
	_leon = leon;
	_orderReads = _leon->_orderReads;
	//_bloom = _leon->_bloom;
	_kmerSize = _leon->_kmerSize;
	
	if (_orderReads){
		for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
			//cerr << "debug AbstractDnaCoder::AbstractDnaCoder - _anchorAddressModel.push_back(Order0Model(256)) 1;" << endl;	
			_anchorKmerTypeModel.push_back(Order0Model(256));
			_anchorPosModel.push_back(Order0Model(256));
			_noAnchorReadSizeValueModel.push_back(Order0Model(256));
			_readSizeValueModel.push_back(Order0Model(256));
			//_nbReadsPerAnchorModel.push_back(Order0Model(256));
			_NposModel.push_back(Order0Model(256));
			_leftErrorPosModel.push_back(Order0Model(256));
			_numericModel.push_back(Order0Model(256));
			_leftErrorModel.push_back(Order0Model(256));
		}
	}
	else{
		for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
			//cerr << "debug AbstractDnaCoder::AbstractDnaCoder - _anchorAddressModel.push_back(Order0Model(256)) 2;" << endl;
			_anchorAddressModel.push_back(Order0Model(256));
			//_isPrevReadAnchorablePosModel.push_back(Order0Model(256));
			_anchorPosModel.push_back(Order0Model(256));
			_noAnchorReadSizeValueModel.push_back(Order0Model(256));
			_readSizeValueModel.push_back(Order0Model(256));
			_NposModel.push_back(Order0Model(256));
			_leftErrorPosModel.push_back(Order0Model(256));
			//_rightErrorPosModel.push_back(Order0Model(256));
			_numericModel.push_back(Order0Model(256));
			_leftErrorModel.push_back(Order0Model(256));
			//_rightErrorModel.push_back(Order0Model(256));
		}
	}

	//cerr << "debug AbstractDnaCoder::AbstractDnaCoder - seg falt : here (1) ?" << endl;
}

void AbstractDnaCoder::startBlock(){

	//cerr << "AbstractDnaCoder::startBlock() - begin" << endl;

	if (_orderReads){
		//cerr << "AbstractDnaCoder::startBlock() - if (_orderReads)" << endl;
		for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){

			_anchorKmerTypeModel[i].clear();
			//_isPrevReadAnchorablePosModel[i].clear();
			_anchorPosModel[i].clear();
			_noAnchorReadSizeValueModel[i].clear();
			_readSizeValueModel[i].clear();
			_NposModel[i].clear();
			_leftErrorPosModel[i].clear();
			//_rightErrorPosModel[i].clear();
			_numericModel[i].clear();
			_leftErrorModel[i].clear();
			//_rightErrorModel[i].clear();
		}
	}
	else{
		//cerr << "AbstractDnaCoder::startBlock() - if (!_orderReads)" << endl;
		for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
			//cerr << "AbstractDnaCoder::startBlock() - seg not here" << endl;
			_anchorAddressModel[i].clear();
			//cerr << "AbstractDnaCoder::startBlock() - seg here !" << endl;
			//_isPrevReadAnchorablePosModel[i].clear();
			_anchorPosModel[i].clear();
			_noAnchorReadSizeValueModel[i].clear();
			_readSizeValueModel[i].clear();
			_NposModel[i].clear();
			_leftErrorPosModel[i].clear();
			//_rightErrorPosModel[i].clear();
			_numericModel[i].clear();
			_leftErrorModel[i].clear();
			//_rightErrorModel[i].clear();
		}
	}
	//cerr << "AbstractDnaCoder::startBlock() - after .clears" << endl;

	_readTypeModel.clear();
	//_isPrevReadAnchorableModel.clear();
	_noAnchorReadModel.clear();
	_bifurcationModel.clear();
	_bifurcationBinaryModel.clear();
	_readAnchorRevcompModel.clear();
	_readSizeDeltaTypeModel.clear();
	_anchorPosDeltaTypeModel.clear();
	_anchorAddressDeltaTypeModel.clear();
	_NposDeltaTypeModel.clear();
	_errorPosDeltaTypeModel.clear();
	_prevReadSize = 0;
	_prevAnchorPos = 0;
	_prevAnchorAddress = 0;
	_prevNpos = 0;
	_prevErrorPos = 0;
	
	_processedSequenceCount = 0;

	//cerr << "AbstractDnaCoder::startBlock() - end" << endl;
}

void AbstractDnaCoder::endRead(){
	_processedSequenceCount += 1;
}
/*
void AbstractDnaCoder::codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right){
	return codeSeedNT(model, kmer, Leon::bin2nt(nt), right);
}

void AbstractDnaCoder::codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right){
	string kmerStr = kmer->toString(_kmerSize);
	//cout << kmerStr << " " << nt << " " << right << endl;
	if(right){
		kmerStr += nt;
		kmerStr.erase(kmerStr.begin());
		*kmer = model->codeSeed(kmerStr.c_str(), Data::ASCII);
		
	}
	else{		
		kmerStr.insert(kmerStr.begin(), nt);
		kmerStr.pop_back();
		*kmer = model->codeSeed(kmerStr.c_str(), Data::ASCII);
	}
		
		//cout << "\t" << kmerStr << endl;
		//cout << "\t" << kmer->toString(_kmerSize) << endl;
}
*/

void AbstractDnaCoder::codeSeedBin(KmerModel* model, kmer_type* kmer, int nt, bool right){
	//if(nt == 4) nt = 0;
	//string kmerStr = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAN"
	//kmer_type kmer2 = model->codeSeed(kmerStr.c_str(), Data::ASCII);
	//cout << kmer2->toString(_kmerSize) << endl;
	
	if(right)
	{
        /** We initialize the kmer. */
        KmerModel::Kmer tmp;  tmp.set (*kmer);

		*kmer = model->codeSeedRight (tmp, nt, Data::INTEGER).value();
	}
	else
	{
        /** We initialize the canonical kmer. */
        KmerModel::Kmer tmp;  tmp.set (revcomp(*kmer, _kmerSize));

		*kmer = model->codeSeedRight (tmp, binrev[nt], Data::INTEGER).value();
		*kmer = revcomp(*kmer, _kmerSize);
	}
}

void AbstractDnaCoder::codeSeedNT(KmerModel* model, kmer_type* kmer, char nt, bool right){
	//if(nt == 'N') nt = 'A';
	return codeSeedBin(model, kmer, Leon::nt2bin(nt), right);
}


void AbstractDnaCoder::addErrorPos(int pos, bool rightExtend){
	//if(rightExtend)
	//	_rightErrorPos.push_back(pos);
	//else
		_leftErrorPos.push_back(pos);
}

//====================================================================================
// ** DnaEncoder
//====================================================================================
DnaEncoder::DnaEncoder(Leon* leon) :
AbstractDnaCoder(leon), _itKmer(_kmerModel), _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0)//, _MCmultipleNoSolid(0)
{
	//cerr << "DnaEncoder::DnaEncoder(Leon* leon) - begin" << endl;
	cerr << "DnaEncoder::DnaEncoder - _processedSequenceCount : " << _processedSequenceCount << endl;
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);
	
	//cerr << "DnaEncoder::DnaEncoder(Leon* leon) - before \"_orderReads = _leon->_orderReads\" "<< endl;
	_orderReads = _leon->_orderReads;
	//cerr << "DnaEncoder::DnaEncoder(Leon* leon) - after \"_orderReads = _leon->_orderReads\" - orderReads =  " << _orderReads << endl;

#ifdef PRINT_DISTRIB
	_distrib.resize(maxSequences);
	_outDistrib = 0;
#endif

	//cerr << "DnaEncoder::DnaEncoder(Leon* leon) - before if(! leon->_isFasta)" << endl;

	//pour quals
	if(! leon->_isFasta)
	{
	_max_read_size = 10000;
	_nb_solids = (int *) malloc(_max_read_size * sizeof(int) );
	_qualseq = (char *) malloc(_max_read_size*sizeof(char ));
	_bufferQuals_size = Leon::READ_PER_BLOCK* 200;
	_bufferQuals = (char *) malloc(_bufferQuals_size * sizeof(char ));
	_bufferQuals_idx=0;
	
	_trunc_mode = true;
	_smoothing_threshold = 2;
		
	}
	cerr << "DnaEncoder::DnaEncoder - _processedSequenceCount : " << _processedSequenceCount << endl;
	//cerr << "DnaEncoder::DnaEncoder - end" << endl;
}

DnaEncoder::DnaEncoder(const DnaEncoder& copy) :
AbstractDnaCoder(copy._leon), _itKmer(_kmerModel),
 _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0)//, _MCmultipleNoSolid(0)
{

//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - begin" << endl;

#ifdef PRINT_DISTRIB
	_distrib.resize(maxSequences);
	_outDistrib = 0;
#endif

	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - before \"_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1)\"" << endl;
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);
	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - after \"_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1)\"" << endl;
	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - before \"_orderReads = copy._orderReads\"" << endl;
	_orderReads = copy._orderReads;
	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - after \"_orderReads = copy._orderReads\"  - orderReads =  " << _orderReads << endl;
	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - before startBlock();" << endl;
	startBlock();
	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - after startBlock();" << endl;
//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - before if(! _leon->_isFasta)" << endl;
	//for quals
	if(! _leon->_isFasta)
	{
	_max_read_size = 10000;
	_nb_solids = (int *) malloc(_max_read_size * sizeof(int) );
	_qualseq = (char *) malloc(_max_read_size*sizeof(char ));
	_bufferQuals_size = Leon::READ_PER_BLOCK* 200;
	_bufferQuals = (char *) malloc(_bufferQuals_size * sizeof(char ));
	//printf("initial buffer qual size %i \n",_bufferQuals_size );

	_bufferQuals_idx =0;
	
	_trunc_mode = true;
	_smoothing_threshold = 2;
	}
	
	///
	
	#ifdef LEON_PRINT_STAT
		_rangeEncoder1.updateModel = false;
		_rangeEncoder2.updateModel = false;
		_rangeEncoder3.updateModel = false;
		_rangeEncoder4.updateModel = false;
		_rangeEncoder5.updateModel = false;
		_rangeEncoder6.updateModel = false;
	#endif

	//cerr << "DnaEncoder::DnaEncoder(const DnaEncoder& copy) - end" << endl;
}

DnaEncoder::~DnaEncoder(){

	if(_thread_id!=0 && (_seqId+1) % Leon::READ_PER_BLOCK != 0 ){
		cerr << "DnaEncoder::~DnaEncoder() - call writeBlock()" << endl;
		writeBlock();
	}
	int nb_remaining = __sync_fetch_and_add (&_leon->_nb_thread_living, -1);

	//printf("\~ this decoder %lli seq  %lli  mctotal   %lli mltnos %p   tid %i \n",_readCount,_MCtotal,_MCmultipleNoSolid,this,_thread_id);
	__sync_fetch_and_add(&_leon->_readCount, _readCount);
	__sync_fetch_and_add(&_leon->_MCtotal, _MCtotal);
	__sync_fetch_and_add(&_leon->_readWithoutAnchorCount, _readWithoutAnchorCount);
	__sync_fetch_and_add(&_leon->_totalDnaSize, _totalDnaSize);
	__sync_fetch_and_add(&_leon->_MCuniqSolid, _MCuniqSolid);
	__sync_fetch_and_add(&_leon->_MCuniqNoSolid, _MCuniqNoSolid);
	__sync_fetch_and_add(&_leon->_MCnoAternative, _MCnoAternative);
	__sync_fetch_and_add(&_leon->_MCmultipleSolid, _MCmultipleSolid);
	//__sync_fetch_and_add(&_leon->_MCmultipleNoSolid, _MCmultipleNoSolid);
	
	#ifdef LEON_PRINT_STAT
		__sync_fetch_and_add(&_leon->_anchorAdressSize, _rangeEncoder3.getBufferSize());
		__sync_fetch_and_add(&_leon->_anchorPosSize, _rangeEncoder2.getBufferSize());
		__sync_fetch_and_add(&_leon->_readSizeSize, _rangeEncoder1.getBufferSize());
		__sync_fetch_and_add(&_leon->_bifurcationSize, _rangeEncoder4.getBufferSize());
		__sync_fetch_and_add(&_leon->_otherSize, _rangeEncoder6.getBufferSize());
		__sync_fetch_and_add(&_leon->_noAnchorSize, _rangeEncoder5.getBufferSize());

		
		_rangeEncoder1.clear();
		_rangeEncoder2.clear();
		_rangeEncoder3.clear();
		_rangeEncoder4.clear();
		_rangeEncoder5.clear();
		_rangeEncoder6.clear();
	#endif
		
#ifndef SERIAL
	//_leon->_blockwriter->incDone(1);
	_leon->_blockwriter->waitForWriter();
	
	if(nb_remaining==1)
	{
		_leon->_blockwriter->FlushWriter();
	}
	
	
	if(! _leon->_isFasta)
	{
		//pour quals
		_leon->_qualwriter->waitForWriter();
		if(nb_remaining==1)
		{
			_leon->_qualwriter->FlushWriter();
		}
	}
#endif

	//pour quals
	if(! _leon->_isFasta)
	{
		free(_nb_solids);
		free(_qualseq);
		free(_bufferQuals);
	}

	//cerr << "DnaEncoder::~DnaEncoder - _processedSequenceCount : " << _processedSequenceCount << endl;
}

void DnaEncoder::operator()(Sequence& sequence){

#ifdef PRINT_DISTRIB
	if(_sequences.size() > maxSequences){
		_sequences.pop_back();
	}
#endif

	_sequence = &sequence;
	//cout << _sequence->getIndex() << endl;
	_seqId = _sequence->getIndex() ;
	_readSize = _sequence->getDataSize();
	_readseq = _sequence->getDataBuffer();
		
	_totalDnaSize += _readSize ;
	//_lastSequenceIndex = sequence->getIndex();
	
//	if(_sequence->getIndex() % Leon::READ_PER_BLOCK == 0){

	execute();

	//_prevSequences = _sequence;
#ifdef PRINT_DISTRIB
	_sequences.insert(_sequences.begin(), _sequence);
#endif

	//_orderReads = true;//_leon->_orderReads;

	if (!_orderReads){
		if(_processedSequenceCount >= Leon::READ_PER_BLOCK ){
			writeBlock();
			startBlock();
		}
	}

}

void DnaEncoder::reset(){

	_processedSequenceCount = 0;
}

void DnaEncoder::writeBlock(){

cerr << "DnaEncoder::writeBlock() - begin" << endl;

	if(_processedSequenceCount == 0) return;
	cerr << "DnaEncoder::writeBlock() - _processedSequenceCount > 0" << endl;
	if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}
	
	int blockId = (  _seqId / Leon::READ_PER_BLOCK)   ;
	//printf("\nTid %i  WB :  blockid %i sid %llu     size: %llu  _processedSequenceCount %i\n",_thread_id, blockId, _seqId, _rangeEncoder.getBufferSize(),_processedSequenceCount );
	//cerr << "DnaEncoder::writeBlock() - _rangeEncoder.getBuffer() : " << _rangeEncoder.getBuffer() << endl;
	//_leon->_realDnaCompressedSize += _rangeEncoder.getBufferSize();
	_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), _processedSequenceCount,blockId);
	_rangeEncoder.clear();
	
	if(! _leon->_isFasta)
	{
		_leon->writeBlockLena((u_int8_t*) _bufferQuals, _bufferQuals_idx ,_processedSequenceCount, blockId);
		_bufferQuals_idx = 0;
	}
	
#ifdef PRINT_DISTRIB
	cout << "----------------------------------------------------" << endl;
	for(int i=0; i<_distrib.size(); i++){
		cout << i << "    " << _distrib[i] << endl;
	}
	cout << "Adressed:    " << _outDistrib << endl;
#endif

	/*
	cout << "------------------------" << endl;
	cout << "Read count:    " << _leon->_readCount << endl;
	cout << "Anchor kmer count:   " << _leon->_anchorKmerCount << endl;
	cout << "Read per anchor:    " << _leon->_readCount / _leon->_anchorKmerCount << endl;
	cout << "Bit per anchor: " << log2(_leon->_anchorKmerCount);
	cout << endl;
	cout << "Read without anchor: " << (_leon->_readWithoutAnchorCount*100) / _leon->_readCount << "%" << endl;
	cout << "\tWith N: " << (_leon->_noAnchor_with_N_kmer_count*100) / _leon->_readWithoutAnchorCount << "%" << endl;
	cout << "\tFull N: " << (_leon->_noAnchor_full_N_kmer_count*100) / _leon->_readWithoutAnchorCount << "%" << endl;
	cout << endl;*/
	//cout << "Total kmer: " <<  _leon->_total_kmer << endl;
	//cout << "Kmer in bloom:  " << (_leon->_total_kmer_indexed*100) / _leon->_total_kmer << "%" << endl;
	//cout << "Uniq mutated kmer: " << (_leon->_uniq_mutated_kmer*100) / _leon->_total_kmer << "%" << endl;
	//cout << "anchor Hash:   memory_usage: " << _anchorKmers.memory_usage() << "    memory_needed: " << _anchorKmers.size_entry ()*_leon->_anchorKmerCount << " (" << (_anchorKmers.size_entry ()*_leon->_anchorKmerCount*100) / _anchorKmers.memory_usage() << "%)" << endl;
	
	cerr << "DnaEncoder::writeBlock() - end" << endl;

}

void DnaEncoder::execute(){

	//if(_leon->_readCount > 18) return;
	//cout << endl << "\tEncoding seq " << _sequence->getIndex() << endl;
	//cout << "\t\t" << _readseq << endl;
	#ifdef PRINT_DEBUG_ENCODER
		cout << endl << "\tEncoding seq " << _sequence->getIndex() << endl;
		cout << "\t\t" << _readseq << endl;
	#endif
	
	//cout << _readseq << endl;
	
	_readCount +=1;
	_Npos.clear();
	
	if(_readSize < _kmerSize){
		encodeNoAnchorRead();
		endRead();
		return;
	}


	 
	//cout << _leon->_readCount << endl;
	//kmer_type anchorKmer = 0;
	u_int32_t anchorAddress;
	
	buildKmers(); // en profiter ici pour faire la compression des qual ?
	
	if(! _leon->_isFasta)
	{
		storeSolidCoverageInfo();
		smoothQuals();
	}
	

	//_isPrevReadAnchorable = false;
	int anchorPos = findExistingAnchor(&anchorAddress); //unsynch
	if(anchorPos == -1)
		anchorPos = _leon->findAndInsertAnchor(_kmers, &anchorAddress);  //unsynch

	//cout << anchorPos << endl;

	//cerr << "debug DnaEncoder::execute - orderReads = " << this->_orderReads << endl;

	if(anchorPos == -1)
		encodeNoAnchorRead();
	else{
		if (_orderReads){

			cerr << "DnaEncoder::execute() - current sequence : " << _readseq << endl;

			ofstream& unsortedReads = _leon->unsortedReads;

			/*
			//old version
			vector< list< struct ReadInfos > >& anchorsSequences = _leon->anchorsSequences;
			// create and insert ReadInfos in the wright list of the table
			// we keep only the minimal read's info needed to minimize
			// memory use
			
			struct ReadInfos* ri = new ReadInfos{};
			*/

			//ri->sequence = *_sequence;
			int readType = 0;
			//int readSize = _readSize;
			//get read

			//ri->cread = (char*)malloc((_readSize+1)*sizeof(char));
			//strncpy(ri->cread, _sequence->getDataBuffer(), _readSize);
			//ri->cread[_readSize] = '\0';

			//int anchorPos = anchorPos;
			//ri->anchorAddress = anchorAddress;
			kmer_type anchor = _kmers[anchorPos];

			//copy Npos in ri
			/*for(int i=0; i<_Npos.size(); i++){
				ri->Npos.push_back(_Npos[i]);
			}*/
			

			//ri->sread = ri->cread;

			//get the info to encode bifurcation list

			_bifurcations.clear();
			_binaryBifurcations.clear();
			_bifurcationTypes.clear();
			_leftErrorPos.clear();
			
			kmer_type kmer = anchor;

			for(int i=anchorPos-1; i>=0; i--){
				kmer = buildBifurcationList(i, kmer, false);
			}

			kmer = anchor;
			for(int i=anchorPos+_kmerSize; i<_readSize; i++){
				kmer = buildBifurcationList(i, kmer, true);
			}

			bool isRevComp;
			if(anchor == min(anchor, revcomp(anchor, _kmerSize))){
				isRevComp = 0;
				//cerr << "debug DnaEncoder::encodeAnchorRead - revcomp : " << 0 << endl;
			}
			else{
				anchor = revcomp(anchor, _kmerSize);
				isRevComp = 1;
				//cerr << "debug DnaEncoder::encodeAnchorRead - revcomp : " << 1 << endl;
			}
			//cerr << "DnaEncoder::execute() - isRevComp : " << isRevComp << endl;
			//copy bifurcationTypes
			//cerr << "\ndebug DnaEncoder::execute() - bifurcations : " << endl;
			/*for(int i=0; i<_bifurcationTypes.size(); i++){
				ri->bifurcationTypes.push_back(_bifurcationTypes[i]);
			//	cerr << "\ndebug DnaEncoder::execute() - _bifurcationTypes[i] : " << (int) _bifurcationTypes[i] << endl;
			}*/

			//copy bifurcations
			//cerr << "\ndebug DnaEncoder::execute() - bifurcations : " << endl;
			/*for(int i=0; i<_bifurcations.size(); i++){
				ri->bifurcations.push_back(_bifurcations[i]);
				//cerr << "\ndebug DnaEncoder::execute() - _bifurcations[i] : " << (int) _bifurcations[i] << endl;
			}*/

			//copy binaryBifurcations
			//cerr << "\ndebug DnaEncoder::execute() - bifurcations : " << endl;
			/*for(int i=0; i<_binaryBifurcations.size(); i++){
				ri->binaryBifurcations.push_back(_binaryBifurcations[i]);
				//cerr << "\ndebug DnaEncoder::execute() - _binaryBifurcations[i] : " << (int) _binaryBifurcations[i] << endl;
			}*/

			//copy leftErrorPos
			/*for(int i=0; i<_leftErrorPos.size(); i++){
				ri->leftErrorPos.push_back(_leftErrorPos[i]);
			}*/

			//ri->revcomp =
			//ri->revAnchor =
			//ri->revcomp =
			//ri->NposCount =
			//ri->nbLeftError =

			//insert the read in the ordered table

			//cerr << "\tdebug DnaEncoder::execute - anchorAddress = " << anchorAddress << endl;

			/*if (anchorAddress >= anchorsSequences.size()){
				anchorsSequences.resize(anchorsSequences.size() + _leon->_coverage);
			}*/
			//anchorsSequences[anchorAddress].push_back(*ri);


			//write on the file to sort	

			unsortedReads << anchor << ";" << 
							isRevComp << ":" <<
							_readSize << ":" <<
							anchorPos << ":" /*<<
							anchorAddress << ":"*/;
							
							/*
							<< ":" <<
							int revcomp not necessary, determined at compression
							*/
			//cerr << "\tdebug DnaEncoder::execute - test unsortedReads = " << endl;
			//Save N positions
			unsortedReads << _Npos.size() << ":";
			int i = 0;
			for(; i<((int)_Npos.size()-1); ++i){
				//cerr << "\tdebug DnaEncoder::execute - Npos.size() = " << Npos.size() << endl;
			 	unsortedReads << _Npos[i] << ",";
			}
			if (i != 0){
				unsortedReads << _Npos[i];
			}
			unsortedReads << ":";
		
			//Save left errors
			unsortedReads << _leftErrorPos.size() << ":";
			i = 0;
			for(; i< ((int)_leftErrorPos.size()-1); ++i){
				unsortedReads << _leftErrorPos[i] << ",";
				//cerr << "debug DnaEncoder::encodeReadsInfos - _prevErrorPos : " << _leftErrorPos[i] << endl;
			}
			if (i != 0){
				unsortedReads << _leftErrorPos[i];
			}
			unsortedReads << ":";
			
			
			//Save bifurcation types
			u_int64_t bifType0 = 0;
			u_int64_t bifType1 = 0;
			//cerr << "\ndebug DnaEncoder::encodeReadsInfos - bifurcations : \n"<< endl <<
			//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< endl <<
			//"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			for(int i=0; i<_bifurcationTypes.size(); i++){

				u_int8_t type = _bifurcationTypes[i];
				//to save some place, type is merged with bif type here
				//unsortedReads << (int) type;
				//cerr << "debug DnaEncoder::encodeReadsInfos -  _bifurcationTypes[i] : " <<  (int) _bifurcationTypes[i] << endl;

				if(type == 0){
					unsortedReads << (int) _bifurcations[bifType0];
					//cerr << "debug DnaEncoder::encodeReadsInfos - _bifurcations[bifType0] : " << (int) _bifurcations[bifType0] << endl;
					bifType0 += 1;
				}
				else{

					unsortedReads << (int) (4+_binaryBifurcations[bifType1]);
					//cerr << "debug DnaEncoder::encodeReadsInfos - _binaryBifurcations[bifType1] : " << (int) _binaryBifurcations[bifType1] << endl;
					bifType1 += 1;
				}
			}
			unsortedReads << endl;
			//unsortedReads << ri->cread << endl;
			//cerr << "\t DnaEncoder::encodeAnchorRead - _processedSequenceCount : " << _processedSequenceCount << endl;

		}
		else{
		encodeAnchorRead(anchorPos, anchorAddress);
		}
	}
	//}
	
	if (!_orderReads){
		endRead();
	}

}




double DnaEncoder::char2proba(char c)
{
	int phred = c -33;
	
	double proba =  exp(-phred* log(10)/10);
	return proba;
	//Q 10 : 10% err
	//Q 20 : 1% err
	//Q 30  0.1% err ..
}


char DnaEncoder::char2phred(char c)
{
	return c -33;
}


void DnaEncoder::smoothQuals()
{
	strcpy (_qualseq, _sequence->getQuality().c_str());  // copy the qual sequence of this read in _qualseq

	if(! _leon->_lossless)
	{
		for (int ii=0; ii< _readSize; ii++)
		{
			if ((_nb_solids[ii]>= _smoothing_threshold) || (((int) _qualseq[ii] > (int) '@') && _trunc_mode ))
			{
				apply_smoothing_at_pos (ii);
			}
		}
	}
	
	_qualseq[_readSize]='\n';
	_qualseq[_readSize+1]='\0';
	
	if( (_bufferQuals_idx+ _readSize+1 ) >= _bufferQuals_size)
	{
		//printf("_bufferQuals_size %i  _bufferQuals_idx %i  seqid %zu \n",_bufferQuals_size,_bufferQuals_idx,_sequence->getIndex()  );
		_bufferQuals_size = _bufferQuals_size * 2;
		_bufferQuals = (char *) realloc(_bufferQuals,_bufferQuals_size * sizeof(char) );
	}
	
	strcpy(_bufferQuals + _bufferQuals_idx, _qualseq);

	_bufferQuals_idx += _readSize+1 ; // with last '\n'

	//fprintf(_leon->_testQual,"%s",_qualseq); //st_qualseq.c_str()

}



bool DnaEncoder::apply_smoothing_at_pos(int pos)
{
	if(char2phred(_qualseq[pos])==0 || char2phred(_qualseq[pos])==2 )
		return false;
	
	bool ok_to_smooth= true;
	
	int diff = ('@' -  _qualseq[pos]);
	if(  diff > 10   )
	{
		if(_nb_solids[pos]>(diff-5))
			ok_to_smooth =true;
		else
			ok_to_smooth = false;
	}
	
	if(ok_to_smooth)
	{
		_qualseq[pos] = '@'; //smooth qual
		return true;
	}
	else return false;
	
}



void DnaEncoder::storeSolidCoverageInfo()
{
	kmer_type kmer, kmerMin;
	Node node;
	
	if(_readSize >= _max_read_size)
	{
		_max_read_size = _readSize + 1000;
		_nb_solids = (int *) realloc(_nb_solids,_max_read_size * sizeof(int) );
		_qualseq = (char *) realloc(_qualseq,_max_read_size*sizeof(char ));

	}
	memset(_nb_solids,0,_max_read_size * sizeof(int) );

	for(int ii=0; ii<_kmers.size(); ii++){
		kmer = _kmers[ii];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		node = Node(Node::Value(kmerMin));

		if(_leon->_graph.contains(node))
		{
			//increments all pos covered by the solid kmer
			for (int jj=0; jj< _kmerSize ; jj++)
			{
				_nb_solids[ii+jj] ++ ;
			}
		}
	}
}


void DnaEncoder::buildKmers(){

	
	
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] == 'N'){
			_Npos.push_back(i);
			_readseq[i] = 'A';
		}
	}
		
	/*
	bool lala = false;
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] == 'N'){
			lala = true;
			break;
		}
	}
	if(lala){
		cout << "lala" << endl;
		cout << string(_readseq).substr(0, _readSize) << endl;
		
		for(int i=0; i<_readSize; i++){
			if(_readseq[i] == 'N'){
				_Npos.push_back(i);
				_readseq[i] = 'A';
				lala = true;
			}
		}
		
		_itKmer.setData(_sequence->getData());
		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			cout << (*_itKmer).toString(_kmerSize) << endl;
		}
		

	
	}
	*/
	_itKmer.setData(_sequence->getData());
	
	_kmers.clear();
	for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
		//cout << (*_itKmer).toString(_kmerSize) << endl;
		_kmers.push_back(_itKmer->value());
	}
	

#ifdef PRINT_DISTRIB
	unordered_set<u_int64_t> H;
	for(kmer_type kmer : _kmers){
		kmer_type kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		H.insert(kmerMin.getVal());
	}

	for(int i=0; i<_sequences.size(); i++){
		Sequence* sequence = _sequences[i];

		_itKmer.setData(sequence->getData());
		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			kmer_type kmerMin2 = min(_itKmer->value(), revcomp(_itKmer->value(), _kmerSize));
			if(H.find(kmerMin2.getVal()) != H.end()){
				_distrib[i] += 1;
				return;
			}
		}
	}

	_outDistrib += 1;
#endif

	//_sequences

	//if(_sequence->getIndex() == 53445) cout << _Npos.size() << endl;
	/*
	if(_sequence->getIndex() == 29){
		cout << string(_readseq).substr(0, _readSize) << endl;
		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			cout << (*_itKmer).toString(_kmerSize) << endl;
		}
	}*/
	
}

int DnaEncoder::findExistingAnchor(u_int32_t* anchorAddress){

	kmer_type kmer, kmerMin;

	/*
	if(_prevSequences != NULL){

		unordered_map<u_int64_t, u_int64_t> H;
		for(int i=0; i<_kmers.size(); i++){
			kmer = _kmers[i];
			kmerMin = min(kmer, revcomp(kmer, _kmerSize));
			H[kmerMin.getVal()] = i;
		}

		_itKmer.setData(_prevSequences->getData());
		int i = 0;
		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			kmer = _itKmer->value();
			kmerMin = min(kmer, revcomp(kmer, _kmerSize));
			if(H.find(kmerMin.getVal()) != H.end()){
				_isPrevReadAnchorable = true;
				_isPrevReadAnchorablePos = i;
				return H[kmerMin.getVal()];
			}
			i += 1;
		}

	}*/

	
	for(int i=0; i<_kmers.size(); i++){
		kmer = _kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		//cout << "kmer : " << kmerMin << endl;
		if(_leon->anchorExist(kmerMin, anchorAddress)){
			return i;
		}
	}
	return -1;
}

bool DnaEncoder::isReadAnchorable(){
	int nbKmerSolid = 0;
	kmer_type kmer, kmerMin;
	Node node;

	for(int i=0; i<_kmers.size(); i++){

		kmer = _kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		node = Node(Node::Value(kmerMin));

		if(_leon->_graph.contains(node)){
			nbKmerSolid += 1;
			i += _kmerSize;
		}

		if(nbKmerSolid >= 2) return true;
	}

	return nbKmerSolid >= 2;

}


//old version of ordering reads
/*
void DnaEncoder::encodeReadsInfos(vector< list< struct ReadInfos > > anchorsSequences){

	for (int index=0; index<anchorsSequences.size(); index++){

		list<struct ReadInfos> riList = anchorsSequences[index];
		for (list<struct ReadInfos>::iterator ri=riList.begin(); ri != riList.end(); ++ri){
   			
   			kmer_type anchor = ri->anchor;
   			int anchorPos = ri->anchorPos;
   			u_int64_t anchorAddress = ri->anchorAddress;


			_rangeEncoder.encode(_readTypeModel, 0);
			//cerr << "debug DnaEncoder::encodeReadsInfos - _rangeEncoder.encode(_readTypeModel, 0);" << endl;
			
			//Encode read size x
			CompressionUtils::encodeNumeric(_rangeEncoder, _readSizeValueModel, _readSize);
			//cerr << "debug DnaEncoder::encodeReadsInfos - encodeNumeric(_rangeEncoder, _readSizeValueModel, _readSize)" << endl;
	
			//Encode anchor pos x
			CompressionUtils::encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);
			//cerr << "debug DnaEncoder::encodeReadsInfos - encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);" << endl;

			//Encode anchor address x
			CompressionUtils::encodeNumeric(_rangeEncoder, _anchorAddressModel, anchorAddress);
			
			//Encode a bit that says if the anchor is normal or revcomp
			if(anchor == min(anchor, revcomp(anchor, _kmerSize))){
				_rangeEncoder.encode(_readAnchorRevcompModel, 0);
				//cerr << "debug DnaEncoder::encodeReadsInfos - revcomp : " << 0 << endl;
			}
			else{
				_rangeEncoder.encode(_readAnchorRevcompModel, 1);
				//cerr << "debug DnaEncoder::encodeReadsInfos - revcomp : " << 1 << endl;
			}


			//Do this before the sort
			//just keep the info to compress
			
			//_bifurcations.clear();
			//_binaryBifurcations.clear();
			//_bifurcationTypes.clear();
			//_leftErrorPos.clear();
			
			//kmer_type kmer = anchor;
			//for(int i=anchorPos-1; i>=0; i--){
			//	kmer = buildBifurcationList(i, kmer, false);
			//}

			//kmer = anchor;
			//for(int i=anchorPos+_kmerSize; i<_readSize; i++){
			//	kmer = buildBifurcationList(i, kmer, true);
			//}
				
			vector<int> _Npos = ri->Npos;
			vector<int> leftErrorPos = ri->leftErrorPos;
			vector<u_int8_t> bifurcations = ri->bifurcations;
			vector<u_int8_t> binaryBifurcations = ri->binaryBifurcations;
			vector<u_int8_t> bifurcationTypes = ri->bifurcationTypes;

			//Encode N positions
			_prevNpos = 0;
			CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _Npos.size());
			for(int i=0; i<_Npos.size(); i++){
				CompressionUtils::encodeNumeric(_rangeEncoder, _NposModel, _Npos[i]-_prevNpos);
				_prevNpos = _Npos[i];
				//cerr << "debug DnaEncoder::encodeReadsInfos - _prevNpos : " << _Npos[i] << endl;
			}
			
			//Encode left errors
			CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorModel, _leftErrorPos.size());
			sort(_leftErrorPos.begin(), _leftErrorPos.end());
			_prevErrorPos = 0;
			for(int i=0; i<_leftErrorPos.size(); i++){
				CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorPosModel, _leftErrorPos[i]-_prevErrorPos);
				_prevErrorPos = _leftErrorPos[i];
				//cerr << "debug DnaEncoder::encodeReadsInfos - _prevErrorPos : " << _leftErrorPos[i] << endl;
			}
			
			//encode bifurcation types
			u_int64_t bifType0 = 0;
			u_int64_t bifType1 = 0;
			//cerr << "\ndebug DnaEncoder::encodeReadsInfos - bifurcations : " << endl;
			for(int i=0; i<bifurcationTypes.size(); i++){
				//cerr << "lol1" << endl;
				u_int8_t type = bifurcationTypes[i];
				//cerr << "lol2" << endl;
				if(type == 0){
					_rangeEncoder.encode(_bifurcationModel, bifurcations[bifType0]);
					bifType0 += 1;
			//		cerr << "debug DnaEncoder::encodeReadsInfos - bifType0 : " << bifType0 << endl;
				}
				else{
					_rangeEncoder.encode(_bifurcationBinaryModel, binaryBifurcations[bifType1]);
					bifType1 += 1;
			//		cerr << "debug DnaEncoder::encodeReadsInfos - bifType1 : " << bifType1 << endl;
				}
			}
			
			//endRead();

			
   		}
	}

	writeBlock();
	startBlock();

}
*/


void DnaEncoder::encodeAnchorRead(int anchorPos, u_int32_t anchorAddress){
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\tEncode anchor read" << endl;
	#endif
	//printf("encode  anchor read \n");

	//encode read type (0: read with anchor, 1: read without anchor)
	#ifdef LEON_PRINT_STAT
		_rangeEncoder6.encode(_readTypeModel, 0);
	#endif
	_rangeEncoder.encode(_readTypeModel, 0);
	//cerr << "debug DnaEncoder::encodeAnchorRead - _rangeEncoder.encode(_readTypeModel, 0);" << endl;
	
	u_int64_t deltaValue;
	u_int8_t deltaType;
	
	//Encode read size
	//deltaType = CompressionUtils::getDeltaValue(_readSize, _prevReadSize, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder1.encode(_readSizeDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder1, _readSizeValueModel, _readSize);
	#endif
	//_rangeEncoder.encode(_readSizeDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _readSizeValueModel, _readSize);
	//cerr << "debug DnaEncoder::encodeAnchorRead - encodeNumeric(_rangeEncoder, _readSizeValueModel, _readSize)" << endl;
	//_prevReadSize = _readSize;
	//printf("read size %i  deltaValue %i\n",_readSize,deltaValue);

	//Encode anchor pos
	//deltaType = CompressionUtils::getDeltaValue(anchorPos, _prevAnchorPos, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder2.encode(_anchorPosDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder2, _anchorPosModel, anchorPos);
	#endif
	//_rangeEncoder.encode(_anchorPosDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);
	//cerr << "debug DnaEncoder::encodeAnchorRead - encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);" << endl;
	//_prevAnchorPos = anchorPos;
	//printf("anchor pos %i \n",anchorPos);

	//Encode anchor address
	//deltaType = CompressionUtils::getDeltaValue(anchorAddress, _prevAnchorAddress, &deltaValue);
	#ifdef LEON_PRINT_STAT
		//_rangeEncoder3.encode(_anchorAddressDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder3, _anchorAddressModel, anchorAddress);
	#endif
	//_rangeEncoder.encode(_anchorAddressDeltaTypeModel, deltaType);
	//if(_isPrevReadAnchorable){
		//	_rangeEncoder.encode(_isPrevReadAnchorableModel, 0);
		//	CompressionUtils::encodeNumeric(_rangeEncoder, _isPrevReadAnchorablePosModel, _isPrevReadAnchorablePos);
		//}
		//else{
		//_rangeEncoder.encode(_isPrevReadAnchorableModel, 1);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorAddressModel, anchorAddress);
		//CompressionUtils::encodeNumeric(_rangeEncoder, _isPrevReadAnchorablePosModel, _isPrevReadAnchorablePos);
		//}
	//_prevAnchorAddress = anchorAddress;
	//printf("anchor adress %i \n",anchorAddress);

	
	kmer_type anchor = _kmers[anchorPos];
	
	//Encode a bit that says if the anchor is normal or revcomp
	if(anchor == min(anchor, revcomp(anchor, _kmerSize))){
		#ifdef LEON_PRINT_STAT
			_rangeEncoder6.encode(_readAnchorRevcompModel, 0);
		#endif
		_rangeEncoder.encode(_readAnchorRevcompModel, 0);
		//cerr << "debug DnaEncoder::encodeAnchorRead - revcomp : " << 0 << endl;
	}
	else{
		#ifdef LEON_PRINT_STAT
			_rangeEncoder6.encode(_readAnchorRevcompModel, 1);
		#endif
		_rangeEncoder.encode(_readAnchorRevcompModel, 1);
		//cerr << "debug DnaEncoder::encodeAnchorRead - revcomp : " << 1 << endl;
	}

	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tAnchor pos: " << anchorPos << endl;
		cout << "\t\t\tAnchor: " << _kmers[anchorPos].toString(_kmerSize) << endl;
	#endif

	_bifurcations.clear();
	_binaryBifurcations.clear();
	_bifurcationTypes.clear();
	_leftErrorPos.clear();
	//_rightErrorPos.clear();

	
	kmer_type kmer = anchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = buildBifurcationList(i, kmer, false);
		//i = buildBifurcationList(i, false);
		//cout << kmer.toString(_kmerSize) << endl;
	}

	kmer = anchor;
	for(int i=anchorPos+_kmerSize; i<_readSize; i++){
		//cout << "Pos: " << i << endl;
		kmer = buildBifurcationList(i, kmer, true);
		//i = buildBifurcationList(i, true);
	//for(int i=anchorPos; i<_kmers.size()-1; i++)
		//cout << kmer.toString(_kmerSize) << endl;
	}
		

	//Encode N positions
	_prevNpos = 0;
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder6, _numericModel, _Npos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _Npos.size());
	for(int i=0; i<_Npos.size(); i++){
		//deltaType = CompressionUtils::getDeltaValue(_Npos[i], _prevNpos, &deltaValue);
		//_rangeEncoder.encode(_NposDeltaTypeModel, deltaType);
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder6, _NposModel, _Npos[i]-_prevNpos);
		#endif
		CompressionUtils::encodeNumeric(_rangeEncoder, _NposModel, _Npos[i]-_prevNpos);
		_prevNpos = _Npos[i];
		//cerr << "debug DnaEncoder::encodeAnchorRead - _prevNpos : " << _Npos[i] << endl;
	}
	
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder6, _leftErrorModel, _leftErrorPos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorModel, _leftErrorPos.size());
	sort(_leftErrorPos.begin(), _leftErrorPos.end());
	_prevErrorPos = 0;
	for(int i=0; i<_leftErrorPos.size(); i++){
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder6, _leftErrorPosModel, _leftErrorPos[i]-_prevErrorPos);
		#endif
		CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorPosModel, _leftErrorPos[i]-_prevErrorPos);
		_prevErrorPos = _leftErrorPos[i];
		//cerr << "debug DnaEncoder::encodeAnchorRead - _prevErrorPos : " << _leftErrorPos[i] << endl;
	}

	/*
	//Encode the positions of sequencing errors
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder4, _leftErrorModel, _leftErrorPos.size());
		CompressionUtils::encodeNumeric(_rangeEncoder4, _rightErrorModel, _rightErrorPos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorModel, _leftErrorPos.size());
	CompressionUtils::encodeNumeric(_rangeEncoder, _rightErrorModel, _rightErrorPos.size());
	_prevErrorPos = anchorPos-1;
	for(int i=0; i<_leftErrorPos.size(); i++){
		u_int64_t errorPos = _leftErrorPos[i];
		CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorPosModel, _prevErrorPos-errorPos);
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder4, _leftErrorPosModel, _prevErrorPos-errorPos);
		#endif
		_prevErrorPos = errorPos;
	}
	_prevErrorPos = anchorPos+_kmerSize;
	for(int i=0; i<_rightErrorPos.size(); i++){
		u_int64_t errorPos = _rightErrorPos[i];
		#ifdef LEON_PRINT_STAT
			CompressionUtils::encodeNumeric(_rangeEncoder4, _rightErrorPosModel, errorPos-_prevErrorPos);
		#endif
		CompressionUtils::encodeNumeric(_rangeEncoder, _rightErrorPosModel, errorPos-_prevErrorPos);
		_prevErrorPos = errorPos;
	}*/

	/*
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder4, _numericSizeModel, _numericModel, _leftErrorPos.size()+_rightErrorPos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _leftErrorPos.size()+_rightErrorPos.size());
	for(int i=0; i<_errorPos.size(); i++){
		deltaType = CompressionUtils::getDeltaValue(_errorPos[i], _prevErrorPos, &deltaValue);
		#ifdef LEON_PRINT_STAT
			_rangeEncoder4.encode(_errorPosDeltaTypeModel, deltaType);
			CompressionUtils::encodeNumeric(_rangeEncoder4, _errorPosSizeModel, _errorPosModel, deltaValue);
		#endif
		_rangeEncoder.encode(_errorPosDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder, _errorPosModel, deltaValue);
		_prevErrorPos = _errorPos[i];
	}*/
	
	u_int64_t bifType0 = 0;
	u_int64_t bifType1 = 0;
	//cerr << "\ndebug DnaEncoder::encodeAnchorRead - bifurcations : " << endl;
	//cout << _bifurcationTypes.size() << " " << _bifurcations.size() << " " << _binaryBifurcations.size() << endl;
	for(int i=0; i<_bifurcationTypes.size(); i++){
		u_int8_t type = _bifurcationTypes[i];
		if(type == 0){
			#ifdef LEON_PRINT_STAT
				_rangeEncoder4.encode(_bifurcationModel, _bifurcations[bifType0]);
			#endif
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationModel, _bifurcations[bifType0]);
			bifType0 += 1;
			//cerr << "debug DnaEncoder::encodeAnchorRead - bifType0 : " << bifType0 << endl;
		}
		else{
			#ifdef LEON_PRINT_STAT
				_rangeEncoder4.encode(_bifurcationBinaryModel, _binaryBifurcations[bifType1]);
			#endif
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationBinaryModel, _binaryBifurcations[bifType1]);
			bifType1 += 1;
			//cerr << "debug DnaEncoder::encodeAnchorRead - bifType1 : " << bifType1 << endl;
		}
	}

	
	
}

void DnaEncoder::encodeSortedFileAnchor(kmer_type anchor){

	u_int64_t anchor_uint64t = anchor.getVal();
	cerr << "\tDnaEncoder::encodeSortedFileAnchor(kmer_type anchor) - anchor_uint64t : " << anchor_uint64t << endl;
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorKmerTypeModel, anchor_uint64t);
	//endRead();
	//cerr << "\tDnaEncoder::encodeSortedFileAnchor(kmer_type anchor) - _processedSequenceCount : " << _processedSequenceCount << endl;
}

void DnaEncoder::encodeSortedFileRead(kmer_type anchor, int isRevComp, int readSize, int anchorPos, /*int anchorAddress,*/ vector<int> Npos,
									vector<int> leftErrorPos, vector<u_int8_t> bifurcations, 
									vector<u_int8_t> binaryBifurcations, vector<u_int8_t> bifurcationTypes){

	//readTypeModel is not needed anymore for sorted file
	//_rangeEncoder.encode(_readTypeModel, 0);

	CompressionUtils::encodeNumeric(_rangeEncoder, _readSizeValueModel, readSize);
	cerr << "\tDnaEncoder::encodeSortedFileRead - readSize : " << readSize << endl;
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorPosModel, anchorPos);
	cerr << "\tDnaEncoder::encodeSortedFileRead - anchorPos : " << anchorPos << endl;

	_rangeEncoder.encode(_readAnchorRevcompModel, isRevComp);

	int prevNpos = 0;
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, Npos.size());
	for(int i=0; i<Npos.size(); i++){
		CompressionUtils::encodeNumeric(_rangeEncoder, _NposModel, Npos[i]-prevNpos);
		prevNpos = Npos[i];
		//cerr << "\tDnaEncoder::encodeSortedFileRead - prevNpos : " << Npos[i] << endl;
	}

	CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorModel, leftErrorPos.size());
	sort(leftErrorPos.begin(), leftErrorPos.end());
	int prevErrorPos = 0;
	for(int i=0; i<leftErrorPos.size(); i++){
		CompressionUtils::encodeNumeric(_rangeEncoder, _leftErrorPosModel, leftErrorPos[i]-prevErrorPos);
		prevErrorPos = leftErrorPos[i];
		//cerr << "\tDnaEncoder::encodeSortedFileRead - prevErrorPos : " << leftErrorPos[i] << endl;
	}

	u_int64_t bifType0 = 0;
	u_int64_t bifType1 = 0;
	//cerr << "\tDnaEncoder::encodeSortedFileRead - bifurcations : " << endl;
	cerr << "\n\tDnaEncoder::encodeSortedFileRead - nb bifurcations : " << bifurcationTypes.size() << " nb normal bif : " << bifurcations.size() << " nb binary bif : " << binaryBifurcations.size() << endl;
	for(int i=0; i<bifurcationTypes.size(); i++){
		u_int8_t type = bifurcationTypes[i];
		if(type == 0){
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationModel, bifurcations[bifType0]);
			bifType0 += 1;
			//cerr << "\tDnaEncoder::encodeSortedFileRead - bifType0 : " << bifType0 << endl;
		}
		else{
			//cout << Leon::nt2bin(_bifurcations[i]) << " ";
			_rangeEncoder.encode(_bifurcationBinaryModel, binaryBifurcations[bifType1]);
			bifType1 += 1;
			//cerr << "\tDnaEncoder::encodeSortedFileRead - bifType1 : " << bifType1 << endl;
		}
	}

	cerr << "\tDnaEncoder::encodeSortedFileRead - buffer size : " << _rangeEncoder.getBufferSize() << endl;
	endRead();
	cerr << "\tDnaEncoder::encodeSortedFileRead - _processedSequenceCount : " << _processedSequenceCount << endl;
}

void DnaEncoder::encodeSortedFileWriteBlock(){
	cerr << "DnaEncoder::operator() - call writeBlock()" << endl;
	writeBlock();
	//startBlock();
}
	
kmer_type DnaEncoder::buildBifurcationList(int pos, kmer_type kmer, bool rightExtend){
		
	char nextNt = _readseq[pos];
	int nextNtBin = Leon::nt2bin(nextNt);

	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;
		//return pos;
	}
	
	//kmer_type kmerMin;
	kmer_type uniqKmer;
	bool firstSolidKmer = false;
	int uniqNt;
	//u_int8_t binNt2;
	bool isKmerSolid = false;
	
	//kmer_type kmer;

	//if(rightExtend)
	//	kmer = _kmers[pos-_kmerSize];
	//else
	//	kmer = _kmers[pos+1];
	
	//cout << "Kmer: " << kmer.toString(_kmerSize) << endl;
	
	
	
	
	//cout << "Real next nt: " << nextNt << endl;
	
	//cout << nextNt << endl;
	
	//kmer = _kmers[pos];
	
	//if(kmer.toString(_kmerSize).find("N") !=  string::npos) cout << "lala" << endl;

	int indexedKmerCount = 0;
	
	
	
	std::bitset<4> res4;
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		Node node;
		//kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		kmer_type mutatedKmerMin = mutatedKmer;
		node = Node(Node::Value(mutatedKmerMin));
		
		//mutatedKmer.printASCII(_kmerSize);
		
		res4[nt] = _leon->_graph.contains(node);
	}
	
	
	
//	std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);
	for(int nt=0; nt<4; nt++){
		
		
		if(res4[nt]){

			indexedKmerCount += 1;

			if(!firstSolidKmer){
				firstSolidKmer = true;
				uniqKmer = kmer;
				codeSeedBin(&_kmerModel, &uniqKmer, nt, rightExtend);
			}

			
			if(nt == nextNtBin){
				isKmerSolid = true;
			}
		}
		
	}

	
	
	
	
	
	_MCtotal +=1;
	

	if(isKmerSolid){

		if(indexedKmerCount == 1){
			_MCuniqSolid += 1;
			return uniqKmer;
		}
		else if(indexedKmerCount == 2){

			char nt1 = -1;
			char nt2 = -1;

			for(int nt=0; nt<4; nt++){
				if(res4[nt]){
					//cout << "\t" << nt << endl;
					if(nt1 == -1)
						nt1 = nt;
					else if(nt2 == -1)
						nt2 = nt;
					else break;
				}
			}


			if(nt1 == nextNtBin){
				//cout << "\t0" << endl;
				_binaryBifurcations.push_back(0);
				_bifurcationTypes.push_back(1);
				_MCmultipleSolid += 1;
			}
			else if(nt2 == nextNtBin){
				//cout << "\t1" << endl;
				_binaryBifurcations.push_back(1);
				_bifurcationTypes.push_back(1);
				_MCmultipleSolid += 1;
			}
			else{

				//if(_sequence->getIndex() < 20)
				//	cout << "\tallo" << endl;
				//_MCuniqNoSolid += 1;
				//nextNt = Leon::bin2nt(nt1);
				//_bifurcations.push_back(nextNtBin);
				//_errorPos.push_back(pos);
				//_bifurcationTypes.push_back(0);
				//return uniqKmer;

				/*
				_MCuniqNoSolid += 1;
				//nextNt = Leon::bin2nt(nt1);

				_bifurcationTypes.push_back(0);
				_bifurcations.push_back(nextNtBin);
				_errorPos.push_back(pos);

				nextNtBin = getBestPath(pos, kmer, res4, rightExtend);
				//cout << (int)nextNtBin << endl;
				if(nextNtBin == -1){
					nextNtBin = nt1;
				}
				nextNt = Leon::bin2nt(nextNtBin);
				_bifurcations.push_back(nextNtBin);
				_bifurcationTypes.push_back(0);*/

			}
			//cout << "PROBLEME IN BUILD BINARY BIFURCATION (DnaEncoder - buildBifurcationList)" << endl;

			//if(_sequence->getIndex() < 10)
			//	cout << (char) Leon::bin2nt(nextNtBin) << endl;;

			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;



		}
		else{

			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;
		}

	}
	else{


		if(indexedKmerCount == 0){
			//cout << "PAF_BREAK    " << pos << endl;
			_MCnoAternative += 1;
		}
		else if(indexedKmerCount == 1){
			//cout << "PAF_UNIQ    " << pos << endl;
			_MCuniqNoSolid += 1;

			//_leon->_readWithAnchorMutationChoicesSize += 0.25;

			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			addErrorPos(pos, rightExtend);
			//_errorPos.push_back(pos);
			return uniqKmer;
		}
		else if(indexedKmerCount == 2){
			//cout << "PAF_MULTIPLE    " << pos << endl;
			_MCuniqNoSolid += 1;

			//nextNt = Leon::bin2nt(nt1);
			//_bifurcations.push_back(nextNtBin);
			/*
			int memo = nextNtBin;
			nextNtBin = getBestPath(pos, kmer, res4, rightExtend);
			if(nextNtBin == -1){
				nextNtBin = memo;
				nextNt = Leon::bin2nt(nextNtBin);
			}
			else{
				//_errorPos.push_back(pos);
				//_bifurcations.push_back(nextNtBin);
				//_bifurcationTypes.push_back(0);
				nextNt = Leon::bin2nt(nextNtBin);
			}*/

			//encode error
			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			addErrorPos(pos, rightExtend);

			//get the first path in bufurcation
			for(int nt=0; nt<4; nt++){
				if(res4[nt]){
					nextNtBin = nt;
					nextNt = Leon::bin2nt(nt);
					break;
				}
			}

			codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
			return kmer;
			//return uniqKmer;
			//_bifurcationTypes.push_back(0);

			/*_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);

			return uniqKmer;*/
		}
		else{
			_MCuniqNoSolid += 1;
			//cout << "PAF_MULTIPLE    " << pos << endl;
			/*
			_errorPos.push_back(pos);
			_bifurcations.push_back(nextNtBin);
			_bifurcationTypes.push_back(0);
			_errorPos.push_back(pos);
			return uniqKmer;*/

		}

		/*
		if(indexedKmerCount != 0){
			int memo = nextNtBin;
			nextNtBin = getBestPath(pos, kmer, res4, rightExtend);
			if(nextNtBin == -1){
				nextNtBin = memo;
				nextNt = Leon::bin2nt(nextNtBin);
			}
			else{
				_errorPos.push_back(pos);
				_bifurcations.push_back(nextNtBin);
				_bifurcationTypes.push_back(0);
				nextNt = Leon::bin2nt(nextNtBin);
			}

			//if(indexedKmerCount == 1 && memo != nextNtBin){
			//	cout << "ALLO" << endl;
			//}
		}*/


		//_leon->_readWithAnchorMutationChoicesSize += 0.25;
		_bifurcations.push_back(nextNtBin);
		_bifurcationTypes.push_back(0);
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;

	}



	/*

	if(indexedKmerCount == 2){

		//_bifurcations.push_back(nextNtBin);

		//cout << "Bin bifurc " << (int)nextNtBin << endl;



	}
	else{
		if(indexedKmerCount == 0){
			_MCnoAternative += 1;
		}
		else{
			if(isKmerSolid){
				_MCmultipleSolid += 1;
				
			}
			else{
				//_MCmultipleNoSolid += 1;


				int memo = nextNtBin;
				nextNtBin = getBestPath(pos, kmer, res4, rightExtend);
				if(nextNtBin == -1){
					nextNtBin = memo;
					nextNt = Leon::bin2nt(nextNtBin);
				}
				else{
					nextNt = Leon::bin2nt(nextNtBin);
				}

			}
		}
		
	}*/
	
	//return pos;
	
}



int DnaEncoder::voteMutations(int pos, int depth, bool rightExtend){
	kmer_type kmer;
	int maxScore = 0;
	int votes[4];
	
	int bestNt;
	bool isValid[4];
	
	//kmer_type mutatedKmers[4];
	//vector<int> mutations;
	//bool isMutateChoice[4];
	kmer = _kmers[pos];
	/*
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		//mutated_kmer.printASCII(_kmerSize);
		
		if(_bloom->contains(mutatedKmerMin)){
			mutations.push_back(nt);
			mutatedKmers[nt] = mutatedKmer;
			votes[nt] = 0;
			//isMutateChoice[nt] = true;
		}
	}
	*/
	//kmer_type currentKmer;
	//cout << _readseq << endl;
	//cout << pos << ": " << kmer.toString(_kmerSize) << endl;
	
	/*
	for(int nt=0; nt<4; nt++){

		//mutatedKmer.printASCII(_kmerSize);

		if(res4[nt]){

			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);

			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;


			if(Leon::bin2nt(nt) == nextNt){
				isKmerSolid = true;
			}
		}

	}


	//for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
	for(int j=0; j<depth; j++){
		//char nextNt;
		//int kmerPos;
		if(rightExtend){
			if(pos+1+_kmerSize+j >= _readSize) break;
			nextNt = _readseq[pos+1+_kmerSize+j];
		}
		else{
			if(pos-2-j < 0) break;
			nextNt = _readseq[pos-2-j];
		}

		std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);

		for(int nt=0; nt<4; nt++){

			if(res4[nt]){

			}
		}

		//cout << j << ": " << nextNt << endl;
		codeSeedNT(&_kmerModel, &mutatedKmer, nextNt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));

		if(_bloom->contains(mutatedKmerMin)){
			votes[nt] += 1;
			if(votes[nt] > maxScore){
				maxScore = votes[nt];
				bestNt = nt;
			}
		}
	}
	//}
	
	
	if(maxScore == 0){
		//cout << "No best NT" << endl;
		bestNt = -1;
	}
	//else
		//cout << "Best nt: " << bin2NT[bestNt] << endl;
	*/
	return bestNt;
	
	
}

void DnaEncoder::encodeNoAnchorRead(){
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\tEncode no anchor read" << endl;
	#endif

	//printf("encode no anchor read \n");
	//Reinsert N because they can be encoded by the coder
	for(int i=0; i<_Npos.size(); i++){
		_readseq[_Npos[i]] = 'N';
	}
	
	#ifdef LEON_PRINT_STAT
		_rangeEncoder6.encode(_readTypeModel, 1);
	#endif
	_rangeEncoder.encode(_readTypeModel, 1);
	
	//_leon->_readWithoutAnchorSize += _readSize*0.375;
	_readWithoutAnchorCount +=1;
	
	/*
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] == 'N'){
			_leon->_noAnchor_with_N_kmer_count += 1;
			break;
		}
	}
	
	bool full_N = true;
	for(int i=0; i<_readSize; i++){
		if(_readseq[i] != 'N'){
			full_N = false;
			break;
		}
	}
	if(full_N){
		_leon->_noAnchor_full_N_kmer_count += 1;
	}*/
	
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder5, _noAnchorReadSizeValueModel, _readSize);
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _noAnchorReadSizeValueModel, _readSize);
	
			
	for(int i=0; i<_readSize; i++){
		
		#ifdef LEON_PRINT_STAT
			_rangeEncoder5.encode(_noAnchorReadModel, Leon::nt2bin(_readseq[i]));
		#endif
		_rangeEncoder.encode(_noAnchorReadModel, Leon::nt2bin(_readseq[i]));
		
	}
	
}


QualDecoder::QualDecoder(Leon* leon, const string& inputFilename)
{
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	
	//cout << " qualdecoder will read from " << inputFilename  << endl;
	_leon = leon;

	_inbuffer = NULL;
}


QualDecoder::~QualDecoder(){

	free(_inbuffer);
	delete _inputFile;
}


void QualDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount){
	
	_processedSequenceCount = 0;
	
	_inputFile->seekg(blockStartPos, _inputFile->beg);
	
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;

	
	//_leon->_progress_decode->inc(1);

	//cout << "\t-----------------------" << endl;
	//cout << "\tDecoding Qual block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	
	_inbuffer = (char * ) realloc(_inbuffer, blockSize* sizeof(char));
	
	_sequenceCount = sequenceCount;
}



void QualDecoder::execute(){
	

//	printf("execute qual decoder _blockStartPos %llu  _blockSize %llu \n",_blockStartPos,_blockSize);

	
	_inputFile->read(_inbuffer,_blockSize );
	
	//printf("----Begin decomp of Block     ----\n");

	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	
	//deflateinit2 to be able to gunzip it fro mterminal
	
	//if (inflateInit2(&zs, (15+32)) != Z_OK)
	if (inflateInit (&zs) != Z_OK)
		throw    Exception ("inflate Init failed while decompressing.");
	
	zs.next_in = (Bytef*) _inbuffer ;
	zs.avail_in = _blockSize ;    // set the z_stream's input
	
	int ret;
	char outbuffer[32768];
	
	// retrieve the compressed bytes blockwise
	do {
		zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
		zs.avail_out = sizeof(outbuffer);
		
		ret = inflate (&zs, Z_SYNC_FLUSH); //ou Z_FINISH ? Z_SYNC_FLUSH
		
		if (_buffer.size() < zs.total_out) {
			// append the block to the output string
			_buffer.append(outbuffer,
							 zs.total_out - _buffer.size());
		}
		if (ret != Z_OK)
		{
			//printf("ret val %i  _blockStartPos %llu \n",ret,_blockStartPos);
			break;
		}
		else
		{
			//printf("-----block ret ok  _blockStartPos %llu  ----\n",_blockStartPos);
		}
	} while (ret == Z_OK);
	
	inflateEnd(&zs);

	_finished = true;
	
	//printf("Should be done decompressing block, in size   %llu   out size    %lu  \n",_blockSize,_buffer.size() );

}


//====================================================================================
// ** DnaDecoder
//====================================================================================

DnaDecoder::DnaDecoder(Leon* leon, const string& inputFilename) :
AbstractDnaCoder(leon)
{
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	
	_anchorDictFile = new ifstream(_leon->_anchorDictFilename.c_str(), ios::in);
	
}

DnaDecoder::DnaDecoder(Leon* leon, Requests* req, const string& inputFilename) :
AbstractDnaCoder(leon)
{
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_finished = false;
	
	_anchorDictFile = new ifstream(_leon->_anchorDictFilename.c_str(), ios::in);

	_decodeReq = true;
	_requests = req;
	
}

DnaDecoder::~DnaDecoder(){
	//delete _rangeDecoder;
	//delete _outputFile;
	delete _inputFile;
	delete _anchorDictFile;
}

void DnaDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount){

	//cout << "debug - DnaDecoder - setup - blockStartPos : " << blockStartPos << endl;
	//cout << "debug - DnaDecoder - setup - blockSize : " << blockSize << endl;
	//cout << "debug - DnaDecoder - setup - sequenceCount : " << sequenceCount << endl;


	startBlock();
	_rangeDecoder.clear();
	
	_inputFile->seekg(blockStartPos, _inputFile->beg);
	_rangeDecoder.setInputFile(_inputFile);
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t-----------------------" << endl;
		cout << "\tDecoding Dna block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	#else
		;//_leon->_progress_decode->inc(1);
	#endif
	
	_sequenceCount = sequenceCount;
		
}

//TODO : remove dupplication code with decodeAnchorRead
bool DnaDecoder::getNextReadInfos(struct ReadInfos* ri){
	//cerr << "\tdebug DnaDecoder::getNextReadInfos - BEGIN " << endl;
	if (_processedSequenceCount < _sequenceCount){
		ri->readType = _rangeDecoder.nextByte(_readTypeModel);
		//cerr << "\tdebug DnaDecoder::getNextReadInfos - ri->readType :  " << (int) ri->readType << endl;
		if(ri->readType == 0){
			/*
			ri->readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
			cerr << "\tdebug DnaDecoder::getNextReadInfos - ri->readSize :  " << (int) ri->readSize << endl;
			ri->anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
			cerr << "\tdebug DnaDecoder::getNextReadInfos - ri->anchorPos :  " << (int) ri->anchorPos << endl;
			ri->anchorAddress = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
			cerr << "\tdebug DnaDecoder::getNextReadInfos - ri->anchorAddress :  " << (int) ri->anchorAddress << endl;
			//TODO
			//suppress the requests method "getAnchor"...
			kmer_type anchor;
			if (_decodeReq){//cerr << "\tdebug DnaDecoder::getNextReadInfos - before anchor" << endl;
				ri->anchor = _requests->getAnchor(_anchorDictFile, ri->anchorAddress); 
			}
			else{
				ri->anchor = _leon->getAnchor(_anchorDictFile, ri->anchorAddress); //laa
			}
			//Decode the bit that says if the anchor is revcomp or not
			if((ri->revcomp = _rangeDecoder.nextByte(_readAnchorRevcompModel)) == 1){
				ri->revAnchor = ri->anchor;
				ri->anchor = revcomp(ri->anchor, _kmerSize);
				//anchor = ri->anchor;
			}
			else{
				ri->revAnchor = revcomp(ri->anchor, _kmerSize);
				//anchor = ri->anchor;
			}
			
			_currentSeq = ri->anchor.toString(_kmerSize);	
			ri->leftErrorPos.clear();
			ri->Npos.clear();
			//cerr << "\tdebug DnaDecoder::getNextReadInfos - anchor :  " << ri->anchor.toString(_kmerSize) << endl;

			//Decode N pos
			_prevNpos = 0;
			//cerr << "\tdebug DnaDecoder::getNextReadInfos - npos count before :  " << ri->NposCount << endl;

			ri->NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);

//cerr << "\tdebug DnaDecoder::getNextReadInfos - npos count :  " << ri->NposCount << endl;
			for(int i=0; i<ri->NposCount; i++){
				u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
				ri->Npos.push_back(nPos);
				_prevNpos = nPos;
			}
			
			//Decode error pos
			ri->nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);

			_prevErrorPos = 0;
			for(int i=0; i<ri->nbLeftError; i++){
				u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
				ri->leftErrorPos.push_back(errorPos);
				_prevErrorPos = errorPos;
			}

			//Extend anchor to the left
			kmer_type kmer = ri->anchor;
			for(int i=ri->anchorPos-1; i>=0; i--){
				kmer = extendAnchor(kmer, i, false);
			}

			//Extend anchor to the right
			kmer = ri->anchor;
			for(int i=ri->anchorPos+_kmerSize; i<ri->readSize; i++){
				kmer = extendAnchor(kmer, i, true);

			}
			cerr << "\tdebug DnaDecoder::getNextReadInfos - GOOD " << endl;
			
			//Inject N in the decoded read sequence
			for(int i=0; i<ri->NposCount; i++){
				cerr << "\tri->NposCount : " << ri->NposCount << endl;
				cerr << "\t_currentSeq.size() : " << _currentSeq.size() << endl;
				cerr << "\ti : " << i << endl;
				cerr << "\tri->Npos[i] : " << ri->Npos[i] << endl;
				_currentSeq[ri->Npos[i]] = 'N';
				cerr << "\tdebug DnaDecoder::getNextReadInfos - CRASH " << endl;
			}
			cerr << "\tdebug DnaDecoder::getNextReadInfos - CRASH " << endl;

			*/

			u_int8_t deltaType;
			u_int64_t deltaValue;
			
			//printf("Decode anchor read \n");

			//Decode read size
			//cerr << "DnaDecoder::decodeAnchorRead() - readSize before" << endl;
			_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
			ri->readSize = _readSize;
			// cerr << "DnaDecoder::decodeAnchorRead() - readSize after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - readSize : " << _readSize << endl;
			
			//Decode anchor pos
			//cerr << "DnaDecoder::decodeAnchorRead() - anchorPos before" << endl;
			int anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
			ri->anchorPos = anchorPos;
			// cerr << "DnaDecoder::decodeAnchorRead() - anchorPos after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - anchorPos : " << anchorPos << endl;
			
			//Decode anchor address
			//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress before" << endl;
			u_int64_t anchorAddress = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
			ri->anchorAddress = anchorAddress;
			//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress : " << anchorAddress << endl;

			
			//TODO
			//suppress the requests method "getAnchor"...
			kmer_type anchor;
			if (_decodeReq){
			//	cerr << "debug dnadecoder - getAnchor (requests) before" << endl;
				anchor = _requests->getAnchor(_anchorDictFile, anchorAddress); //laa
			}
			else{
			//	cerr << "debug dnadecoder - getAnchor (leon) before" << endl;
				anchor = _leon->getAnchor(_anchorDictFile, anchorAddress); //laa
			}
			//cerr << "DnaDecoder::decodeAnchorRead() - getAnchor after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - anchor : " << anchor.toString(_kmerSize) << endl;
			
			//Decode the bit that says if the anchor is revcomp or not
			if((ri->revcomp = _rangeDecoder.nextByte(_readAnchorRevcompModel)) == 1){
				anchor = revcomp(anchor, _kmerSize);
			}

			ri->anchor = anchor;
			ri->revAnchor = revcomp(anchor, _kmerSize);
			
			_currentSeq = anchor.toString(_kmerSize);	
			_leftErrorPos.clear();
			_Npos.clear();
			
			//cout << _readSize << " " << anchorPos << " " << anchorAddress << endl;
			
			//Decode N pos
			_prevNpos = 0;
			//cerr << "DnaDecoder::decodeAnchorRead() - NposCount before" << endl;
			u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
			ri->NposCount = NposCount;
			//cerr << "DnaDecoder::decodeAnchorRead() - NposCount after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - NposCount : " << (int) NposCount << endl;
			for(int i=0; i<NposCount; i++){
				//cerr << "DnaDecoder::decodeAnchorRead() - nPos before" << endl;
				u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
				//cerr << "DnaDecoder::decodeAnchorRead() - nPos after" << endl;
				//cerr << "DnaDecoder::decodeAnchorRead() - nPos : " << (int) nPos << endl;
				_Npos.push_back(nPos);
				ri->Npos.push_back(nPos);
				//cerr << nPos << endl;
				_prevNpos = nPos;
				//_Npos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
			}
			
			//Decode error pos
			//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError before" << endl;
			u_int64_t nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);
			ri->nbLeftError = nbLeftError;
			//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError after" << endl;
			//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError : " << (int) nbLeftError << endl;
			_prevErrorPos = 0;
			for(int i=0; i<nbLeftError; i++){
				u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
				addErrorPos(errorPos, true);
				ri->leftErrorPos.push_back(errorPos);
				_prevErrorPos = errorPos;
			}
		
			//Extend anchor to the left
			//cerr << "DnaDecoder::decodeAnchorRead() - Extend anchor to the left start" << endl;
			kmer_type kmer = ri->anchor;
			for(int i=anchorPos-1; i>=0; i--){
				kmer = extendAnchor(kmer, i, false);
			}
			//cerr << "DnaDecoder::decodeAnchorRead() - Extend anchor to the left end" << endl;
			//Extend anchor to the right
			kmer = anchor;
			for(int i=anchorPos+_kmerSize; i<_readSize; i++){
				kmer = extendAnchor(kmer, i, true);
				//cerr << "\t" << kmer.toString(_kmerSize) << endl;
			}

			
			//Inject N in the decoded read sequence
			//printf("npos s %i currseq %s \n",_Npos.size(),_currentSeq.c_str());
			for(int i=0; i<_Npos.size(); i++){
				_currentSeq[_Npos[i]] = 'N';
			}

			//cerr << "DnaDecoder::decodeAnchorRead() - dnadecoder END" << endl;

		}
		//if not anchored read
		else if(ri->readType == 1){
				
			_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _noAnchorReadSizeValueModel);

			for(int i=0; i<_readSize; i++){
				_currentSeq += Leon::bin2nt(_rangeDecoder.nextByte(_noAnchorReadModel));
			}

			ri->readSize = _readSize;
			ri->anchorPos = -1;
			ri->anchorAddress = -1;
			ri->anchor;
			ri->revAnchor;
			ri->revcomp = -1;
			ri->NposCount = -1;
			ri->Npos;
			ri->leftErrorPos;
			ri->nbLeftError = -1;

		}

		endRead();

		ri->sread = _buffer;

		_buffer.clear();

		return true;
	}

	else{

		_finished = true;

		return false;
	}
}

bool DnaDecoder::getNextOrderedReadsInfosBLock(struct OrderedReadsInfosBlock* orib){
	
	if (_processedSequenceCount < _sequenceCount){
	
		Hash16<kmer_type, u_int32_t >  * anchorKmersSorted = _requests->_anchorKmersSortedD;

		//cerr << "\tDnaDecoder::execute() - before decoding anchor" << endl;
		u_int64_t anchor_uint64t = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorKmerTypeModel);
		//cerr << "\tDnaDecoder::getNextOrderedReadsInfosBLock() - anchor_uint64t : " << anchor_uint64t << endl;

		_anchor.setVal(anchor_uint64t);
		orib->anchor.setVal(anchor_uint64t);
		//cerr << "\tDnaDecoder::getNextOrderedReadsInfosBLock() - anchor : " << _anchor.toString(_kmerSize) << endl;

		//get the number of reads encoded with actual anchor
		//and verify if anchor is revcomp

		u_int32_t nbReads;
		//cerr << "\tDnaDecoder::getNextOrderedReadsInfosBLock() - anchorKmersSorted->contains(_anchor) : " << 	anchorKmersSorted->contains(_anchor) << endl; 
		//cerr << "\tDnaDecoder::getNextOrderedReadsInfosBLock() - before getting nbReads" << endl; 
		anchorKmersSorted->get(_anchor, &nbReads);
		orib->nbReads = nbReads;
		//cerr << "\tDnaDecoder::getNextOrderedReadsInfosBLock() - nbReads : " << orib->nbReads << endl;
		return true;

	}
	else{

		_finished = true;

		return false;
	}
}

bool DnaDecoder::getNextOrderedReadInfos(struct OrderedReadInfos* ori){
	
	if (_processedSequenceCount < _sequenceCount){

		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq : " << _currentSeq << endl;

		_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
		ori->readSize = _readSize;
		//cerr << "\tDnaEncoder::getNextOrderedReadInfos - readSize : " << ori->readSize << endl;
		
		int anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
		ori->anchorPos = anchorPos;
		//cerr << "\tDnaEncoder::getNextOrderedReadInfos - anchorPos : " << ori->anchorPos << endl;

		int isRevComp = _rangeDecoder.nextByte(_readAnchorRevcompModel);
		ori->revcomp = isRevComp;
		//cerr << "\tDnaEncoder::getNextOrderedReadInfos - revcomp : " << ori->revcomp << endl;

		if(ori->revcomp){
			ori->anchor = revcomp(_anchor, _kmerSize);
			ori->revAnchor = _anchor;
		}
		else{
			ori->anchor = _anchor;
			ori->revAnchor = revcomp(_anchor, _kmerSize);
		}
		
		_currentSeq = ori->anchor.toString(_kmerSize);
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq = current anchor" << endl;
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq : " << _currentSeq << endl;
		_leftErrorPos.clear();
		_Npos.clear();	

		_prevNpos = 0;
		u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - NposCount : " << (int) NposCount << endl;
		for(int i=0; i<NposCount; i++){

			u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
			//cerr << "DnaDecoder::getNextOrderedReadInfos() - nPos : " << (int) nPos << endl;
			_Npos.push_back(nPos);
			ori->Npos.push_back(nPos);
			_prevNpos = nPos;
		}

		u_int64_t nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);
		//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError after" << endl;
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - nbLeftError : " << (int) nbLeftError << endl;
		_prevErrorPos = 0;
		for(int i=0; i<nbLeftError; i++){
			u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
			addErrorPos(errorPos, true);
			ori->leftErrorPos.push_back(errorPos);
			_prevErrorPos = errorPos;
		}

		//Extend anchor to the left
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - Extend anchor to the left" << endl;
		kmer_type kmer = ori->anchor;
		for(int i=anchorPos-1; i>=0; i--){
			kmer = extendAnchor(kmer, i, false);
		}
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq : " << _currentSeq << endl;
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - Extend anchor to the right" << endl;
		//Extend anchor to the right
		kmer = ori->anchor;
		for(int i=anchorPos+_kmerSize; i<_readSize; i++){
			kmer = extendAnchor(kmer, i, true);
			//cerr << "\t" << kmer.toString(_kmerSize) << endl;
		}
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq : " << _currentSeq << endl;

			
		//Inject N in the decoded read sequence
		//cerr << "DnaDecoder::getNextOrderedReadInfos() - Inject Ns" << endl;
		//printf("npos s %i currseq %s \n",_Npos.size(),_currentSeq.c_str());
		for(int i=0; i<_Npos.size(); i++){
			_currentSeq[_Npos[i]] = 'N';
		}

		//cerr << "DnaDecoder::getNextOrderedReadInfos() - _currentSeq : " << _currentSeq << endl;

		endRead();

		ori->sread = _buffer;

		_buffer.clear();
		//cerr << "\tDnaEncoder::encodeSortedFileRead - _processedSequenceCount : " << _processedSequenceCount << endl;

		return true;

	}
	else{

		_finished = true;

		return false;
	}
}

int _nbTests;

void DnaDecoder::execute(){

	//decodeFirstHeader();
	cerr << "\tDnaDecoder::execute() - start" << endl;

	//cerr << "\tDnaDecoder::execute() - _sequenceCount : " << _sequenceCount << endl;
	if (_orderReads){
		//cerr << "\tDnaDecoder::execute() - reading sorted file" << endl;

		//TODO : temporarily using _requests dict, until dict is encoded

		//Hash16<kmer_type, u_int32_t >  * anchorKmersSorted = _requests->_anchorKmersSorted;
		Hash16<kmer_type, u_int32_t >  * anchorKmersSorted = _requests->_anchorKmersSortedD;

		bool decodingNewAnchor = true;
		u_int32_t nbReads = 0; 
		//int nbTestsMax = 2000;
		_nbTests = 0;
		while(_processedSequenceCount < _sequenceCount /*&& _nbTests < nbTestsMax*/){
			if (decodingNewAnchor){
				//cerr << "\n\n\tDnaDecoder::execute() - NB TEST : " << _nbTests << endl;
				//cerr << "\n\n\tDnaDecoder::execute() - decode anchor" << endl;

				u_int64_t anchor_uint64t = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorKmerTypeModel);
				//cerr << "\tDnaDecoder::execute() - anchor_uint64t : " << anchor_uint64t << endl;

				_anchor.setVal(anchor_uint64t);
				//cerr << "\tDnaDecoder::execute() - anchor : " << _anchor.toString(_kmerSize) << endl;
				decodingNewAnchor = false;
				++_nbTests;
			}
			else{
				//get the number of reads encoded with actual anchor
				//and verify if anchor is revcomp

				u_int32_t nbReads; 
				anchorKmersSorted->get(_anchor, &nbReads);
				//cerr << "\n\tDnaDecoder::execute() - nbReads to decode : " << nbReads << endl;

				for (int i=0; i<nbReads; ++i){
					//if (_nbTests < nbTestsMax){
						//cerr << "\n\n\tDnaDecoder::execute() - NB TEST : " << _nbTests << endl;
						//cerr << "\n\tDnaDecoder::execute() - decode read nb " << i+1 << endl;
						//cerr << "\tDnaDecoder::execute() - decodeSortedAnchorRead()" << endl;
						decodeSortedAnchorRead();
						endRead();
						++_nbTests;	
					//}
				}
				decodingNewAnchor = true;		
			}
			//cerr << "\tDnaDecoder::execute() - iteration nb : " << _nbTests << endl;
			//cerr << "\tDnaDecoder::execute() - _processedSequenceCount : " << _processedSequenceCount << endl;
			//cerr << "\tDnaDecoder::execute() - _sequenceCount : " << _sequenceCount << endl;
			//cerr << "\tDnaDecoder::execute() - _buffer : " << _buffer.c_str() << endl; 
		}
	}
	else{
		while(_processedSequenceCount < _sequenceCount){
		
			//cout << "lala" << endl;
			//int i=0;
			//while(i < Leon::READ_PER_BLOCK){
			//while(_inputFile->tellg() <= _blockStartPos+_blockSize){
			//if(_leon->_readCount > 1) return;



			cerr << "\tDnaDecoder::execute() - reading unsorted file" << endl;
			cerr << "\tDnaDecoder::execute() - readType before" << endl;
			u_int8_t readType = _rangeDecoder.nextByte(_readTypeModel);
			cerr << "\tDnaDecoder::execute() - readType after" << endl;
			cerr << "\tDnaDecoder::execute() - readType : " << (int) readType << endl;
			cerr << "Read type: " << (int)readType << endl;

			cerr << "\tDnaDecoder::execute() - before decodeAnchor" << endl;
			if(readType == 0){
				cerr << "\tDnaDecoder::execute() - decodeAnchorRead()" << endl;
				decodeAnchorRead(); //ici		
			}
			else if(readType == 1){
				cerr << "\tDnaDecoder::execute() - decodeNoAnchorRead();" << endl;
				decodeNoAnchorRead();
			}
			

			cerr << "\tDnaDecoder::execute() - after decodeAnchor" << endl;	
			endRead();
			cerr << _inputFile->tellg() << " " << _blockStartPos+_blockSize << endl;
			
		}	
		/*
		string trueSeq = string((*_leon->_testBankIt)->getDataBuffer());
		trueSeq = trueSeq.substr(0, _readSize);
		//cout << trueSeq << endl;
		//cout << _currentSeq << endl;
		if(trueSeq != _currentSeq){
			cout << (*_leon->_testBankIt)->getIndex() << "\t\tseq different !!" << endl;
			cout << "\t\t" << trueSeq << endl;
			cout << "\t\t" << _currentSeq << endl;
			_leon->_readCount += 1;
			return;
		}
		_leon->_testBankIt->next();
		*/
		#ifdef PRINT_DEBUG_DECODER
			_readCount += 1;
			cout << _leon->_readCount << ": " << _currentSeq << endl;
		#endif
			
		//i++;
		//_leon->_readCount += 1;
		//if(i == 1) return;
		//_currentSeq.clear();
			
		//cout << (int)(_inputFile->tellg() < _blockStartPos+_blockSize) << endl;
		
		
	}
	
	cerr << "\tDnaDecoder::execute() - end" << endl;
	_finished = true;
	
}

void DnaDecoder::decodeAnchorRead(){


	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecode anchor read" << endl;
	#endif

	u_int8_t deltaType;
	u_int64_t deltaValue;
	
	//printf("Decode anchor read \n");

	//Decode read size
	//deltaType = _rangeDecoder.nextByte(_readSizeDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
//	printf("read deltaValue %llu \n",deltaValue);
	//_readSize = CompressionUtils::getValueFromDelta(deltaType, _prevReadSize, deltaValue);
	//_prevReadSize = _readSize;
	//cerr << "DnaDecoder::decodeAnchorRead() - readSize before" << endl;
	_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
	// cerr << "DnaDecoder::decodeAnchorRead() - readSize after" << endl;
	//cerr << "DnaDecoder::decodeAnchorRead() - readSize : " << _readSize << endl;
//	printf("read size %i \n",_readSize);
	
	//Decode anchor pos
	//deltaType = _rangeDecoder.nextByte(_anchorPosDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
	//int anchorPos = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorPos, deltaValue);
	//_prevAnchorPos = anchorPos;
//	printf("anchor pos %i \n",anchorPos);
	//cerr << "DnaDecoder::decodeAnchorRead() - anchorPos before" << endl;
	int anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
	// cerr << "DnaDecoder::decodeAnchorRead() - anchorPos after" << endl;
	// cerr << "DnaDecoder::decodeAnchorRead() - anchorPos : " << anchorPos << endl;
	
	//Decode anchor address
	//deltaType = _rangeDecoder.nextByte(_anchorAddressDeltaTypeModel);
	//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
	//u_int64_t anchorAddress = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorAddress, deltaValue);
	//_prevAnchorAddress = anchorAddress;
	//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress before" << endl;
	u_int64_t anchorAddress = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
	//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress after" << endl;
	//cerr << "DnaDecoder::decodeAnchorRead() - anchorAddress : " << anchorAddress << endl;

	
	//TODO
	//suppress the requests method "getAnchor"...
	kmer_type anchor;
	if (_decodeReq){
	//	cout << "debug dnadecoder - getAnchor (requests) before" << endl;
		anchor = _requests->getAnchor(_anchorDictFile, anchorAddress); //laa
	}
	else{
	//	cout << "debug dnadecoder - getAnchor (leon) before" << endl;
		anchor = _leon->getAnchor(_anchorDictFile, anchorAddress); //laa
	}
	//cerr << "DnaDecoder::decodeAnchorRead() - getAnchor after" << endl;
	//cerr << "DnaDecoder::decodeAnchorRead() - anchor : " << anchor.toString(_kmerSize) << endl;
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tRead size: " << _readSize << endl;
		cout << "\t\t\tAnchor pos: " << anchorPos << endl;
		cout << "\t\t\tAnchor adress: " << anchorAddress << endl;
		cout << "\t\t\tAnchor: " << anchor.toString(_kmerSize) << endl;
	#endif
	
	//Decode the bit that says if the anchor is revcomp or not

	if(_rangeDecoder.nextByte(_readAnchorRevcompModel) == 1)
		anchor = revcomp(anchor, _kmerSize);
	
	_currentSeq = anchor.toString(_kmerSize);	
	_leftErrorPos.clear();
	//_rightErrorPos.clear();
	_Npos.clear();
	
	cout << _readSize << " " << anchorPos << " " << anchorAddress << endl;
	//Decode N pos
	_prevNpos = 0;
	//cout << "DnaDecoder::decodeAnchorRead() - NposCount before" << endl;
	u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	//cerr << "DnaDecoder::decodeAnchorRead() - NposCount after" << endl;
	//cerr << "DnaDecoder::decodeAnchorRead() - NposCount : " << (int) NposCount << endl;
	for(int i=0; i<NposCount; i++){
		//deltaType = _rangeDecoder.nextByte(_NposDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel);
		//u_int64_t nPos = CompressionUtils::getValueFromDelta(deltaType, _prevNpos, deltaValue);
		//cerr << "DnaDecoder::decodeAnchorRead() - nPos before" << endl;
		u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
		//cerr << "DnaDecoder::decodeAnchorRead() - nPos after" << endl;
		//cerr << "DnaDecoder::decodeAnchorRead() - nPos : " << (int) nPos << endl;
		_Npos.push_back(nPos);
		//cout << nPos << endl;
		_prevNpos = nPos;
		//_Npos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	
	//Decode error pos
	//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError before" << endl;
	u_int64_t nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);
	//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError after" << endl;
	//cerr << "DnaDecoder::decodeAnchorRead() - nbLeftError : " << (int) nbLeftError << endl;
	_prevErrorPos = 0;
	for(int i=0; i<nbLeftError; i++){
		u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
		addErrorPos(errorPos, true);
		_prevErrorPos = errorPos;
	}

	/*
	u_int64_t nbRightError = CompressionUtils::decodeNumeric(_rangeDecoder, _rightErrorModel);
	_prevErrorPos = anchorPos-1;
	for(int i=0; i<nbLeftError; i++){
		int deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel);
		int errorPos = _prevErrorPos-deltaValue;
		//deltaType = _rangeDecoder.nextByte(_errorPosDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _errorPosModel); //reprise
		//u_int64_t errorPos = CompressionUtils::getValueFromDelta(deltaType, _prevErrorPos, deltaValue);
		addErrorPos(errorPos, false);
		//_errorPos.push_back(errorPos);
		_prevErrorPos = errorPos;
		//_errorPos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	_prevErrorPos = anchorPos+_kmerSize;
	for(int i=0; i<nbRightError; i++){
		int deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _rightErrorPosModel);
		int errorPos = _prevErrorPos+deltaValue;
		//deltaType = _rangeDecoder.nextByte(_errorPosDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _errorPosModel); //reprise
		//u_int64_t errorPos = CompressionUtils::getValueFromDelta(deltaType, _prevErrorPos, deltaValue);
		addErrorPos(errorPos, true);
		//_errorPos.push_back(errorPos);
		_prevErrorPos = errorPos;
		//_errorPos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}*/

	
	//Extend anchor to the left
	//cerr << "DnaDecoder::decodeAnchorRead() - Extend anchor to the left start" << endl;
	kmer_type kmer = anchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = extendAnchor(kmer, i, false);
	}
	//cout << "DnaDecoder::decodeAnchorRead() - Extend anchor to the left end" << endl;
	//Extend anchor to the right
	kmer = anchor;
	for(int i=anchorPos+_kmerSize; i<_readSize; i++){
		kmer = extendAnchor(kmer, i, true);
		//cout << "\t" << kmer.toString(_kmerSize) << endl;
	}
	
	//Inject N in the decoded read sequence
	//printf("npos s %i currseq %s \n",_Npos.size(),_currentSeq.c_str());
	for(int i=0; i<_Npos.size(); i++){
		_currentSeq[_Npos[i]] = 'N';
	}

	//cerr << "DnaDecoder::decodeAnchorRead() - dnadecoder END" << endl;
}

void DnaDecoder::decodeSortedAnchorRead(){


	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecode anchor read" << endl;
	#endif

	//printf("Decode anchor read \n");

	//Decode read type
	/*cerr << "debug dnadecoder - readType before" << endl;
	u_int8_t readType = _rangeDecoder.nextByte(_readTypeModel);
	cerr << "debug dnadecoder - readType after" << endl;
	cerr << "debug dnadecoder - readType : " << (int) readType << endl;*/

	//Decode read size
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - readSize before" << endl;
	_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - readSize after" << endl;
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - readSize : " << _readSize << endl;
	
	//Decode anchor pos
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - anchorPos before" << endl;
	int anchorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - anchorPos after" << endl;
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - anchorPos : " << anchorPos << endl;
	
	//Decode the bit that says if the anchor is revcomp or not

	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - revcomp before" << endl;
	int isRevcomp = _rangeDecoder.nextByte(_readAnchorRevcompModel);
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - isRevComp : " << isRevcomp << endl;
	/*if (_nbTests == 6){
		isRevcomp = 1;
	}*/

	kmer_type currentAnchor;
	if(isRevcomp){
		currentAnchor = revcomp(_anchor, _kmerSize);
	}
	else{
		currentAnchor = _anchor;
	}
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - currentAnchor : " << currentAnchor << endl;
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - revcomp after" << endl;
	
	_currentSeq = currentAnchor.toString(_kmerSize);	
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - _currentSeq : " << _currentSeq << endl;
	_leftErrorPos.clear();
	//_rightErrorPos.clear();
	_Npos.clear();

	//cout << _readSize << " " << anchorPos << " " << anchorAddress << endl;
	//Decode N pos
	_prevNpos = 0;
	//cout << "debug dnadecoder - NposCount before" << endl;
	u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - NposCount : " << NposCount << endl;
	//cout << "debug dnadecoder - NposCount after" << endl;
	//cout << "NposCount : "  << (int) NposCount << endl;
	//cout << "debug dnadecoder - NposCount : " << (int) NposCount << endl;
	for(int i=0; i<NposCount; i++){
		//deltaType = _rangeDecoder.nextByte(_NposDeltaTypeModel);
		//deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel);
		//u_int64_t nPos = CompressionUtils::getValueFromDelta(deltaType, _prevNpos, deltaValue);
	//	cout << "debug dnadecoder - nPos before" << endl;
		u_int64_t nPos = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel) + _prevNpos;
	//		cerr << "\tDnaDecoder::decodeSortedAnchorRead() - nPos : " << nPos << endl;
	//	cout << "debug dnadecoder - nPos after" << endl;
	//	cout << "debug dnadecoder - nPos : " << (int) nPos << endl;
		_Npos.push_back(nPos);
		//cout << nPos << endl;
		_prevNpos = nPos;
		//_Npos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	
	//Decode error pos
	//cout << "debug dnadecoder - nbLeftError before" << endl;
	u_int64_t nbLeftError = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorModel);
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - nbLeftError : " << nbLeftError << endl;
	//cout << "debug dnadecoder - nbLeftError after" << endl;
	//cout << "debug dnadecoder - nbLeftError : " << (int) nbLeftError << endl;
	_prevErrorPos = 0;
	for(int i=0; i<nbLeftError; i++){
		u_int64_t errorPos = CompressionUtils::decodeNumeric(_rangeDecoder, _leftErrorPosModel) + _prevErrorPos;
	//	cerr << "\tDnaDecoder::decodeSortedAnchorRead() - errorPos : " << errorPos << endl;
		addErrorPos(errorPos, true);
		_prevErrorPos = errorPos;
	}

	//Extend anchor to the left
	//cout << "debug Extend anchor to the left start" << endl;
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - extending anchor : " << endl;
	kmer_type kmer = currentAnchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = extendAnchor(kmer, i, false);
	}
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - after left extend : " << kmer.toString(_kmerSize) << endl;
	//cout << "debug Extend anchor to the left end" << endl;
	//Extend anchor to the right
	kmer = currentAnchor;
	for(int i=anchorPos+_kmerSize; i<_readSize; i++){
		kmer = extendAnchor(kmer, i, true);
		//cout << "\t" << kmer.toString(_kmerSize) << endl;
	}
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - after rigth extend : " << kmer.toString(_kmerSize) << endl;
	//Inject N in the decoded read sequence
	//printf("npos s %i currseq %s \n",_Npos.size(),_currentSeq.c_str());
	for(int i=0; i<_Npos.size(); i++){
		_currentSeq[_Npos[i]] = 'N';
	}
	//cerr << "\tDnaDecoder::decodeSortedAnchorRead() - after writing Ns : " << _currentSeq << endl;

	//cout << "debug dnadecoder END" << endl;
}

kmer_type DnaDecoder::extendAnchor(kmer_type kmer, int pos, bool rightExtend){
	
	//cerr << endl << "DnaDecoder::extendAnchor - extendAnchor BEGIN" << endl;

	u_int8_t nextNt;
	//int nextNtBin;
	kmer_type resultKmer;
		
	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
		//cerr << "\tDnaDecoder::extendAnchor - next nt is N " << endl;
		nextNt = 'A';
		if(rightExtend){
			_currentSeq += nextNt;
		}
		else{
			_currentSeq.insert(_currentSeq.begin(), nextNt);
		}
		//cout << _currentSeq << endl;
		//if(nextNt == 'N') nextN
		
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
		//cout << kmer.toString(_kmerSize) << endl;
		//cout << resultKmer.toString(_kmerSize) << endl;
		return resultKmer;
	}
	


	/*
	if(rightExtend){
		if(std::find(_rightErrorPos.begin(), _rightErrorPos.end(), pos) != _rightErrorPos.end()){
			nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));

			if(rightExtend)
				_currentSeq += nextNt;
			else
				_currentSeq.insert(_currentSeq.begin(), nextNt);

			std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);

			for(int nt=0; nt<4; nt++){

				if(res4[nt]){
					kmer_type mutatedKmer = kmer;
					codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
					return mutatedKmer;
				}
			}
		}
	}
	else{*/
		if(std::find(_leftErrorPos.begin(), _leftErrorPos.end(), pos) != _leftErrorPos.end()){

			nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));
			//cerr << "\tDnaDecoder::extendAnchor - next nt is an error : " << nextNt << endl;

			if(rightExtend)
				_currentSeq += nextNt;
			else
				_currentSeq.insert(_currentSeq.begin(), nextNt);

			std::bitset<4> res4;

			for(int nt=0; nt<4; nt++){

				kmer_type mutatedKmer = kmer;
				codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
				Node node;
				//kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
				kmer_type mutatedKmerMin = mutatedKmer;
				node = Node(Node::Value(mutatedKmerMin));
		
		//mutatedKmer.printASCII(_kmerSize);
		
				res4[nt] = _leon->_graph.contains(node);

				if(res4[nt]){
					kmer_type mutatedKmer = kmer;
					codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
					return mutatedKmer;
				}
			}
		}
	//}

	
		
	//cout << kmer.toString(_kmerSize) << endl;
	kmer_type uniqKmer, mutatedSolidKmer;
	int uniqNt;
	//bool isKmerSolid = false;
	

	//kmer = _kmers[pos];

	int indexedKmerCount = 0;
	
	//cout << kmer.toString(_kmerSize) << endl;
	
	
	
	std::bitset<4> res4;
	for(int nt=0; nt<4; nt++){

		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		Node node;
		//kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		kmer_type mutatedKmerMin = mutatedKmer;
		node = Node(Node::Value(mutatedKmerMin));
		
		//mutatedKmer.printASCII(_kmerSize);

		//cerr << "DnaDecoder::extendAnchor - extendAnchor before graph copntains" << endl;
		res4[nt] = _leon->_graph.contains(node);
		//cerr << "DnaDecoder::extendAnchor - extendAnchor after graph copntains" << endl;

		if(res4[nt]){
			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
			
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
		}
	}
		//cout << "debug extendAnchor test1" << endl;

	
	/*
	for(int nt=0; nt<4; nt++){
		//if(nt == original_nt){
		//	continue;
		//}
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		//mutatedKmer.printASCII(_kmerSize);
		
		if(_bloom->contains(mutatedKmerMin)){
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
		}
		
	}
 */
	
	if(indexedKmerCount == 1){
		nextNt = Leon::bin2nt(uniqNt);
		//cerr << "\tDnaDecoder::extendAnchor - next nt is unique : " << nextNt << endl;
		//cout << "case 1         " << nextNt << endl;
		resultKmer = uniqKmer;
		//cerr << "\tDnaDecoder::extendAnchor - next1 nt is : " << nextNt << endl;
	}
	else if(indexedKmerCount == 2){

		char nt1 = -1;
		char nt2 = -1;

		for(int nt=0; nt<4; nt++){
			if(res4[nt]){
				if(nt1 == -1)
					nt1 = nt;
				else if(nt2 == -1)
					nt2 = nt;
				else break;
			}
		}

		u_int8_t nextBinaryNt = _rangeDecoder.nextByte(_bifurcationBinaryModel);

		//cout << (int)nextBinaryNt << endl;
		if(nextBinaryNt == 0)
			nextNt = Leon::bin2nt(nt1);
		else
			nextNt = Leon::bin2nt(nt2);
		//cerr << "\tDnaDecoder::extendAnchor - next nt binary bif : " << nextNt << endl;
		//cout << nextNt << endl;
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
		//cerr << "\tDnaDecoder::extendAnchor - next2 nt is : " << nextNt << endl;


	}
	else{
		nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_bifurcationModel));
		//cout << "case 2          "<< nextNt << endl;
		//cerr << "\tDnaDecoder::extendAnchor - next nt is normal bif : " << nextNt << endl;
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
		//cerr << "\tDnaDecoder::extendAnchor - nextn nt is : " << nextNt << endl;
	}
	
	//cout << nextNt << endl;
	
	//if(nextNt == 'N') cout << "lala" << endl;
	
	if(rightExtend){
		_currentSeq += nextNt;
	}
	else{
		_currentSeq.insert(_currentSeq.begin(), nextNt);
	}
	
	//cout << _currentSeq << endl;
	//cout << resultKmer.toString(_kmerSize) << endl;
	return resultKmer;
		
}

void DnaDecoder::decodeNoAnchorRead(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tDecode no anchor read" << endl;
	#endif
	
	_readSize = CompressionUtils::decodeNumeric(_rangeDecoder, _noAnchorReadSizeValueModel);
	//cout << "\tRead size: " << _readSize << endl;
	for(int i=0; i<_readSize; i++){
		_currentSeq += Leon::bin2nt(_rangeDecoder.nextByte(_noAnchorReadModel));
	}
	//endSeq();
	//cout << read << endl;
}
	
void DnaDecoder::endRead(){
	AbstractDnaCoder::endRead();
	
	//cerr << "DnaDecoder::endRead() - _buffer max size : " << _buffer.max_size() << endl;
	//cerr << "DnaDecoder::endRead() - _buffer size before add : " << _buffer.length() << endl;

	_buffer += _currentSeq + '\n';

	//cerr << "DnaDecoder::endRead() - _buffer size after add : " << _buffer.length() << endl;
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tRead: " << _currentSeq << endl;
	#endif
	//_outputFile->write(_currentSeq.c_str(), _currentSeq.size());
	_currentSeq.clear();
}




