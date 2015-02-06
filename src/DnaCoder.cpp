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
_noAnchorReadModel(5), _mutationModel(5), //5value: A, C, G, T, N
_readAnchorRevcompModel(2),
_readSizeDeltaTypeModel(3),
_anchorPosDeltaTypeModel(3),
_anchorAddressDeltaTypeModel(3),
_NposDeltaTypeModel(3),
_errorPosDeltaTypeModel(3)
{
	_leon = leon;
	_bloom = _leon->_bloom;
	_kmerSize = _leon->_kmerSize;
	
	
	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_anchorAddressModel.push_back(Order0Model(256));
		_anchorPosModel.push_back(Order0Model(256));
		_noAnchorReadSizeValueModel.push_back(Order0Model(256));
		_readSizeValueModel.push_back(Order0Model(256));
		_NposModel.push_back(Order0Model(256));
		_errorPosModel.push_back(Order0Model(256));
		_numericModel.push_back(Order0Model(256));
	}
	
	
}

void AbstractDnaCoder::startBlock(){
	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_anchorAddressModel[i].clear();
		_anchorPosModel[i].clear();
		_noAnchorReadSizeValueModel[i].clear();
		_readSizeValueModel[i].clear();
		_NposModel[i].clear();
		_errorPosModel[i].clear();
		_numericModel[i].clear();
	}
	_readTypeModel.clear();
	_noAnchorReadModel.clear();
	_mutationModel.clear();
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


//====================================================================================
// ** DnaEncoder
//====================================================================================
DnaEncoder::DnaEncoder(Leon* leon) :
AbstractDnaCoder(leon), _itKmer(_kmerModel), _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0), _MCmultipleNoSolid(0)
{
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

}

DnaEncoder::DnaEncoder(const DnaEncoder& copy) :
AbstractDnaCoder(copy._leon), _itKmer(_kmerModel),
 _totalDnaSize(0), _readCount(0), _MCtotal(0), _readWithoutAnchorCount(0),
_MCuniqSolid (0), _MCuniqNoSolid(0), _MCnoAternative(0), _MCmultipleSolid(0), _MCmultipleNoSolid(0)
{
	_thread_id = __sync_fetch_and_add (&_leon->_nb_thread_living, 1);

	startBlock();
	//_leon = copy._leon;
	//_bloom = copy._bloom;
	
	#ifdef LEON_PRINT_STAT
		_rangeEncoder1.updateModel = false;
		_rangeEncoder2.updateModel = false;
		_rangeEncoder3.updateModel = false;
		_rangeEncoder4.updateModel = false;
		_rangeEncoder5.updateModel = false;
	#endif
}

DnaEncoder::~DnaEncoder(){

	if(_thread_id!=0 && (_seqId+1) % Leon::READ_PER_BLOCK != 0){
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
	__sync_fetch_and_add(&_leon->_MCmultipleNoSolid, _MCmultipleNoSolid);
	
	#ifdef LEON_PRINT_STAT
		__sync_fetch_and_add(&_leon->_anchorAdressSize, _rangeEncoder3.getBufferSize());
		__sync_fetch_and_add(&_leon->_anchorPosSize, _rangeEncoder2.getBufferSize());
		__sync_fetch_and_add(&_leon->_readSizeSize, _rangeEncoder1.getBufferSize());
		__sync_fetch_and_add(&_leon->_bifurcationSize, _rangeEncoder4.getBufferSize());
		__sync_fetch_and_add(&_leon->_noAnchorSize, _rangeEncoder5.getBufferSize());
		
		_rangeEncoder1.clear();
		_rangeEncoder2.clear();
		_rangeEncoder3.clear();
		_rangeEncoder4.clear();
		_rangeEncoder5.clear();
	#endif
		
#ifndef SERIAL
	//_leon->_blockwriter->incDone(1);
	_leon->_blockwriter->waitForWriter();
	
	if(nb_remaining==1)
	{
		_leon->_blockwriter->FlushWriter();
	}
#endif

	
	
}

void DnaEncoder::operator()(Sequence& sequence){
	_sequence = &sequence;
	//cout << _sequence->getIndex() << endl;
	_seqId = _sequence->getIndex() ;
	_readSize = _sequence->getDataSize();
	_readseq = _sequence->getDataBuffer();
		
	_totalDnaSize += _readSize ;
	//_lastSequenceIndex = sequence->getIndex();
	
//	if(_sequence->getIndex() % Leon::READ_PER_BLOCK == 0){

	execute();
	
	
	if(_processedSequenceCount >= Leon::READ_PER_BLOCK ){
		
		writeBlock();
		startBlock();
	}
	
}

void DnaEncoder::writeBlock(){
	if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}
	
	int blockId = (  _seqId / Leon::READ_PER_BLOCK)   ;
	//printf("\nTid %i  WB :  blockid %i sid %llu     size: %llu \n",_thread_id, blockId, _seqId, _rangeEncoder.getBufferSize() );

	//_leon->_realDnaCompressedSize += _rangeEncoder.getBufferSize();
	_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), _processedSequenceCount,blockId);
	_rangeEncoder.clear();
	
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
	
	if(_readSize < _kmerSize){
		encodeNoAnchorRead();
		endRead();
		return;
	}
	 
	//cout << _leon->_readCount << endl;
	//kmer_type anchorKmer = 0;
	u_int32_t anchorAddress;
	
	buildKmers();
	int anchorPos = findExistingAnchor(&anchorAddress); //unsynch
	
	if(anchorPos == -1)
		anchorPos = _leon->findAndInsertAnchor(_kmers, &anchorAddress);  //unsynch
	
	//cout << anchorPos << endl;
	
	if(anchorPos == -1)
		encodeNoAnchorRead();
	else{
		encodeAnchorRead(anchorPos, anchorAddress);
	}
	
	endRead();

}

void DnaEncoder::buildKmers(){

	
	_Npos.clear();
	
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
	
	for(int i=0; i<_kmers.size(); i++){
		kmer = _kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		if(_leon->anchorExist(kmerMin, anchorAddress)){
			return i;
		}
	}
	return -1;
}


void DnaEncoder::encodeAnchorRead(int anchorPos, u_int32_t anchorAddress){
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\tEncode anchor read" << endl;
	#endif
	//printf("encode  anchor read \n");

	//encode read type (0: read with anchor, 1: read without anchor)
	_rangeEncoder.encode(_readTypeModel, 0);
	
	u_int64_t deltaValue;
	u_int8_t deltaType;
	
	//Encode read size
	deltaType = CompressionUtils::getDeltaValue(_readSize, _prevReadSize, &deltaValue);
	#ifdef LEON_PRINT_STAT
		_rangeEncoder1.encode(_readSizeDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder1, _readSizeValueModel, deltaValue);
	#endif
	_rangeEncoder.encode(_readSizeDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _readSizeValueModel, deltaValue);
	_prevReadSize = _readSize;
	//printf("read size %i  deltaValue %i\n",_readSize,deltaValue);

	//Encode anchor pos
	deltaType = CompressionUtils::getDeltaValue(anchorPos, _prevAnchorPos, &deltaValue);
	#ifdef LEON_PRINT_STAT
		_rangeEncoder2.encode(_anchorPosDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder2, _anchorPosModel, deltaValue);
	#endif
	_rangeEncoder.encode(_anchorPosDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorPosModel, deltaValue);
	_prevAnchorPos = anchorPos;
	//printf("anchor pos %i \n",anchorPos);

	//Encode anchor address
	deltaType = CompressionUtils::getDeltaValue(anchorAddress, _prevAnchorAddress, &deltaValue);
	#ifdef LEON_PRINT_STAT
		_rangeEncoder3.encode(_anchorAddressDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder3, _anchorAddressModel, deltaValue);
	#endif
	_rangeEncoder.encode(_anchorAddressDeltaTypeModel, deltaType);
	CompressionUtils::encodeNumeric(_rangeEncoder, _anchorAddressModel, deltaValue);
	_prevAnchorAddress = anchorAddress;
	//printf("anchor adress %i \n",anchorAddress);


	
	kmer_type anchor = _kmers[anchorPos];
	
	//Encode a bit that says if the anchor is normal or revcomp
	if(anchor == min(anchor, revcomp(anchor, _kmerSize)))
		_rangeEncoder.encode(_readAnchorRevcompModel, 0);
	else
		_rangeEncoder.encode(_readAnchorRevcompModel, 1);
			
	#ifdef PRINT_DEBUG_ENCODER
		cout << "\t\t\tAnchor pos: " << anchorPos << endl;
		cout << "\t\t\tAnchor: " << _kmers[anchorPos].toString(_kmerSize) << endl;
	#endif
	
	_bifurcations.clear();
	_errorPos.clear();

	
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
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _Npos.size());
	for(int i=0; i<_Npos.size(); i++){
		deltaType = CompressionUtils::getDeltaValue(_Npos[i], _prevNpos, &deltaValue);
		_rangeEncoder.encode(_NposDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder, _NposModel, deltaValue);
		_prevNpos = _Npos[i];
	}
	
	//Encode the positions of sequencing errors
	#ifdef LEON_PRINT_STAT
		CompressionUtils::encodeNumeric(_rangeEncoder4, _numericSizeModel, _numericModel, _errorPos.size());
	#endif
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _errorPos.size());
	for(int i=0; i<_errorPos.size(); i++){
		deltaType = CompressionUtils::getDeltaValue(_errorPos[i], _prevErrorPos, &deltaValue);
		#ifdef LEON_PRINT_STAT
			_rangeEncoder4.encode(_errorPosDeltaTypeModel, deltaType);
			CompressionUtils::encodeNumeric(_rangeEncoder4, _errorPosSizeModel, _errorPosModel, deltaValue);
		#endif
		_rangeEncoder.encode(_errorPosDeltaTypeModel, deltaType);
		CompressionUtils::encodeNumeric(_rangeEncoder, _errorPosModel, deltaValue);
		_prevErrorPos = _errorPos[i];
	}
	
	for(int i=0; i<_bifurcations.size(); i++){
		#ifdef LEON_PRINT_STAT
			_rangeEncoder4.encode(_mutationModel, Leon::nt2bin(_bifurcations[i]));
		#endif
		//cout << Leon::nt2bin(_bifurcations[i]) << " ";
		_rangeEncoder.encode(_mutationModel, Leon::nt2bin(_bifurcations[i]));
	}
	//cout << endl;
	
	
}
	
kmer_type DnaEncoder::buildBifurcationList(int pos, kmer_type kmer, bool rightExtend){
		
	char nextNt = _readseq[pos];
		
	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;
		//return pos;
	}
	
	kmer_type kmerMin, uniqKmer;
	int uniqNt;
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
	
	
	
	
	std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);
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
	
	
	/*
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		//kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		kmer_type mutatedKmerMin = mutatedKmer;
		
		//mutatedKmer.printASCII(_kmerSize);
		
		if(_bloom->contains(mutatedKmerMin)){
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
			
			
			if(Leon::bin2nt(nt) == nextNt){
				isKmerSolid = true;
			}
		}
		
	}*/
	
	_MCtotal +=1;
	
	if(indexedKmerCount == 1){
		if(isKmerSolid){
			_MCuniqSolid += 1;
		}
		else{
			_MCuniqNoSolid += 1;
			//_leon->_readWithAnchorMutationChoicesSize += 0.25;
			_bifurcations.push_back(nextNt);
			_errorPos.push_back(pos);
		}
		//codeSeedNT(&_kmerModel, &kmer, uniqNt, rightExtend);
		return uniqKmer;
	}
	else{
		if(indexedKmerCount == 0){
			_MCnoAternative += 1;
		}
		else{
			if(isKmerSolid){
				_MCmultipleSolid += 1;
				
				/*
				if(_solidMutaChainLockTime <= 0){
					
					if(extendMutaChain(kmer, pos, rightExtend)){
						#ifdef PRINT_DEBUG_EXTMUTA
							cout << "\t\t" << _solidMutaChainStartPos << " " << _solidMutaChainSize << endl;
						#endif
						if(rightExtend)
							return pos + _solidMutaChainSize;
						else
							return pos - _solidMutaChainSize;
					}
					else{
						_solidMutaChainLockTime = 1;
						//cout << _solidMutaChainLockTime << endl;
						//return pos;*/
						/*
						//cout << _solidMutaChainSize << endl;
						if(rightExtend){
							for(int i=_solidMutaChainStartPos; i<_solidMutaChainStartPos+_solidMutaChainSize-1; i++){
								#ifdef PRINT_DEBUG_EXTMUTA
									cout << _readseq[i] << endl;
								#endif
								_bifurcations.push_back(_readseq[i]);
							}
							return _solidMutaChainStartPos+_solidMutaChainSize-2;
						}
						else{
							for(int i=_solidMutaChainStartPos; i>_solidMutaChainStartPos-_solidMutaChainSize+1; i--){
								#ifdef PRINT_DEBUG_EXTMUTA
									cout << _readseq[i] << endl;
								#endif
								_bifurcations.push_back(_readseq[i]);
							}
							return _solidMutaChainStartPos-_solidMutaChainSize+2;
						}*//*
					}
				}
				else{
					if(_solidMutaChainLockTime > 0) _solidMutaChainLockTime -= 1;
				}
				*/
					
			}
			else{
				_MCmultipleNoSolid += 1;
			}
		}
		
		//_leon->_readWithAnchorMutationChoicesSize += 0.25;
		_bifurcations.push_back(nextNt);
		codeSeedNT(&_kmerModel, &kmer, nextNt, rightExtend);
		return kmer;
	}
	
	//return pos;
	
}

bool DnaEncoder::extendMutaChain(kmer_type kmer, int pos, bool rightExtend){
	return false;
	
	#ifdef PRINT_DEBUG_EXTMUTA
		//if(_lala == 2) return false;
	#endif
	
	_solidMutaChainStartPos = pos;
			
	#ifdef PRINT_DEBUG_EXTMUTA
		cout << "Extend muta chain !!!" << endl;
		cout << "\t" << _readseq << endl;
		cout << "\tStart kmer: " << kmer.toString(_kmerSize) << endl;
	#endif
	
	vector< vector< vector<kmer_type> > > mutaChains;
	vector< vector<kmer_type> > mutas;
	for(int i=0; i<4; i++) mutas.push_back(vector<kmer_type>());
	
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		
		if(_bloom->contains(mutatedKmerMin)){
			
			#ifdef PRINT_DEBUG_EXTMUTA
				cout << "\t\tAlternatives: " << mutatedKmer.toString(_kmerSize) << endl;
			#endif
			
			mutas[nt].push_back(mutatedKmer);
			
			//extendMutaChainRec(mutatedKmer, pos, rightExtend);
		}
	}
	
	mutaChains.push_back(mutas);
	
	return extendMutaChainRec(mutaChains, rightExtend);
}

bool DnaEncoder::extendMutaChainRec(vector< vector< vector<kmer_type> > >& mutaChains, bool rightExtend){
	int pos = mutaChains.size();
	_solidMutaChainSize = mutaChains.size();
	
	if(_solidMutaChainStartPos+pos >= _readSize || _solidMutaChainStartPos-pos < 0){
		#ifdef PRINT_DEBUG_EXTMUTA
			cout << "\t\t\t\t\tRead overflow" << endl;
		#endif
		_solidMutaChainSize += 1 ;
		return false;
	}
		
		
	vector< vector<kmer_type> > mutas;
	for(int i=0; i<4; i++) mutas.push_back(vector<kmer_type>());
	
	vector< vector<kmer_type> > mutatedKmers = mutaChains[pos-1];
	
	for(int NT=0; NT<4; NT++){
		
		for(int i=0; i < mutatedKmers[NT].size(); i++){
			kmer_type kmer = mutatedKmers[NT][i];
			
			for(int nt=0; nt<4; nt++){
				
				kmer_type mutatedKmer = kmer;
				codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
				kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
				
				if(_bloom->contains(mutatedKmerMin)){
					mutas[nt].push_back(mutatedKmer);
					
					#ifdef PRINT_DEBUG_EXTMUTA
						cout << "\t\t\tSuper alternatives: " << mutatedKmer.toString(_kmerSize) << endl;
					#endif
				}
			}
			
			
		}
	}
	mutaChains.push_back(mutas);
	
	_solidMutaChainSize = mutaChains.size();

	char originalNt;
	int originalNtBin;
	
	if(rightExtend)
		originalNt = _readseq[_solidMutaChainStartPos + pos];
	else
		originalNt = _readseq[_solidMutaChainStartPos - pos];
	
	originalNtBin = Leon::nt2bin(originalNt);
	
	#ifdef PRINT_DEBUG_EXTMUTA
		cout << "\t\t\t\tOriginal nt at: " << originalNt << endl;
	#endif
	
	
	if(mutaChains[pos][originalNtBin].size() == 0){
		#ifdef PRINT_DEBUG_EXTMUTA
			cout << "\t\t\t\t\tOriginal nt is not solid" << endl;
		#endif
		return false;
	}
		
	
	bool fullUniq = true;
	bool uniq = false;
	//mutas = mutaChains[mutaChains.size()-1];
	
	for(int nt=0; nt<4; nt++){
		
		if(mutas[nt].size() == 1){
			if(nt == originalNtBin){
				uniq = true;
			}
			else{
				fullUniq = false;
			}
		}
		
	}
		
	//cout << _sequence->getIndex() <<  "  Pos: " << _solidMutaChainStartPos + pos << endl;

	
	//cout << uniq << " " <<  fullUniq << endl;
	if(uniq && fullUniq){
		_solidMutaChainSize = pos;
		#ifdef PRINT_DEBUG_EXTMUTA
			cout << "\t\t\t\t\tMuta chain success" << endl;
			_lala += 1;
		#endif
		return true;
	}
	else{
		if(pos > 10){
			#ifdef PRINT_DEBUG_EXTMUTA
				cout << "\t\t\t\t\tToo much extend" << endl;
			#endif
			return false;
		}
		return extendMutaChainRec(mutaChains, rightExtend);
	}
	
	
	
}

int DnaEncoder::voteMutations(int pos, bool rightExtend){
	kmer_type kmer;
	int maxScore = 0;
	int votes[4];
	
	int bestNt;
	
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
	
	for(int nt=0; nt<4; nt++){
		
		kmer_type mutatedKmer = kmer;
		codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
		kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
		
		for(int j=0; j<4; j++){
			char nextNt;
			//int kmerPos;
			if(rightExtend){
				if(pos+1+_kmerSize+j >= _readSize) break;
				nextNt = _readseq[pos+1+_kmerSize+j];
			}
			else{
				if(pos-2-j < 0) break;
				nextNt = _readseq[pos-2-j];
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
	}
	
	/*
	cout << "---------------------" << endl;
	for(int i=0; i<mutations.size(); i++){
		int nt = mutations[i];
		cout << bin2NT[nt] << " " << votes[nt] << endl;
	}*/
	
	if(maxScore == 0){
		//cout << "No best NT" << endl;
		bestNt = -1;
	}
	//else
		//cout << "Best nt: " << bin2NT[bestNt] << endl;
	
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

DnaDecoder::~DnaDecoder(){
	//delete _rangeDecoder;
	//delete _outputFile;
	delete _inputFile;
	delete _anchorDictFile;
}

void DnaDecoder::setup(u_int64_t blockStartPos, u_int64_t blockSize, int sequenceCount){
	startBlock();
	_rangeDecoder.clear();
	
	_inputFile->seekg(blockStartPos, _inputFile->beg);
	_rangeDecoder.setInputFile(_inputFile);
	
	_blockStartPos = blockStartPos;
	_blockSize = blockSize;
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t-----------------------" << endl;
		cout << "\tDecoding block " << _blockStartPos << " - " << _blockStartPos+_blockSize << endl;
	#else
		_leon->_progress_decode->inc(1);
	#endif
	
	_sequenceCount = sequenceCount;
}

void DnaDecoder::execute(){
	
	//decodeFirstHeader();
		
	while(_processedSequenceCount < _sequenceCount){
		
		//cout << "lala" << endl;
	//int i=0;
	//while(i < Leon::READ_PER_BLOCK){
	//while(_inputFile->tellg() <= _blockStartPos+_blockSize){
		//if(_leon->_readCount > 1) return;
	
	
		u_int8_t readType = _rangeDecoder.nextByte(_readTypeModel);
		//cout << "Read type: " << (int)readType << endl;

		if(readType == 0)
			decodeAnchorRead(); //ici
		else if(readType == 1)
			decodeNoAnchorRead();
			
		endRead();
		//cout << _inputFile->tellg() << " " << _blockStartPos+_blockSize << endl;
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
	
	//cout << "endooo" << endl;
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
	deltaType = _rangeDecoder.nextByte(_readSizeDeltaTypeModel);
	deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _readSizeValueModel);
//	printf("read deltaValue %llu \n",deltaValue);
	
	_readSize = CompressionUtils::getValueFromDelta(deltaType, _prevReadSize, deltaValue);
	_prevReadSize = _readSize;
	
//	printf("read size %i \n",_readSize);
	
	//Decode anchor pos
	deltaType = _rangeDecoder.nextByte(_anchorPosDeltaTypeModel);
	deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosModel);
	int anchorPos = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorPos, deltaValue);
	_prevAnchorPos = anchorPos;
//	printf("anchor pos %i \n",anchorPos);

	
	//Decode anchor address
	deltaType = _rangeDecoder.nextByte(_anchorAddressDeltaTypeModel);
	deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _anchorAddressModel);
	u_int64_t anchorAddress = CompressionUtils::getValueFromDelta(deltaType, _prevAnchorAddress, deltaValue);
	_prevAnchorAddress = anchorAddress;
	
	kmer_type anchor = _leon->getAnchor(_anchorDictFile, anchorAddress); //laa
	
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
	_errorPos.clear();
	_Npos.clear();
	
	//Decode N pos
	u_int64_t NposCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(int i=0; i<NposCount; i++){
		
		deltaType = _rangeDecoder.nextByte(_NposDeltaTypeModel);
		deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _NposModel); //reprise
		u_int64_t nPos = CompressionUtils::getValueFromDelta(deltaType, _prevNpos, deltaValue);
		_Npos.push_back(nPos);
		_prevNpos = nPos;
		//_Npos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	
	//Decode error pos
	u_int64_t errorPosCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(int i=0; i<errorPosCount; i++){
		
		deltaType = _rangeDecoder.nextByte(_errorPosDeltaTypeModel);
		deltaValue = CompressionUtils::decodeNumeric(_rangeDecoder, _errorPosModel); //reprise
		u_int64_t errorPos = CompressionUtils::getValueFromDelta(deltaType, _prevErrorPos, deltaValue);
		_errorPos.push_back(errorPos);
		_prevErrorPos = errorPos;
		//_errorPos.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _anchorPosSizeModel, _anchorPosModel));
	}
	
	//Extend anchor to the left
	kmer_type kmer = anchor;
	for(int i=anchorPos-1; i>=0; i--){
		kmer = extendAnchor(kmer, i, false);
	}
	
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
}

kmer_type DnaDecoder::extendAnchor(kmer_type kmer, int pos, bool rightExtend){
	
	u_int8_t nextNt;
	kmer_type resultKmer;
		
	if(std::find(_Npos.begin(), _Npos.end(), pos) != _Npos.end()){
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
	
	if(std::find(_errorPos.begin(), _errorPos.end(), pos) != _errorPos.end()){
		nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_mutationModel));
		//cout << "tap 0     " << nextNt << endl;
		//_errorPos.erase(_errorPos.begin());
		if(rightExtend){
			_currentSeq += nextNt;
		}
		else{
			_currentSeq.insert(_currentSeq.begin(), nextNt);
		}
		//cout << _currentSeq << endl;
		//if(nextNt == 'N') nextN
		//resultKmer = kmer;
		
		for(int nt=0; nt<4; nt++){
			//if(nt == original_nt){
			//	continue;
			//}
			
			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
			kmer_type mutatedKmerMin = min(mutatedKmer, revcomp(mutatedKmer, _kmerSize));
			
			//mutatedKmer.printASCII(_kmerSize);
			
			if(_bloom->contains(mutatedKmerMin)){
				return mutatedKmer;
			}
			
		}
	
		//codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
		//return resultKmer;
	}
	
		
	//cout << kmer.toString(_kmerSize) << endl;
	kmer_type uniqKmer, mutatedSolidKmer;
	int uniqNt;
	//bool isKmerSolid = false;
	

	//kmer = _kmers[pos];

	int indexedKmerCount = 0;
	
	//cout << kmer.toString(_kmerSize) << endl;
	
	
	
	std::bitset<4> res4  = _bloom->contains4(kmer,rightExtend);
	for(int nt=0; nt<4; nt++){
		if(res4[nt]){
			kmer_type mutatedKmer = kmer;
			codeSeedBin(&_kmerModel, &mutatedKmer, nt, rightExtend);
			
			indexedKmerCount += 1;
			uniqNt = nt;
			uniqKmer = mutatedKmer;
			}
		}
		

	
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
		//cout << "case 1         " << nextNt << endl;
		resultKmer = uniqKmer;
	}
	else{
		nextNt = Leon::bin2nt(_rangeDecoder.nextByte(_mutationModel));
		//cout << "case 2          "<< nextNt << endl;
		resultKmer = kmer;
		codeSeedNT(&_kmerModel, &resultKmer, nextNt, rightExtend);
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
	
	_buffer += _currentSeq + '\n';
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\t\tRead: " << _currentSeq << endl;
	#endif
	//_outputFile->write(_currentSeq.c_str(), _currentSeq.size());
	_currentSeq.clear();
}




