
#include "DnaCoder.hpp"

const u_int64_t ANCHOR_KMERS_HASH_SIZE = 10000000;
/*
#define PRINT_DEBUG_ENCODER

#define PRINT_DEBUG_DECODER
*/
		
//====================================================================================
// ** AbstractDnaCoder
//====================================================================================
AbstractDnaCoder::AbstractDnaCoder(Leon* leon)
{
	_leon = leon;
	_bloom = _leon->_bloom;
}


//====================================================================================
// ** DnaEncoder
//====================================================================================
DnaEncoder::DnaEncoder(Leon* leon) :
AbstractDnaCoder(leon), _anchorKmers(ANCHOR_KMERS_HASH_SIZE)
{
}

DnaEncoder::DnaEncoder(const DnaEncoder& copy) :
AbstractDnaCoder(copy._leon), _anchorKmers(ANCHOR_KMERS_HASH_SIZE)
{
	//_leon = copy._leon;
	//_bloom = copy._bloom;
}

DnaEncoder::~DnaEncoder(){
	writeBlock();
}

void DnaEncoder::operator()(Sequence& sequence){
	_sequence = &sequence;
	//cout << _sequence->getIndex() << endl;
	_readSize = _sequence->getDataSize();
	_readseq = _sequence->getDataBuffer();
	_kmerSize = _leon->_kmerSize;
	
	_leon->_totalDnaSize += _readSize; //unsynch
	
	
	//_lastSequenceIndex = sequence->getIndex();
	
	if(_sequence->getIndex() % Leon::READ_PER_BLOCK == 0){
		writeBlock();
		startBlock();
	}
	
	execute();
}

void DnaEncoder::writeBlock(){
	/*
	if(_rangeEncoder.getBufferSize() > 0){
		_rangeEncoder.flush();
	}
	
	_leon->_totalDnaCompressedSize += _rangeEncoder.getBufferSize();
	_leon->writeBlock(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize());
	_rangeEncoder.clear();
	*/
	
	cout << "------------------------" << endl;
	cout << "Read count:    " << _leon->_readCount << endl;
	cout << "Anchor kmer count:   " << _leon->_anchorKmerCount << endl;
	cout << "Read per anchor:    " << _leon->_readCount / _leon->_anchorKmerCount << endl;
	cout << "Bit per anchor: " << log2(_leon->_anchorKmerCount);
	cout << endl;
	cout << "Read without anchor: " << (_leon->_readWithoutAnchorCount*100) / _leon->_readCount << "%" << endl;
	cout << "\tWith N: " << (_leon->_noAnchor_with_N_kmer_count*100) / _leon->_readWithoutAnchorCount << "%" << endl;
	cout << "\tFull N: " << (_leon->_noAnchor_full_N_kmer_count*100) / _leon->_readWithoutAnchorCount << "%" << endl;
	cout << endl;
	//cout << "Total kmer: " <<  _leon->_total_kmer << endl;
	//cout << "Kmer in bloom:  " << (_leon->_total_kmer_indexed*100) / _leon->_total_kmer << "%" << endl;
	//cout << "Uniq mutated kmer: " << (_leon->_uniq_mutated_kmer*100) / _leon->_total_kmer << "%" << endl;
	cout << "anchor Hash:   memory_usage: " << _anchorKmers.memory_usage() << "    memory_needed: " << _anchorKmers.size_entry ()*_leon->_anchorKmerCount << " (" << (_anchorKmers.size_entry ()*_leon->_anchorKmerCount*100) / _anchorKmers.memory_usage() << "%)" << endl;
	
}


void DnaEncoder::startBlock(){
}

void DnaEncoder::execute(){
	_leon->_readCount += 1;
	
	bool anchorExist = false;
	kmer_type anchorKmer = 0;
	kmer_type current_kmer, current_kmer_min;
	KmerModel model (_kmerSize, KMER_DIRECT);
	KmerModel::Iterator itKmer (model);
	
	int i=0;
	itKmer.setData(_sequence->getData());
	for (itKmer.first(); !itKmer.isDone(); itKmer.next(), i++){
		current_kmer = *itKmer;
		current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		if(_anchorKmers.has_key(current_kmer_min)){
			anchorKmer = current_kmer_min;
			anchorExist = true;
			break;
		}
	}
	
	i=0;
	if(!anchorExist){
		itKmer.setData (_sequence->getData());
		for (itKmer.first(); !itKmer.isDone(); itKmer.next(), i++){
			current_kmer = *itKmer;
			current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
			if(_bloom->contains(current_kmer_min)){
				anchorKmer = current_kmer_min;
				_anchorKmers.insert(anchorKmer, 1);
				_leon->_anchorKmerCount += 1;
				anchorExist = true;
				break;
			}
		}
	}
		
	if(!anchorExist){
		_leon->_readWithoutAnchorSize += _readSize*0.375;
		_leon->_readWithoutAnchorCount += 1;
		
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
		}
		
		return;
	}
	
	/*
	i=0;
	itKmer.setData (_sequence->getData());
	for (itKmer.first(); !itKmer.isDone(); itKmer.next(), i++){
		if(i >= _readSize-_kmerSize) continue;
		
		_bloocoo->_total_kmer +=1;
		
		current_kmer = *itKmer;
		current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		
		if(_bloom->contains(current_kmer_min)){
			_bloocoo->_total_kmer_indexed += 1;
		}
		else{
			continue;
		}
			
		
		//cout << i  << " ";
		//current_kmer.printASCII(_kmerSize);
		//cout << current_kmer << endl;
		//kmer_type current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		
		//int original_nt = (_readseq[i]>>1)&3;
		int indexedKmerCount = 0;
		
		for(int nt=0; nt<4; nt++)
		{
			//if(nt == original_nt){
			//	continue;
			//}
			
			kmer_type mutated_kmer = current_kmer;
			codeSeedBin(&model, &mutated_kmer, nt, RIGHT);
			kmer_type mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
			
			//mutated_kmer.printASCII(_kmerSize);
			
			if(_bloom->contains(mutated_kmer_min)){
				indexedKmerCount += 1;
			}
			
		}
		
		if(indexedKmerCount == 1){
			_bloocoo->_uniq_mutated_kmer += 1;
		}
	}*/
	
}
