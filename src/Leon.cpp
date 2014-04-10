//
//  Leon.cpp
//  leon
//
//  Created by Guillaume Rizk on 17/01/2014.
//  Copyright (c) 2014 rizk. All rights reserved.
//

/*
TODO
* Header coder:
* 	Checker les numeric en comparant leur valeur plutot que ses character
* 	Nombre de fields a réduire à la valeur de misIndex plutot qu'au véritable nombre de field
* 	Missing: encodage des zéro au début d'un numéric
* 	Remplacer la méthode strtoul par une methode string to u_int64_t dans CompressionUtils
* 
* Dna coder:
* 	Sécurité pour ne pas saturer la OAHash _anchorKmers, remplacer map par unordered_map
* 
* Optimisation:
* 	methode anchorExist: 2 acces a la map (key exist, puis operator [])
* 
* Bug:
* 	OAHash : faux positif ?
 * */
#include "Leon.hpp"
#include <DSK.hpp>

using namespace std;

#define SERIAL
#define PRINT_DEBUG



const u_int64_t ANCHOR_KMERS_HASH_SIZE = 500000000;
const char* Leon::STR_COMPRESS = "-c";
const char* Leon::STR_DECOMPRESS = "-d";
    




Leon::Leon () :
Tool("leon"),
_generalModel(256), _numericSizeModel(8),// _anchorKmers(ANCHOR_KMERS_HASH_SIZE),
_anchorDictModel(5) //5value: A, C, G, T, N
{
    //_kmerSize(27)
    /** We get an OptionsParser for DSK. */
    OptionsParser parserDSK = DSK::getOptionsParser();
    getParser()->add (parserDSK);
    
    getParser()->push_back (new OptionNoParam (Leon::STR_COMPRESS, "compress", false));
    getParser()->push_back (new OptionNoParam (Leon::STR_DECOMPRESS, "decompress", false));

    /** We add options specific to this tool. */

    
}


void Leon::execute()
{
    _solidFile = "tmp.solid"; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    bool compress = false;
    bool decompress = false;
    if(getParser()->saw (Leon::STR_COMPRESS)) compress = true;
    if(getParser()->saw (Leon::STR_DECOMPRESS)) decompress = true;
	if((compress && decompress) || (!compress && !decompress)){
		cout << "Choose one option among -c (compress) or -d (decompress)" << endl;
		return;
	}

	
	//setup global
	for(int i=0; i<8; i++){
		_numericModel.push_back(Order0Model(256));
	}
	
	if(compress)
		executeCompression();
	else
		executeDecompression();

	

        
    
    //outputFile->flush();

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    //getInfo()->add (1, "result");
    //getInfo()->add (2, "nb solid kmers in reads", "%ld", total_nb_solid_kmers_in_reads);
    
    if(!compress){
		delete _inputFile;
		delete _descInputFile;
	}
    
}

void Leon::createBloom ()
{
    printf("----- create bloom ----\n");
    TIME_INFO (getTimeInfo(), "fill bloom filter");
    
    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    NBITS_PER_KMER = 12;
    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_count));
    
    u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
    
    
    //printf("raw solidFileSize %llu fsize %llu    %lu %lu \n",System::file().getSize(_solidFile), solidFileSize,sizeof (kmer_type),sizeof (kmer_count));
    
    /** We create the kmers iterator from the solid file. */
    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
                                                                new IteratorFile<kmer_count>(_solidFile),
                                                                solidFileSize,
                                                                "fill bloom filter"
                                                                );
    LOCAL (itKmers);
    

	
    /** We instantiate the bloom object. */
    //BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    cout << "ESTIMATED:" << estimatedBloomSize << endl;
    _bloomSize = estimatedBloomSize;
    BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CacheCoherent,getInput()->getInt(STR_NB_CORES));
    _bloom = builder.build (itKmers);
}

void Leon::createKmerAbundanceHash(){
	cout << "\tFilling kmer abundance hash" << endl;
	
	vector<int> thresholds {200, 50, 20, 10};
	u_int64_t size = 0;
	u_int64_t maxSize = 500000000;
    _kmerAbundance = new OAHash<kmer_type>(maxSize);
    
    KmerModel model(_kmerSize);
    
    IteratorFile<kmer_count> it(_solidFile);
    
    for(int i=0; i<thresholds.size()-1; i++){
		for (it.first(); !it.isDone(); it.next()){
			kmer_count count = (*it);
			
			if(count.abundance >= thresholds[i]){
				if(i == 0 || count.abundance < thresholds[i-1]){
					_kmerAbundance->insert(count.value, count.abundance);
					//size += OAHash<kmer_type>::size_entry();
					size += sizeof(kmer_type) + sizeof(u_int32_t) + sizeof(int)*2;
				}
			}
	
			//if(size > maxSize-100000) break;
			if(size > maxSize) break;
		}
		
		if(size > maxSize) break;
	}

	cout << "\t\tNeeded memory: " << size << endl;
	cout << "\t\tAllocated memory: " << _kmerAbundance->memory_usage() << endl;
}














void Leon::executeCompression(){
	cout << "Start compression" << endl;
	
    _kmerSize = getInput()->getInt (STR_KMER_SIZE);
    //_solidFile = getInput()->getStr (STR_KMER_SOLID); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	_nks = getInput()->getInt("-nks");
    
	//Create and fill bloom
    createBloom ();
    LOCAL (_bloom);
	
    _inputFilename = getInput()->getStr(STR_URI_FILE);
    
    
    cout << "\tInput filename: " << _inputFilename << endl;
	
	u_int8_t infoByte = 0;
	
	//guess filename extension
	if(_inputFilename.find(".fa") !=  string::npos || _inputFilename.find(".fasta") !=  string::npos){
		cout << "\tInput format: Fasta" << endl;
		infoByte |= 0x01;
		_isFasta = true;
	}
	else if(_inputFilename.find(".fq") !=  string::npos || _inputFilename.find(".fastq") !=  string::npos){
		cout << "\tInput format: Fastq" << endl;
		_isFasta = false;
	} 
	else{
		cout << "\tUnknown input extension. Input extension must be one among fasta (.fa, .fasta) or fastq (.fq, .fastq)" << endl;
		return;
	}
	
	_rangeEncoder.encode(_generalModel, infoByte);
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _kmerSize);
	
    _inputBank = BankRegistery::singleton().getFactory()->createBank(_inputFilename);



    
    /*************************************************/
    // We create the modified file
    /*************************************************/
    
    
    string dir = System::file().getDirectory(_inputFilename);
    string prefix = System::file().getBaseName(_inputFilename);
    _outputFilename = dir + "/" + prefix + ".txt"; //".leon"
	_outputFile = System::file().newFile(_outputFilename, "wb");
    cout << "\tOutput filename: " << _outputFilename << endl;


    u_int64_t total_nb_solid_kmers_in_reads = 0;
    int nb_threads_living;
    _nb_cores = getInput()->getInt(STR_NB_CORES);
    
    //Compression
	startHeaderCompression();
	startDnaCompression();
	
	endCompression();
}
		
void Leon::writeBlock(u_int8_t* data, u_int64_t size){
	if(size <= 0) return;
	
	//cout << "\t-----------------------" << endl;
	//cout << "\tWrite block " << _blockSizes.size() << endl;
	//cout << "\tSequence " << encoder->_lastSequenceIndex-READ_PER_BLOCK << " - " << encoder->_lastSequenceIndex << endl;
	//cout << "Thread id: " << thread_id << endl;
	//cout << "\tEncoded size (byte): " << size << endl;
	
	_outputFile->fwrite(data, size, 1);
	//int thread_id = encoder->getId();
	
	_blockSizes.push_back(size);
		
	/*
	int thread_id = encoder->getId();
	cout << "Incoming thread: " << thread_id << endl;
	pthread_mutex_lock(&writer_mutex);
	while(thread_id > _next_writer_thread_id){
		thread_id = encoder->getId();
		cout << thread_id << " " << _next_writer_thread_id << endl;
		pthread_cond_wait(&buffer_full_cond, &writer_mutex);
	}
	pthread_mutex_unlock(&writer_mutex);
	
	cout << "-----------------------" << endl;
	cout << "Thread id: " << thread_id << endl;
	cout << "Buffer size: " << encoder->getBufferSize() << endl;
	
	_outputFile->fwrite(encoder->getBuffer(), encoder->getBufferSize(), 1);
	
	pthread_mutex_lock(&writer_mutex);
	cout << encoder->_lastSequenceIndex << endl;
	_next_writer_thread_id = thread_id + 1;
	if(_next_writer_thread_id == _nb_cores){
		_next_writer_thread_id = 0;
	}
	pthread_mutex_unlock(&writer_mutex);
	
	pthread_cond_broadcast(&buffer_full_cond);*/
	
}
		
void Leon::endCompression(){
	_rangeEncoder.flush();
	_outputFile->fwrite(_rangeEncoder.getBuffer(true), _rangeEncoder.getBufferSize(), 1);
	_outputFile->flush();
	
	u_int64_t inputFileSize = System::file().getSize(_inputFilename.c_str());
	cout << endl;
	cout << "\tInput: " << endl;
	cout << "\t\tFilename: " << _inputFilename << endl;
	cout << "\t\tSize: " << inputFileSize << endl;
	
	u_int64_t outputFileSize = System::file().getSize(_outputFilename);
	cout << "\tOutput: " << endl;
	cout << "\t\tFilename: " << _outputFilename << endl;
	cout << "\t\tSize: " << outputFileSize << endl;
	cout << "\tCompression rate: " << (float)((double)outputFileSize / (double)inputFileSize) << endl;
	cout << "\t\tHeader: " << (float)_headerCompRate << endl;
	cout << "\t\tDna: " << (float)_dnaCompRate << endl;
}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
void Leon::startHeaderCompression(){
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          _inputBank->estimateNbSequences(),
                                                          "Compressing headers"
                                                          );
    LOCAL(itSeq);
    
    
	_totalHeaderSize = 0;
	_totalHeaderCompressedSize = 0;
	
    cout << endl << "Start header compression" << endl;
    
    //write first header to file and store it in _firstHeader variable
	ifstream inputFileTemp(getInput()->getStr(STR_URI_FILE), ios::in);
	getline(inputFileTemp, _firstHeader);
	inputFileTemp.close();
	_firstHeader.erase(_firstHeader.begin());
	cout << "\tFirst Header: " << _firstHeader << endl;
	cout << "\tSize: " << _firstHeader.size() << endl;
	_totalHeaderSize += _firstHeader.size();
	
	//encode the size of the first header on 2 byte and the header itself
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _firstHeader.size());
	for(int i=0; i < _firstHeader.size(); i++){
		_rangeEncoder.encode(_generalModel, _firstHeader[i]);
	}
	_rangeEncoder.flush();
	_totalHeaderCompressedSize += _rangeEncoder.getBufferSize();
	_outputFile->fwrite(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), 1);
	_rangeEncoder.clear();
	
	//cout << "Block start pos: " << _outputFile->tell() << endl;
	
	//iterate on read sequences and compress headers
	TIME_INFO (getTimeInfo(), "header compression");

	//int nb_threads_living = 0 ;
	
	#ifdef SERIAL
		setDispatcher (new SerialDispatcher());
	#else
		setDispatcher (  new Dispatcher (_nb_cores) );
	#endif

	//getDispatcher()->iterate (itSeq,  HeaderEncoder(this, &nb_threads_living), 10000);
	getDispatcher()->iterate (itSeq,  HeaderEncoder(this), READ_PER_BLOCK);
	endHeaderCompression();
}

void Leon::endHeaderCompression(){
	//u_int64_t descriptionStartPos = _outputFile->tell();
	//cout << "Description start pos: " << descriptionStartPos << endl;
	
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _blockSizes.size());
	for(int i=0; i<_blockSizes.size(); i++){
		//cout << "block size: " << _blockSizes[i] << endl;
		CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _blockSizes[i]);
	}
	
	_blockSizes.clear();
	//Encode description start position in the output file
	//CompressionUtils::encodeFixedNumeric(_rangeEncoder, _generalModel, descriptionStartPos, 8);
	
	
	
	//_totalHeaderCompressedSize += _rangeEncoder.getBufferSize();

	
	//string descStr = to_string(descriptionStartPos);
	//descStr += (u_int8_t)descStr.size();
	//cout << descStr.size() << " " << descStr << endl;
	//_outputFile->fwrite(descStr.c_str(), descStr.size(), 1);
	_headerCompRate = ((double)_totalHeaderCompressedSize / _totalHeaderSize);
	cout << "End header compression" << endl;
	cout << "\tBlock count: " << _blockSizes.size() << endl;
	//cout << "\tBlock data size: " << _rangeEncoder.getBufferSize() << endl;
	cout << "\tTotal header size: " << _totalHeaderSize << endl;
	cout << "\tTotal header compressed size: " << _totalHeaderCompressedSize << endl;
	cout << "\tCompression rate: " << (float)(_headerCompRate) << endl;
	
	//_rangeEncoder.clear();
}



































void Leon::startDnaCompression(){
    cout << endl << "Start dna compression" << endl;
    
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          _inputBank->estimateNbSequences(),
                                                          "Compressing dna"
                                                          );
    LOCAL(itSeq);
    
	_anchorAdress = 0;
	_totalDnaSize = 0;
	_totalDnaCompressedSize = 0;
	_realDnaCompressedSize = 0;
	
	createKmerAbundanceHash();
    
	//iterate on read sequences and compress headers
	TIME_INFO (getTimeInfo(), "DNA compression");

	//int nb_threads_living = 0 ;
	
	#ifdef SERIAL
		setDispatcher (new SerialDispatcher());
	#else
		setDispatcher (  new Dispatcher (_nb_cores) );
	#endif

	//getDispatcher()->iterate (itSeq,  HeaderEncoder(this, &nb_threads_living), 10000);
	getDispatcher()->iterate (itSeq,  DnaEncoder(this), READ_PER_BLOCK);
	endDnaCompression();
	
	//auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (argv[1])); 
}

void Leon::endDnaCompression(){
	/*
		u_int64_t _noAnchor_full_N_kmer_count;
		u_int64_t _noAnchor_with_N_kmer_count;
		u_int64_t _readCount;
		double _anchorKmerCount;
		double _readWithoutAnchorCount;
		double _total_kmer_indexed;
		double _uniq_mutated_kmer;
		u_int64_t _total_kmer;
		*/
	cout << endl;
	cout << endl;
	cout << System::file().getBaseName(getInput()->getStr(STR_URI_FILE)) << endl;
	cout << "\tRead count: " << _readCount << endl;
	cout << "\tDna size: " << _totalDnaSize << endl;
	cout << "\tRead per anchor: " << _readCount / _anchorKmerCount << endl;
	cout << "\tBit per anchor: " << log2(_anchorKmerCount);
	cout << endl;
	cout << "\tRead without anchor: " << (_readWithoutAnchorCount*100) / _readCount << "%" << endl;
	cout << "\t\tWith N: " << (_noAnchor_with_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	cout << "\t\tFull N: " << (_noAnchor_full_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	//cout << "\tAnchor hash size: " << _anchorKmerCount*(sizeof(kmer_type)+sizeof(int)) << endl;
	//cout << "Total kmer: " <<  _leon->_total_kmer << endl;
	//cout << "Kmer in bloom:  " << (_leon->_total_kmer_indexed*100) / _leon->_total_kmer << "%" << endl;
	//cout << "Uniq mutated kmer: " << (_leon->_uniq_mutated_kmer*100) / _leon->_total_kmer << "%" << endl;
	//cout << "anchor Hash:   memory_usage: " << _anchorKmers.memory_usage() << "    memory_needed: " << _anchorKmers.size_entry ()*_anchorKmerCount << " (" << (_anchorKmers.size_entry ()*_leon->_anchorKmerCount*100) / _anchorKmers.memory_usage() << "%)" << endl;
	
	cout << endl;
	u_int64_t readWithAnchorCount = _readCount - _readWithoutAnchorCount;
	
	//@ anchor kmer
	_readWithAnchorSize += ((double)readWithAnchorCount * log2(_anchorKmerCount)) / 8;
	_readWithAnchorSize += readWithAnchorCount*2;
	_readWithAnchorSize += _readWithAnchorMutationChoicesSize;
	
	//anchor dict
	_totalDnaCompressedSize += _anchorKmerCount*sizeof(kmer_type);
	//bloom
    _totalDnaCompressedSize += _bloom->getSize();
	//read without anchor (read_size * 0.375) 3 bit/nt
	_totalDnaCompressedSize += _readWithoutAnchorSize;
	_totalDnaCompressedSize += _readWithAnchorSize;
	
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _blockSizes.size());
	for(int i=0; i<_blockSizes.size(); i++){
		//cout << "block size: " << _blockSizes[i] << endl;
		CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _blockSizes[i]);
	}
	
	writeBloom();
	writeAnchorDict();
	
	_dnaCompRate = ((double)_realDnaCompressedSize / _totalDnaSize);
	cout << "\tCompression rate: " << ((double)_totalDnaCompressedSize / _totalDnaSize) << "    " << (float)_dnaCompRate << endl;
	cout << "\t\tBloom: " << ((_bloom->getSize()*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\t\tAnchor kmer dictionnary: " << ((_anchorKmerCount*sizeof(kmer_type)*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\t\tRead with anchor: " << ((_readWithAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\t\t\t@anchor: " << (((double)readWithAnchorCount * log2(_anchorKmerCount)*100 / 8) / (double)_readWithAnchorSize) << "%" << endl;
	cout << "\t\t\tMutations choices: " << ((_readWithAnchorMutationChoicesSize*100) / (double)_readWithAnchorSize)  << endl;
	cout << "\t\t\tOther: " << ((readWithAnchorCount*2*100)/(double)_readWithAnchorSize) << "%" << endl;
	cout << "\t\tRead without anchor: " << ((_readWithoutAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << endl;
	cout << "\tTotal mutations: " << _MCtotal << endl;
	cout << "\t\tUniq origNT solid: " << ((_MCuniqSolid*100)/_MCtotal) << endl;
	cout << "\t\tMultiple origNT solid: " << ((_MCmultipleSolid*100)/_MCtotal) << endl;
	cout << "\t\tNo alternative: " << ((_MCnoAternative*100)/_MCtotal) << endl;
	cout << "\t\tUniq origNT no solid: " << ((_MCuniqNoSolid*100)/_MCtotal) << endl;
	cout << "\t\tMultiple origNT no solid: " << ((_MCmultipleNoSolid*100)/_MCtotal) << endl;
	
}

void Leon::writeBloom(){
	//_bloom->save(_outputFilename + ".bloom.temp");
	
	_realDnaCompressedSize += _bloom->getSize();
	//_outputFile->fwrite(_anchorRangeEncoder.getBuffer(), size, 1);

	u_int64_t size = _bloom->getSize();
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, size);
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, _bloomSize);
	
	//u_int8_t*& bloomData = _bloom->getArray();
	_outputFile->fwrite(_bloom->getArray(), size, 1);
}

void Leon::writeAnchorDict(){
	_anchorRangeEncoder.flush();
	_realDnaCompressedSize += _anchorRangeEncoder.getBufferSize();
	
	u_int64_t size = _anchorRangeEncoder.getBufferSize();
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, size);
	
	//cout << "Anchor dict size: " << size << endl;
	//cout << "\t pos: " << _outputFile->tell() << endl;
	//cout << "count: " << _anchorKmerCount << endl;
	
	_outputFile->fwrite(_anchorRangeEncoder.getBuffer(), size, 1);
	

}

bool Leon::anchorExist(const kmer_type& kmer, u_int32_t* anchorAdress){
	if(_anchorKmers.find( kmer ) != _anchorKmers.end()){ 
		*anchorAdress = _anchorKmers[kmer];
		return true;
	}
	return false;
	//return _anchorKmers.get(kmer, (int*)anchorAdress); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!changer OAHash in to u_int32_t
	//u_int32_t anchorAdress;
	//bool exist = _anchorKmers.get(kmer, &anchorAdress);
	//return _anchorKmers.has_key(kmer);
}

int Leon::findAndInsertAnchor(const vector<kmer_type>& kmers, u_int32_t* anchorAdress){
	

		
	//cout << "\tSearching and insert anchor" << endl;
	int maxAbundance = -1;
	int bestPos;
	kmer_type bestKmer;
	//int bestPos = -1;
	

	kmer_type kmer, kmerMin;
	
	/*
	////////////
	for(int i=0; i<kmers.size(); i++){
		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		if(_bloom->contains(kmerMin)){
			encodeInsertedAnchor(kmerMin);
			_anchorKmers.insert(kmerMin, _anchorAdress);
			*anchorAdress = _anchorAdress;
			_anchorKmerCount += 1;
			_anchorAdress += 1;
			return i;
		}
	}
	return -1;
	/////////////////////*/
	
	for(int i=0; i<kmers.size(); i++){
		kmer = kmers[i];
		kmerMin = min(kmer, revcomp(kmer, _kmerSize));
		
		int abundance;
		if(_kmerAbundance->get(kmerMin, &abundance)){
			if(abundance > maxAbundance){
				maxAbundance = abundance;
				bestKmer = kmerMin;
				bestPos = i;
				//cout << bestPos << " " << abundance << " " << kmer.toString(_kmerSize) << " " << revcomp(kmer, _kmerSize).toString(_kmerSize) << endl;
			}
			//cout << abundance << endl;
		}
		else if(maxAbundance == -1 && _bloom->contains(kmerMin)){
			maxAbundance = _nks;
			bestKmer = kmerMin;
			bestPos = i;
		}
		/*
		if(_bloom->contains(kmerMin)){
			encodeInsertedAnchor(kmerMin);
			_anchorKmers.insert(kmerMin, _anchorAdress);
			*anchorAdress = _anchorAdress;
			_anchorKmerCount += 1;
			_anchorAdress += 1;
			return i;
		}
		* */
	}
	
	if(maxAbundance == -1)
		return -1;

	encodeInsertedAnchor(bestKmer);
	
	_anchorKmers[bestKmer] = _anchorAdress;
	//_anchorKmers.insert(bestKmer, _anchorAdress);
	
	
	
	*anchorAdress = _anchorAdress;
	_anchorKmerCount += 1;
	_anchorAdress += 1;
	return bestPos;
}

void Leon::encodeInsertedAnchor(const kmer_type& kmer){
	//static int i = 0;
	
	string kmerStr = kmer.toString(_kmerSize);
	//if(i<10) cout << "\t\t" << kmerStr << endl;
	for(int i=0; i<kmerStr.size(); i++){
		_anchorRangeEncoder.encode(_anchorDictModel, Leon::nt2bin(kmerStr[i]));
	}
	
	//i+= 1;
}




















void Leon::executeDecompression(){

	cout << "Start decompression" << endl;
    string inputFilename = getInput()->getStr(STR_URI_FILE);
    //string inputFilename = prefix + ".txt"; //".leon"
	//_outputFile = System::file().newFile(outputFilename, "wb");
	cout << "\tInput filename: " << inputFilename << endl;
	
	_descInputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	//outputFile.open(outputFilenameDecompressed.c_str(), ios::out);
		
	_rangeDecoder.setInputFile(_descInputFile);
	
	//Decode the first byte of the compressed file which is an info byte
	u_int8_t infoByte = _rangeDecoder.nextByte(_generalModel);
	//the first bit holds the file format. 0: fastq, 1: fasta
	_isFasta = ((infoByte & 0x01) == 0x01);
	if(_isFasta)
		cout << "\tInput format: Fasta" << endl;
	else
		cout << "\tInput format: Fastq" << endl;
	
	
	
	
	////////// test bank
	string dir = System::file().getDirectory(inputFilename);
    string prefix = System::file().getBaseName(inputFilename);
	if(_isFasta)
		_testBank = BankRegistery::singleton().getFactory()->createBank(dir + "/" + prefix + ".fasta");
	else
		_testBank = BankRegistery::singleton().getFactory()->createBank(dir + "/" + prefix + ".fastq");
    //_testBank = BankRegistery::singleton().getFactory()->createBank("/local/gbenoit/data/raw/ecoli/SRR765913.fasta");
    //_testBank = BankRegistery::singleton().getFactory()->createBank("/local/gbenoit/data/raw/ecoli/SRR747737.fasta");
	_testBankIt = _testBank->iterator();
	_testBankIt->first();
	//cout << "lala     " <<  << endl;
	//_testBankIt->next();
	//cout << "lala     " << (*_testBankIt)->getComment() << endl;
	//////////
	
	
	
	
	//Get kmer size
	_kmerSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	cout << "\tKmer size: " << _kmerSize << endl;
	
	//Decode the first header
	u_int16_t firstHeaderSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	for(int i=0; i<firstHeaderSize; i++){
		_firstHeader += _rangeDecoder.nextByte(_generalModel);
	}
	cout << _firstHeader << endl;
	_filePos = _descInputFile->tellg();
	//cout << "Block start pos: " << blockStartPos << endl;
	
	//Go to the end of the file to decode blocks informations, data are read in reversed order (from right to left in the file)
	//The first number is the number of data blocks
	_descInputFile->seekg(0, _descInputFile->end);
	_rangeDecoder.setInputFile(_descInputFile, true);

	//HEADER DECOMPRESSION
	cout << "Start Header decompression" << endl;
	setupNextComponent();
	
	HeaderDecoder headerDecoder(this, _inputFile, NULL);
	for(int i=0; i<_blockSizes.size(); i++){
		u_int64_t blockSize = _blockSizes[i];
		//headerDecoder.setup(_filePos, blockSize);
		_filePos += blockSize;
	}
	
	//DNA DECOMPRESSION
	cout << "Start DNA decompression" << endl;
	setupNextComponent();

	decodeBloom();
	decodeAnchorDict();
	
	
	_inputFile->seekg(_filePos, _inputFile->beg);
	DnaDecoder dnaDecoder(this, _inputFile, NULL);
	for(int i=0; i<_blockSizes.size(); i++){
		//cout << "lala" << endl;
	//for(int i=0; i<1; i++){
		u_int64_t blockSize = _blockSizes[i];
		//dnaDecoder.setup(_filePos, blockSize);
		_filePos += blockSize;
		//cout << _readCount << endl;
	}
	
}


void Leon::setupNextComponent(){
	//Go to the block position
	_inputFile->seekg(_filePos, _inputFile->beg);
	
	_blockSizes.clear();
	//u_int64_t size = 0;
	
	u_int64_t blockCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	for(int i=0; i<blockCount; i++){
		u_int64_t blockSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
		_blockSizes.push_back(blockSize);
		//size += blockSize;
	}
	
	cout << "block count: " << blockCount << endl;
	for(int i=0; i<_blockSizes.size(); i++){
		cout << _blockSizes[i] << " ";
	}
	cout << endl;
	
}

void Leon::decodeBloom(){
	cout << "\tDecode bloom filter" << endl;

	u_int64_t dictPos = _inputFile->tellg();
	for(int i=0; i<_blockSizes.size(); i++){
		dictPos += _blockSizes[i];
	}
	//cout << "Anchor dict pos: " << dictPos << endl;
	
	_inputFile->seekg(dictPos, _inputFile->beg);
	
	
	//createBloom();
	u_int64_t bloomSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	
	char* buffer = new char[bloomSize];
    _inputFile->read(buffer, bloomSize);
    
    cout << bloomSize << endl;
    
	u_int64_t tai = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	//u_int64_t tai = (bloomSize) * 8LL;
	
	cout << tai << endl;
	
	//_bloom->contains("lala");
	//nchar  = (1+tai/8LL);
	_bloom = new BloomCacheCoherent<kmer_type>(tai, 7);
	
	memcpy(_bloom, buffer, _bloom->getSize());
	//u_int8_t*& bloomData = _bloom->getArray();
	//_outputFile->fwrite(bloomData, _bloom->getSize(), 1);
	
	_inputFile->seekg(bloomSize, _inputFile->cur);
}

void Leon::decodeAnchorDict(){
	cout << "\tDecode anchor dict" << endl;
	
	u_int64_t anchorDictSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	//cout << anchorDictSize << endl;
	

	_anchorRangeDecoder.setInputFile(_inputFile);
	string anchorKmer = "";
	//int lala =0;
	
	//_vecAnchorKmers.push_back(kmer_type()); //Insert a random kmer at first pos because adresses start at 1 (OAHash)
	u_int64_t dictPos = _inputFile->tellg();
	
	KmerModel model(_kmerSize, KMER_DIRECT);
	//int i=0;
	//cout << _inputFile->tellg()  << " " << dictPos+anchorDictSize << endl;
	while(_inputFile->tellg() <= dictPos+anchorDictSize){
		u_int8_t c = _anchorRangeDecoder.nextByte(_anchorDictModel);
		anchorKmer += Leon::bin2nt(c);
		if(anchorKmer.size() == _kmerSize){
			
			
			//if(i<=10) cout << anchorKmer << endl;
			//cout << "1: " << anchorKmer << endl;
			kmer_type kmer = model.codeSeed(anchorKmer.c_str(), Data::ASCII);
			//cout << "2: " << model.toString(kmer) << endl;
			//lala += 1;
			_vecAnchorKmers.push_back(kmer);
			
			anchorKmer.clear();
			//i++;
		}
	}
	
	
	cout << "\t\tAnchor count: " << _vecAnchorKmers.size() << endl;
}



kmer_type Leon::getAnchor(u_int32_t adress){
	return _vecAnchorKmers[adress];
}









