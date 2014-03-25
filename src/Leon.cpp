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
 * */
#include "Leon.hpp"
#include <DSK.hpp>

using namespace std;

#define SERIAL
#define PRINT_DEBUG



const char* Leon::STR_COMPRESS = "-c";
const char* Leon::STR_DECOMPRESS = "-d";
    




Leon::Leon () : Tool("leon"), _generalModel(256), _headerNumericSizeModel(10)
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
		_headerNumericModels.push_back(Order0Model(256));
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
    
    
    delete _inputFile;
    
}

Bloom<kmer_type>* Leon::createBloom ()
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
                                                                new IteratorFile<kmer_count> (_solidFile),
                                                                solidFileSize,
                                                                "fill bloom filter"
                                                                );
    LOCAL (itKmers);
    
    /** We instantiate the bloom object. */
    BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CacheCoherent,getInput()->getInt(STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build (itKmers);

    /** We return the created bloom filter. */
    return bloom;
}

void Leon::executeCompression(){
	cout << "Start compression" << endl;
	
    _kmerSize           = getInput()->getInt (STR_KMER_SIZE);
    _solidFile          = getInput()->getStr (STR_KMER_SOLID);
    
	//Create and fill bloom
    _bloom = createBloom ();
    LOCAL (_bloom);

    string inputFilename = getInput()->getStr(STR_URI_FILE);
    //string inputFilename = "/local/gbenoit/data/illumina/SRR959239.fasta";
    
    
    cout << "\tInput filename: " << inputFilename << endl;
	
	u_int8_t infoByte = 0;
	
	//guess filename extension
	if(inputFilename.find(".fa") !=  string::npos || inputFilename.find(".fasta") !=  string::npos){
		cout << "\tInput format: Fasta" << endl;
		infoByte |= 0x01;
		_isFasta = true;
	}
	else if(inputFilename.find(".fq") !=  string::npos || inputFilename.find(".fastq") !=  string::npos){
		cout << "\tInput format: Fastq" << endl;
		_isFasta = false;
	} 
	else{
		cout << "\tUnknown input extension. Input extension must be one among fasta (.fa, .fasta) or fastq (.fq, .fastq)" << endl;
		return;
	}
	
	_rangeEncoder.encode(_generalModel, infoByte);
	
    IBank* inbank = BankRegistery::singleton().getFactory()->createBank(inputFilename);


    /*************************************************/
    // We create a sequence iterator for the bank
    /*************************************************/
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          inbank->iterator(),
                                                          inbank->estimateNbSequences(),
                                                          "Compressing headers"
                                                          );
    LOCAL (itSeq);
    
    /*************************************************/
    // We create the modified file
    /*************************************************/
    
    
    string dir = System::file().getDirectory(inputFilename);
    string prefix = System::file().getBaseName(inputFilename);
    string outputFilename = dir + "/" + prefix + ".txt"; //".leon"
	_outputFile = System::file().newFile(outputFilename, "wb");
    cout << "\tOutput filename: " << outputFilename << endl;


    u_int64_t total_nb_solid_kmers_in_reads = 0;
    int nb_threads_living;
    _nb_cores = getInput()->getInt(STR_NB_CORES);
    
    //Compression
	//startHeaderCompression(itSeq);
	startDnaCompression(itSeq);
}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
void Leon::startHeaderCompression(Iterator<Sequence>* itSeq){
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
	CompressionUtils::encodeFixedNumeric(_rangeEncoder, _headerNumericModels, _firstHeader.size(), 2);
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
	u_int64_t descriptionStartPos = _outputFile->tell();
	cout << "Description start pos: " << descriptionStartPos << endl;
	
	_rangeEncoder.clear();
	

	CompressionUtils::encodeNumeric(_rangeEncoder, _headerNumericSizeModel, _headerNumericModels, _blockSizes.size());
	for(int i=0; i<_blockSizes.size(); i++){
		CompressionUtils::encodeNumeric(_rangeEncoder, _headerNumericSizeModel, _headerNumericModels, _blockSizes[i]);
	}
	
	//Encode description start position in the output file
	//CompressionUtils::encodeFixedNumeric(_rangeEncoder, _generalModel, descriptionStartPos, 8);
	
	_rangeEncoder.flush();
	
	//_totalHeaderCompressedSize += _rangeEncoder.getBufferSize();
	
	_outputFile->fwrite(_rangeEncoder.getBuffer(true), _rangeEncoder.getBufferSize(), 1);
	
	//string descStr = to_string(descriptionStartPos);
	//descStr += (u_int8_t)descStr.size();
	//cout << descStr.size() << " " << descStr << endl;
	//_outputFile->fwrite(descStr.c_str(), descStr.size(), 1);
	
	cout << "End header compression" << endl;
	cout << "\tBlock count: " << _blockSizes.size() << endl;
	cout << "\tBlock data size: " << _rangeEncoder.getBufferSize() << endl;
	cout << "\tTotal header size: " << _totalHeaderSize << endl;
	cout << "\tTotal header compressed size: " << _totalHeaderCompressedSize << endl;
	cout << "\tCompression rate: " << (float)((double)_totalHeaderCompressedSize / _totalHeaderSize) << endl;
	
	_rangeEncoder.clear();
}

void Leon::writeBlock(u_int8_t* data, u_int64_t size){
	if(size <= 0) return;
	
	cout << "\t-----------------------" << endl;
	cout << "\tWrite block " << _blockSizes.size() << endl;
	//cout << "\tSequence " << encoder->_lastSequenceIndex-READ_PER_BLOCK << " - " << encoder->_lastSequenceIndex << endl;
	//cout << "Thread id: " << thread_id << endl;
	cout << "\tEncoded size (byte): " << size << endl;
	
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

































void Leon::startDnaCompression(Iterator<Sequence>* itSeq){
	_totalDnaSize = 0;
	_totalDnaCompressedSize = 0;
	
    cout << endl << "Start dna compression" << endl;
    
    
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
	cout << "------------------------" << endl;
	cout << "Read count:    " << _readCount << endl;
	cout << "Anchor kmer count:   " << _anchorKmerCount << endl;
	cout << "Read per anchor:    " << _readCount / _anchorKmerCount << endl;
	cout << "Bit per anchor: " << log2(_anchorKmerCount);
	cout << endl;
	cout << "Read without anchor: " << (_readWithoutAnchorCount*100) / _readCount << "%" << endl;
	cout << "\tWith N: " << (_noAnchor_with_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	cout << "\tFull N: " << (_noAnchor_full_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	cout << "Anchor hash size: " << _anchorKmerCount*(sizeof(kmer_type)+sizeof(int)) << endl;
	//cout << "Total kmer: " <<  _leon->_total_kmer << endl;
	//cout << "Kmer in bloom:  " << (_leon->_total_kmer_indexed*100) / _leon->_total_kmer << "%" << endl;
	//cout << "Uniq mutated kmer: " << (_leon->_uniq_mutated_kmer*100) / _leon->_total_kmer << "%" << endl;
	//cout << "anchor Hash:   memory_usage: " << _anchorKmers.memory_usage() << "    memory_needed: " << _anchorKmers.size_entry ()*_anchorKmerCount << " (" << (_anchorKmers.size_entry ()*_leon->_anchorKmerCount*100) / _anchorKmers.memory_usage() << "%)" << endl;
	
	cout << endl;
	cout << endl;
	u_int64_t readWithAnchorCount = _readCount - _readWithoutAnchorCount;
	
	//@ anchor kmer
	_readWithAnchorSize += ((double)readWithAnchorCount * log2(_anchorKmerCount)) / 8;
	_readWithAnchorSize += readWithAnchorCount*2;
	
	//anchor dict
	_totalDnaCompressedSize += _anchorKmerCount*sizeof(kmer_type);
	//bloom
    _totalDnaCompressedSize += _bloom->getSize();
	//read without anchor (read_size * 0.375) 3 bit/nt
	_totalDnaCompressedSize += _readWithoutAnchorSize;
	_totalDnaCompressedSize += _readWithAnchorSize;
	
	cout << "Compression rate: " << ((double)_totalDnaCompressedSize / _totalDnaSize) << endl;
	cout << "\tBloom: " << ((_bloom->getSize()*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\tAnchor kmer dictionnary: " << ((_anchorKmerCount*sizeof(kmer_type)*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\tRead with anchor: " << ((_readWithAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\t\t@anchor: " << (((double)readWithAnchorCount * log2(_anchorKmerCount)*100 / 8) / (double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\t\tMutations choices: " << "-" << endl;
	cout << "\t\tOther: " << ((readWithAnchorCount*2*100)/(double)_totalDnaCompressedSize) << "%" << endl;
	cout << "\tRead without anchor: " << ((_readWithoutAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
}

























void Leon::executeDecompression(){
	cout << "Start decompression" << endl;
    string inputFilename = getInput()->getStr(STR_URI_FILE);
    //string inputFilename = prefix + ".txt"; //".leon"
	//_outputFile = System::file().newFile(outputFilename, "wb");
	cout << "\tInput filename: " << inputFilename << endl;
	
	_inputFile = new ifstream(inputFilename.c_str(), ios::in|ios::binary);
	//outputFile.open(outputFilenameDecompressed.c_str(), ios::out);
		
	_rangeDecoder.setInputFile(_inputFile);
	
	//Decode the first byte of the compressed file which is an info byte
	u_int8_t infoByte = _rangeDecoder.nextByte(_generalModel);
	//the first bit holds the file format. 0: fastq, 1: fasta
	_isFasta = ((infoByte & 0x01) == 0x01);
	if(_isFasta)
		cout << "\tInput format: Fasta" << endl;
	else
		cout << "\tInput format: Fastq" << endl;
	
	//Decode the first header
	u_int16_t firstHeaderSize = CompressionUtils::decodeFixedNumeric(_rangeDecoder, _headerNumericModels, 2);
	for(int i=0; i<firstHeaderSize; i++){
		_firstHeader += _rangeDecoder.nextByte(_generalModel);
	}
	cout << _firstHeader << endl;
	u_int64_t blockStartPos = _inputFile->tellg();
	//cout << "Block start pos: " << blockStartPos << endl;
	
	//Go to the end of the file to decode blocks informations, data are read in reversed order (from right to left in the file)
	//The first number is the number of data blocks
	_inputFile->seekg(0, _inputFile->end);
	_rangeDecoder.setInputFile(_inputFile, true);
	u_int64_t blockCount = CompressionUtils::decodeNumeric(_rangeDecoder, _headerNumericSizeModel, _headerNumericModels);
	for(int i=0; i<blockCount; i++){
		_blockSizes.push_back(CompressionUtils::decodeNumeric(_rangeDecoder, _headerNumericSizeModel, _headerNumericModels));
	}
	
	for(int i=0; i<_blockSizes.size(); i++){
		cout << _blockSizes[i] << " ";
	}
	cout << endl;
	
	//Go to the block position
	_inputFile->seekg(blockStartPos, _inputFile->beg);
	
	HeaderDecoder headerDecoder(this, _inputFile, NULL);
	for(int i=0; i<blockCount; i++){
	//for(int i=0; i<1; i++){
		u_int64_t blockSize = _blockSizes[i];
		headerDecoder.setup(blockStartPos, blockSize);
		blockStartPos += blockSize;
	}
	/*
	//Go to the begin of block description and decode blocks informations
	_inputFile->seekg(descriptionStartPos, _inputFile->beg);
	cout << "I am at: " << _inputFile->tellg() << endl;
	_rangeDecoder.setInputFile(_inputFile);
	int headerBlockCount = 2;
	for(int i=0; i<headerBlockCount; i++){
		u_int64_t blockSize = CompressionUtils::decodeNumeric(_rangeDecoder, _headerNumericSizeModel, _headerNumericModels);
		cout << blockSize << endl;
		_blockSizes.push_back(blockSize);
	}
	//CompressionUtils::encodeNumeric(_rangeEncoder, _headerNumericSizeModel, _headerNumericModels, _blockSizes.size());
	*/
}








/*


//main constructor
CompressReads::CompressReads (Bloom<kmer_type>* bloom, Leon * leon, u_int64_t* nb_solids_kmers, int nb_cores, int * nbliving)
: _bloom(bloom), _leon(leon),
_total_nb_solid_kmers_in_reads (nb_solids_kmers), _local_nb_solid_kmers_in_reads(0),
_synchro(System::thread().newSynchronizer()), _nb_living(nbliving)
{
<<<<<<< HEAD
	_thread_id = __sync_fetch_and_add (_nb_living, 1);
=======
    printf("----- create bloom ----\n");
    TIME_INFO (getTimeInfo(), "fill bloom filter");
    
    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    NBITS_PER_KMER = 12;
    u_int64_t solidFileSize = (System::file().getSize(_solidFile) / sizeof (kmer_count));
    
    u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
    
    
    //printf("raw solidFileSize %llu fsize %llu    %lu %lu \n",System::file().getSize(_solidFile), solidFileSize,sizeof (kmer_type),sizeof (kmer_count));
    
    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
                                                                new IteratorFile<kmer_count> (_solidFile),
                                                                solidFileSize,
                                                                "fill bloom filter"
                                                                );
    LOCAL (itKmers);
    
    BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build (itKmers);

	
}

//copy construct
CompressReads::CompressReads(const CompressReads& cr) //called by dispatcher iterate to create N functors
{
	
	//functors share same bloom, bankwriter, leon and synchronizer
	_bloom = cr._bloom;
	_leon = cr._leon;
	_synchro = cr._synchro;
	_nb_living = cr._nb_living;
	_total_nb_solid_kmers_in_reads = cr._total_nb_solid_kmers_in_reads;
	_local_nb_solid_kmers_in_reads =0;
	_thread_id = __sync_fetch_and_add (_nb_living, 1);
}

CompressReads::~CompressReads ()
{
	_thread_id = __sync_fetch_and_add (_total_nb_solid_kmers_in_reads, _local_nb_solid_kmers_in_reads);
	
}

//iterator over sequences, this operator analyzes one sequence
void CompressReads::operator()( Sequence& sequence){
	
	int readlen = sequence.getDataSize();
	char * readseq = sequence.getDataBuffer(); // the nucleotide sequence of the read
	size_t sizeKmer = _leon->_kmerSize;
	
	//cout << sequence.getIndex() << endl;
	
	KmerModel model (sizeKmer,KMER_DIRECT);
	KmerModel::Iterator itKmer (model);
	
	kmer_type current_kmer;
	kmer_type current_kmer_min;
	
	
	itKmer.setData (current_seq.getData());
	

	
	
	int ii=0;

	
	// We iterate the kmers of this sequence
	for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
	{
		
		current_kmer = *itKmer;
		current_kmer_min = min(revcomp(current_kmer, sizeKmer), current_kmer);
		
		
		
		if ( _bloom->contains(current_kmer_min)) //kmer is solid
		{
			
			_local_nb_solid_kmers_in_reads ++;

			
		}
		
		
	} // end of kmers iteration over the read
	

} // end operator () that treats one sequence
*/



