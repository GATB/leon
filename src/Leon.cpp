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

/*
TODO
*  ecrire directement fichier decomp, sans fichier temporatire pour dna, header, qual
* Header coder:
* 	Remplacer la méthode strtoul par une methode string to u_int64_t dans CompressionUtils
* 
* Dna coder:
* 	Sécurité pour ne pas saturer la OAHash _anchorKmers, remplacer map par unordered_map
* 	Le buffer de _anchorRangeEncoder n'est jamais vidé pendant toute l'execution (Solution: utilisé un fichier temporaire) 
*   Centraliser les méthode codeSeedBin, codeSeedNt, nt2bin, bin2nt...
* 
* Optimisation:
* 	methode anchorExist: 2 acces a la map (key exist, puis operator [])
* 	Si les reads on tous la même taille, stocker un bit d'info qui contiendrais l'information et ne pas le stocker à chaque fois	
* 	Compresser le dictionnaire des ancre (codage par dictionnaire, minimizer?)
* 
* Remarque API
* 	Les paramètre du compresseur et decompresseur ne sont pas les même (+ pas besoin de DSK pour le decompresseur)
* 	Comment retrouvé rapidement le chemin du fichier de sortie de dsk (h5) (acutellement code redondant copier depuis dsk)
* 
* Bug:
* 	OAHash: l'algo de decompression plante si _anchorKmers est une OAHash
* 	KmerIt: Si read_size < k, alors le itKmer.setData(seq...) conserve la valeur d'une séquence précédente plutot que d'être remis à 0
* 
* 
 * */
 
 
/*
 * Compressed file :
 * 		
 * 		Headers data blocks
 * 		Bloom
 * 		Anchor Dict
 * 		Reads data blocks
 * 
 * 		Reads blocks position
 * 		Headers blocks position
 * 		First Header
 * 		Info bytes (kmer size, fats aor fastq...)
 * 
 * 
 */
 
 
#include "Leon.hpp"


using namespace std;

//#define SERIAL //this macro is also define in the execute() method
//#define PRINT_DEBUG
//#define PRINT_DEBUG_DECODER



const u_int64_t ANCHOR_KMERS_HASH_SIZE = 500000000;
const char* Leon::STR_COMPRESS = "-c";
const char* Leon::STR_DECOMPRESS = "-d";
const char* Leon::STR_TEST_DECOMPRESSED_FILE = "-test-file";
const char* Leon::STR_DNA_ONLY = "-seq-only";
const char* Leon::STR_NOHEADER = "-noheader";
const char* Leon::STR_NOQUAL = "-noqual";

const int Leon::nt2binTab[128] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0, 0, //69
	0, 3, 0, 0, 0, 0, 0, 0, 4, 0, //79
	0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
	};
const int Leon::bin2ntTab[5] = {'A', 'C', 'T', 'G', 'N'};

/*
const char Leon::_nt2bin[128] = {};
const int Leon::_bin2nt[5] = {};

const char Leon::_nt2bin['A'] = 0;*/
   /*
const char Leon::_nt2bin[] = {};
const char Leon::_nt2bin['A'] = {0};
const char Leon::_nt2bin['C'] = {1};
const char Leon::_nt2bin['T'] = {2};
const char Leon::_nt2bin['G'] = {0};
const char Leon::_nt2bin['N'] = {0};
const char Leon::_bin2nt = {'A', 'C', 'T', 'G', 'N'};
*/

Leon::Leon ( bool compress, bool decompress) :
Tool("leon"),
_generalModel(256),// _anchorKmers(ANCHOR_KMERS_HASH_SIZE),
_anchorDictModel(5),_nb_thread_living(0), _blockwriter(0), //5value: A, C, G, T, N
_readCount (0), _totalDnaSize(0), _compressedSize(0),_MCtotal(0),_MCnoAternative(0),
_MCuniqSolid(0),_MCuniqNoSolid(0),_MCmultipleSolid(0),_readWithoutAnchorCount(0),
_anchorDictSize(0), _anchorAdressSize(0), _anchorPosSize(0), _readSizeSize(0), _bifurcationSize(0), _noAnchorSize(0),
_progress_decode(0), _inputBank(0),_total_nb_quals_smoothed(0),_lossless(false),_input_qualSize(0),_compressed_qualSize(0), _qualwriter(NULL)

{
_isFasta = true;
	std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

    //_kmerSize(27)
	_compress = compress;
	_decompress = decompress;

	/** We don't want default options of Tool (or we want to put them at a specific location). */
	setParser (new OptionsParser ("leon"));

    getParser()->push_back (new OptionOneParam (STR_URI_FILE, "input file (e.g. FASTA/FASTQ for compress or .leon file for decompress)",   true));
    getParser()->push_back (new OptionNoParam  ("-c", "compression",   false));
    getParser()->push_back (new OptionNoParam  ("-d", "decompression", false));
    getParser()->push_back (new OptionOneParam (STR_NB_CORES, "number of cores (default is the available number of cores)", false, "0"));
    getParser()->push_back (new OptionOneParam (STR_VERBOSE,  "verbosity level",                                            false, "1", false));

	
	getParser()->push_back (new OptionNoParam  ("-lossless", "switch to lossless compression for qualities (default is lossy. lossy has much higher compression rate, and the loss is in fact a gain. lossy is better!)",   false));

	
    IOptionsParser* compressionParser = new OptionsParser ("compression");

    /** We add the sorting count options and hide all of them by default and display one some of them. */
    compressionParser->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1, false);

    if (IOptionsParser* input = compressionParser->getParser (STR_URI_INPUT))  {  input->setName (STR_URI_FILE);  }
    if (IOptionsParser* input = compressionParser->getParser (STR_KMER_SIZE))  {  input->setVisible (true);  }

    compressionParser->push_back (new OptionOneParam(STR_KMER_ABUNDANCE, "abundance threshold for solid kmers (default inferred)", false));
	
	compressionParser->push_back (new OptionNoParam (Leon::STR_DNA_ONLY, "store dna seq only, header and quals are discarded, will decompres to fasta (same as -noheader -noqual)", false));

	compressionParser->push_back (new OptionNoParam (Leon::STR_NOHEADER, "discard header", false));
	compressionParser->push_back (new OptionNoParam (Leon::STR_NOQUAL, "discard quality scores", false));

    IOptionsParser* decompressionParser = new OptionsParser ("decompression");
    decompressionParser->push_back (new OptionNoParam (Leon::STR_TEST_DECOMPRESSED_FILE, "check if decompressed file is the same as original file (both files must be in the same folder)", false));
	


	
    getParser()->push_back (compressionParser);
    getParser()->push_back (decompressionParser, 0, false);

	pthread_mutex_init(&findAndInsert_mutex, NULL);
	pthread_mutex_init(&writeblock_mutex, NULL);
}

Leon::~Leon ()
{
	setBlockWriter(0);
	setInputBank (0);
	
	if(_qualwriter != NULL)
		delete(_qualwriter);
	
	if (_progress_decode)  { delete _progress_decode; }
}

void Leon::execute()
{
	_time = clock(); //Used to calculate time taken by decompression
	
	gettimeofday(&_tim, NULL);
	 _wdebut_leon = _tim.tv_sec +(_tim.tv_usec/1000000.0);
	
	
	
	if(getParser()->saw ("-lossless"))
		_lossless = true;
		
	
    //bool compress = false;
    //bool decompress = false;
    //if(getParser()->saw (Leon::STR_COMPRESS)) compress = true;
    //if(getParser()->saw (Leon::STR_DECOMPRESS)) decompress = true;
	if((_compress && _decompress) || (!_compress && !_decompress)){
		cout << "Choose one option among -c (compress) or -d (decompress)" << endl;
		return;
	}

	//getParser()->displayWarnings(); //pb si ici, affiche warnings apres exec dsk ,et prob option -c -d qui sont pas dans le parser car 'globales'
	
    u_int64_t total_nb_solid_kmers_in_reads = 0;
    int nb_threads_living;
    _nb_cores = getInput()->getInt(STR_NB_CORES);
    
    
	//setup global
	for(int i=0; i<CompressionUtils::NB_MODELS_PER_NUMERIC; i++){
		_numericModel.push_back(Order0Model(256));
	}
	
	if(_compress){
		//#define SERIAL
		executeCompression();
	}
	else{
		executeDecompression();
	}
	

        
    
    //outputFile->flush();

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    //getInfo()->add (1, "result");
    //getInfo()->add (2, "nb solid kmers in reads", "%ld", total_nb_solid_kmers_in_reads);
    
    if(_decompress){
		delete _inputFile;
		delete _descInputFile;
	}
    
}

void Leon::createBloom (){
    TIME_INFO (getTimeInfo(), "fill bloom filter");
    
	
	
	//u_int64_t solidFileSize
	
	_auto_cutoff = 0 ;
	u_int64_t nbs = 0 ;
	u_int64_t nb_kmers_infile;
	
	
	//cout << _dskOutputFilename << endl;
	Storage* storage = StorageFactory(STORAGE_HDF5).load (_dskOutputFilename);
	LOCAL (storage);
	
	Partition<kmer_count> & solidCollection = storage->root().getGroup("dsk").getPartition<kmer_count> ("solid");
	
	/** We get the number of solid kmers. */
    u_int64_t solidFileSize = solidCollection.getNbItems();
	
	nb_kmers_infile = solidCollection.getNbItems();
	//(System::file().getSize(_dskOutputFilename) / sizeof (kmer_count)); //approx total number of kmer

	
	
	if( ! getParser()->saw(STR_KMER_ABUNDANCE)){
		
		//retrieve cutoff

		
		Collection<NativeInt64>& cutoff  = storage->getGroup("dsk").getCollection<NativeInt64> ("cutoff");
		Iterator<NativeInt64>* iter = cutoff.iterator();
		LOCAL (iter);
		for (iter->first(); !iter->isDone(); iter->next())  {
			_auto_cutoff = iter->item().toInt();
		}
		//////
		
		//retrieve nb solids
		
		
		Collection<NativeInt64>& storagesolid  = storage->getGroup("dsk").getCollection<NativeInt64> ("nbsolidsforcutoff");
		Iterator<NativeInt64>* iter2 = storagesolid.iterator();
		LOCAL (iter2);
		for (iter2->first(); !iter2->isDone(); iter2->next())  {
			nbs = iter2->item().toInt();
		}
		//////
		
	}
	else
	{
		_auto_cutoff =0;
		nbs  = nb_kmers_infile;

	}
	
	
	


	
	
    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    NBITS_PER_KMER = 12;
    
	
    u_int64_t estimatedBloomSize = (u_int64_t) ((double)nbs * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
    
    
    //printf("raw solidFileSize %llu fsize %llu    %lu %lu \n",System::file().getSize(_solidFile), solidFileSize,sizeof (kmer_type),sizeof (kmer_count));
    
    /** We create the kmers iterator from the solid file. */
//    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
//                                                                new IteratorFile<kmer_count>(_dskOutputFilename),
//                                                                nb_kmers_infile,
//                                                                "fill bloom filter"
//                                                                );
	
	/** We create the kmers iterator from the solid file. */
    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
																solidCollection.iterator(),
																nb_kmers_infile,
																"fill bloom filter"
																);
    LOCAL (itKmers);

	
	

	
    /** We instantiate the bloom object. */
    //BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    //cout << "ESTIMATED:" << estimatedBloomSize << endl;
    //_bloomSize = estimatedBloomSize;
	
	if(_auto_cutoff)
		printf("\tauto cutoff: %i   (total solids %lli)\n",_auto_cutoff,nbs);
	
	//modif ici pour virer les kmers < auto cutoff
    BloomBuilder<> builder (estimatedBloomSize, 7,_kmerSize,tools::misc::BLOOM_NEIGHBOR,getInput()->getInt(STR_NB_CORES),_auto_cutoff);
    _bloom = builder.build (itKmers); // BLOOM_NEIGHBOR // BLOOM_CACHE
	
	//printf ("bloom size %lli  bits per k %f  nbkemrs infile %lli \n",estimatedBloomSize,NBITS_PER_KMER,nbs);
    //BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    //Bloom<kmer_type>* bloom = builder.build (itKmers);
    
	/*
	vector<u_int64_t> tab;
	
	for(int i=0; i<256; i++){
		tab.push_back(0);
	}
	
	for(int i=0; i<_bloom->getSize(); i++){
		tab[_bloom->getArray()[i]] += 1;
	}
	
	cout << tab.size() << " " << _bloom->getSize() << endl;
	for(int i=0; i<tab.size(); i++){
		cout << i << " " << tab[i] << endl;
	}*/
	
	//_rangeEncoder.clear();
	//_generalModel.clear();
	//for(int i=0; i<_bloom->getSize(); i++){
	//	_rangeEncoder.encode(_generalModel, _bloom->getArray()[i]);
	//}

	/*
	z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, Z_BEST_COMPRESSION) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)_bloom->getArray();
    zs.avail_in = _bloom->getSize();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);
    
	cout << "Bloom size: " << _bloom->getSize() << endl;
	cout << "Compressed Bloom size: " << outstring.size() << endl; */
	

}

void Leon::createKmerAbundanceHash(){
	#ifdef PRINT_DEBUG
		cout << "\tFilling kmer abundance hash" << endl;
	#endif

	int ithresholds[4] = {200,50,20,10};
	vector<int> thresholds(ithresholds, ithresholds + sizeof(ithresholds) / sizeof(int));	
	u_int64_t size = 0;
	u_int64_t maxSize = 500000000;
	u_int64_t absoluteMaxSize = (maxSize*3)/4;
    _kmerAbundance = new OAHash<kmer_type>(maxSize);



	
	
    KmerModel model(_kmerSize);
    
    IteratorFile<kmer_count> it(_dskOutputFilename);
    
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
			if(size > absoluteMaxSize) break;
		}
		
		if(size > absoluteMaxSize) break;
	}

	#ifdef PRINT_DEBUG
		cout << "\t\tNeeded memory: " << size << endl;
		cout << "\t\tAllocated memory: " << _kmerAbundance->memory_usage() << endl;
	#endif
}














void Leon::executeCompression(){
	#ifdef PRINT_DEBUG
		cout << "Start compression" << endl;
	#endif
	
    _kmerSize      = getInput()->getInt (STR_KMER_SIZE);
	_nks           = getInput()->getInt (STR_KMER_ABUNDANCE_MIN);
    _inputFilename = getInput()->getStr (STR_URI_FILE);
    
	#ifdef PRINT_DEBUG
		cout << "\tInput filename: " << _inputFilename << endl;
	#endif
	
	u_int8_t infoByte = 0;
	
	
	/** We look for the beginnin of the suffix. */
	int lastindex = _inputFilename.find_last_of (".");
	
	/** We build the result. */
	string extension = _inputFilename.substr(lastindex+1);
	
	_noHeader =false;

	
	if(getParser()->saw (Leon::STR_NOHEADER))
	{
		_noHeader = true;
		infoByte |= 0x02; //no header
	}
	
	if(getParser()->saw (Leon::STR_NOQUAL))
	{
		_isFasta = true;
		infoByte |= 0x01; //fasta mode == no quals
	}
	
	
	if(getParser()->saw (Leon::STR_DNA_ONLY))
	{
		_noHeader = true;
		_isFasta = true;

		infoByte |= 0x02; //no header
		infoByte |= 0x01; //fasta mode == no quals

	}


    //_inputBank = Bank::singleton().createBank(_inputFilename);
	setInputBank (Bank::open(_inputFilename));
	
	//cout << Bank::getType(_inputFilename) << endl;
	
	//guess filename extension
	if( extension.compare("fa")==0 || extension.compare("fasta")==0 )
	//if(_inputFilename.find(".fa") !=  string::npos || _inputFilename.find(".fasta") !=  string::npos)
	{
		//#ifdef PRINT_DEBUG
			cout << "\tInput format: Fasta" << endl;
		//#endif
		infoByte |= 0x01;
		_isFasta = true;
	}
	else //if(_inputFilename.find(".fq") !=  string::npos || _inputFilename.find(".fastq") !=  string::npos)
		if( extension.compare("fq")==0 || extension.compare("fastq")==0 )
	{
			cout << "\tInput format: Fastq";
		
		if(! getParser()->saw (Leon::STR_DNA_ONLY) && ! getParser()->saw (Leon::STR_NOQUAL))
		{
		
		if (_lossless)
			cout << ", compressing qualities in lossless mode" << endl;
			else
		cout << ", compressing qualities in lossy mode (use -lossless for lossless compression)"<< endl;
			
			_isFasta = false;

		}
		
		
	}
	else{
		cout << "\tUnknown input extension. Input extension must be one among fasta (.fa, .fasta) or fastq (.fq, .fastq)" << endl;
		return;
	}
	
	_rangeEncoder.encode(_generalModel, infoByte);
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _kmerSize);
	

    //Redundant from dsk solid file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    _dskOutputFilename = getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) + ".h5"  :
        System::file().getBaseName (_inputFilename) + ".h5"; //_inputFilename instead of prefix GR

#if 1

    /*************************************************/
    // Sorting count part
    /*************************************************/

    {
        Storage* product = StorageFactory(STORAGE_HDF5).create (_dskOutputFilename, true, false);
        LOCAL (product);

        /** We create a DSK instance and execute it. */
        SortingCountAlgorithm<> sortingCount (
            product,
            _inputBank,
            _kmerSize,
            make_pair(_nks,~0),
            getInput()->getInt(STR_MAX_MEMORY),
            getInput()->getInt(STR_MAX_DISK),
            getInput()->getInt(STR_NB_CORES),
            gatb::core::tools::misc::KMER_SOLIDITY_DEFAULT
        );

        sortingCount.getInput()->add (0, STR_VERBOSE, getInput()->getStr(STR_VERBOSE));

        sortingCount.execute();
    }

#endif

    /*************************************************/
    // We create the modified file
    /*************************************************/
    
    string dir = System::file().getDirectory(_inputFilename);
    string prefix = System::file().getBaseName(_inputFilename);
    //_outputFilename = dir + "/" + System::file().getBaseName(prefix) + ".leon";
    _outputFilename = _inputFilename + ".leon";
	_outputFile = System::file().newFile(_outputFilename, "wb");
	
	if(! _isFasta)
	{
		_FileQualname = dir + "/" + System::file().getBaseName(prefix) + ".qual";
		_FileQual = System::file().newFile(_FileQualname, "wb");
		_qualwriter = new OrderedBlocks( _FileQual , _nb_cores ) ;
	}

	setBlockWriter (new OrderedBlocks(_outputFile, _nb_cores ));

	

	
#ifdef PRINT_DEBUG
	cout << "\tOutput filename: " << _outputFilename << endl;
	cout << "prefix " << prefix << endl;
	cout << "dir " << dir << endl;
	cout << "dskout  " << _dskOutputFilename << endl;

#endif

	
	
    //Compression
	if(! _noHeader)
	{
		startHeaderCompression();
	}
	
	

	//_blockWrit
	//setBlockWriter(0);
	//setBlockWriter (new OrderedBlocks(_outputFile, _nb_cores ));

	startDnaCompression();
	

	
	
	endCompression();
}


void Leon::endQualCompression(){
	
	
	//append blocks info at the end
	_FileQual->fwrite(& _qualBlockSizes[0],  sizeof(u_int64_t) , _qualBlockSizes.size());

	u_int64_t val =  _qualBlockSizes.size();
	
	_FileQual->fwrite(& val,  sizeof(u_int64_t) , 1);

	_qualCompRate = ((double)_compressed_qualSize / _input_qualSize);

	
}


void Leon::writeBlockLena(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID){

	
	//printf("write block lena block %i  insize %i \n ",blockID,size);
	//zlib compression for the quals
	
	z_stream zs;
	memset(&zs, 0, sizeof(zs));
	
	//deflateinit2 to be able to gunzip it fro mterminal
	
	//if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
		//			(15+16), 8, Z_DEFAULT_STRATEGY) != Z_OK)
	
	
			if (deflateInit(&zs, Z_BEST_COMPRESSION) != Z_OK)
		throw Exception ("deflateInit failed while compressing.");
	
	zs.next_in = (Bytef*) data ;
	zs.avail_in = size ;           // set the z_stream's input
	
	int ret;
	char outbuffer[32768];
	std::string outstring;
	
	// retrieve the compressed bytes blockwise
	do {
		zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
		zs.avail_out = sizeof(outbuffer);
		
		ret = deflate(&zs, Z_FINISH);
		
		if (outstring.size() < zs.total_out) {
			// append the block to the output string
			outstring.append(outbuffer,
							 zs.total_out - outstring.size());
		}
	} while (ret == Z_OK);
	
	deflateEnd(&zs);
	
	/////////////////
	
	//_qualwriter->insert((u_int8_t*) data, size ,blockID);

	
	_qualwriter->insert((u_int8_t*) outstring.data(), outstring.size() ,blockID);
	_qualwriter->incDone(1);
	
	
	
	
	pthread_mutex_lock(&writeblock_mutex);
	
	_input_qualSize += size;
	_compressed_qualSize +=  outstring.size();
	if ((2*(blockID+1)) > _qualBlockSizes.size() )
	{
		_qualBlockSizes.resize(2*(blockID+1));
	}
	
	_qualBlockSizes[2*blockID] =  outstring.size();
	_qualBlockSizes[2*blockID+1] = encodedSequenceCount;
	
	
	pthread_mutex_unlock(&writeblock_mutex);

	
	
}

void Leon::writeBlock(u_int8_t* data, u_int64_t size, int encodedSequenceCount,u_int64_t blockID){
	if(size <= 0) return;
	
	
	//cout << "\t-----------------------" << endl;
	//cout << "\tWrite block " << _blockSizes.size() << endl;
	//cout << "\tSequence " << encoder->_lastSequenceIndex-READ_PER_BLOCK << " - " << encoder->_lastSequenceIndex << endl;
	//cout << "Thread id: " << thread_id << endl;
	//cout << "\tEncoded size (byte): " << size << endl;
	
	

	
	
	
	
#ifdef SERIAL
	_outputFile->fwrite(data, size, 1);
	
#else
	_blockwriter->insert(data,size,blockID);
	_blockwriter->incDone(1);

#endif
	
	
	pthread_mutex_lock(&writeblock_mutex);

	_compressedSize += size;

	//int thread_id = encoder->getId();

	if ((2*(blockID+1)) > _blockSizes.size() )
	{
		_blockSizes.resize(2*(blockID+1));
	}
	
	_blockSizes[2*blockID] = size ;
	_blockSizes[2*blockID+1] = encodedSequenceCount;
	
	
	pthread_mutex_unlock(&writeblock_mutex);


		
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
	cout << "End compression" << endl;
	
	cout << "\tInput: " << endl;
	cout << "\t\tFilename: " << _inputFilename << endl;
	cout << "\t\tSize: " << inputFileSize <<  "  (" <<  inputFileSize/1024LL/1024LL  << " MB)"<< endl;
	
	u_int64_t outputFileSize = System::file().getSize(_outputFilename);
	if(! _isFasta)
	{
		outputFileSize += System::file().getSize(_FileQualname) ;
	}
	
	cout << "\tOutput: " << endl;
	cout << "\t\tFilename: " << _outputFilename << endl;
	if(! _isFasta)
		cout << "\t\tqualities in : " << _FileQualname << endl;
	
	cout << "\t\tTotal Size: " << outputFileSize <<  "  (" <<  outputFileSize/1024LL/1024LL  << " MB)"<< endl;
	std::cout.precision(4);
	cout << "\tOverall Compression rate: " << (float)((double)outputFileSize / (double)inputFileSize) << " (" << (float)((double)inputFileSize / (double)outputFileSize ) << ")"<< endl;
	if(! _noHeader)
	{
		cout << "\t\tHeader: " << (float)_headerCompRate << " (" << (float) (1/_headerCompRate) << ")" << endl;
	}
	else
	{
		cout << "\t\tHeader completely discarded, in '-seq-only' mode "  << endl;
	}
	cout << "\t\tDna: " << (float)_dnaCompRate << " (" << (float)(1/_dnaCompRate) << ")" << endl ;
	if( _isFasta) cout << endl;
	else
	{
		cout << "\t\tQual: " << (float)_qualCompRate ;
	
		if(_lossless)
			cout << " (in lossless mode) " << endl  << endl;
		else
			cout << " (in lossy mode) " << endl  << endl;
	}

	
	
	
	gettimeofday(&_tim, NULL);
	_wfin_leon  = _tim.tv_sec +(_tim.tv_usec/1000000.0);
	std::cout.precision(2);
	printf("\tTime: %.2fs\n", (  _wfin_leon - _wdebut_leon) );
	printf("\tSpeed: %.2f mo/s\n", (System::file().getSize(_inputFilename)/1000000.0) / (  _wfin_leon - _wdebut_leon) );

	
	//printf("\tTime: %.2fs\n", (double)(clock() - _time)/CLOCKS_PER_SEC);
	//printf("\tSpeed: %.2f mo/s\n", (System::file().getSize(_inputFilename)/1000000.0) / ((double)(clock() - _time)/CLOCKS_PER_SEC));
}
		
		




void Leon::startHeaderCompression(){
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          _inputBank->estimateNbItems(),
                                                          "Compressing headers"
                                                          );
    LOCAL(itSeq);
    
    
	_totalHeaderSize = 0;
	_compressedSize = 0;
	
	#ifdef PRINT_DEBUG
		cout << endl << "Start header compression" << endl;
    #endif
    
    //write first header to file and store it in _firstHeader variable
	//ifstream inputFileTemp(getInput()->getStr(STR_URI_FILE).c_str(), ios::in);
	//getline(inputFileTemp, _firstHeader);   //should be get comment from itseq
	//inputFileTemp.close();
	itSeq->first();
	_firstHeader = itSeq->item().getComment();
	_firstHeader.erase(_firstHeader.begin());
	itSeq->reset();

	#ifdef PRINT_DEBUG
		cout << "\tFirst Header: " << _firstHeader << endl;
		cout << "\tSize: " << _firstHeader.size() << endl;
	#endif
	
	_totalHeaderSize += _firstHeader.size();
	
	//encode the size of the first header on 2 byte and the header itself
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _firstHeader.size());
	for(int i=0; i < _firstHeader.size(); i++){
		_rangeEncoder.encode(_generalModel, _firstHeader[i]);
	}
	//_rangeEncoder.flush();
	//_totalHeaderCompressedSize += _rangeEncoder.getBufferSize();
	//_outputFile->fwrite(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), 1);
	//_rangeEncoder.clear();
	
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
	
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes.size());
	for(int i=0; i<_blockSizes.size(); i++){
		//cout << "block size: " << _blockSizes[i] << endl;
		CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes[i]);
	}
	


		
	_headerCompRate = ((double)_compressedSize / _totalHeaderSize);
	
	//#ifdef PRINT_DEBUG
	cout << "\tEnd header compression" << endl;
	//cout << "\t\tData blocks count: " << _blockSizes.size() << endl;
	//cout << "\tBlock data size: " << _rangeEncoder.getBufferSize() << endl;
	cout << "\t\tHeaders size: " << _totalHeaderSize << endl;
	cout << "\t\tHeaders compressed size: " << _compressedSize << endl;
	std::cout.precision(4);
	cout << "\t\tCompression rate: " << (float)(_headerCompRate) << endl;
		

	//#endif
	//_rangeEncoder.clear();
	_blockSizes.clear();
	//printf("end endHeaderCompression \n");
}



































void Leon::startDnaCompression(){
	#ifdef PRINT_DEBUG
		cout << endl << "Start reads compression" << endl;
    #endif
    
	

	//Create and fill bloom
    createBloom ();
   // LOCAL (_bloom); //now we need it later
    
	int64_t nbestimated = _inputBank->estimateNbItems();
	


//	_anchorKmers = new Hash16<kmer_type, u_int32_t > ( (nbestimated/10) *  sizeof(u_int32_t)  *10LL /1024LL / 1024LL ); // hmm  Hash16 would need a constructor with sizeof main entry //maybe *2 for low coverage dataset
	u_int64_t nbcreated ;
	_anchorKmers = new Hash16<kmer_type, u_int32_t > ( nbestimated/10 , &nbcreated ); //creator with nb entries given
//	printf("asked %lli entries, got %llu \n",nbestimated/10 ,nbcreated);
	
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          _inputBank->iterator(),
                                                          nbestimated,
                                                          "Compressing dna"
                                                          );
    LOCAL(itSeq);
    
	//create a temporary output file to store the anchors dict
	//_dictAnchorFile = System::file().newFile(_outputFilename + ".adtemp", "wb"); 
	_dictAnchorFile = new ofstream((_outputFilename + ".adtemp").c_str(), ios::out|ios::binary);
	
	_lastAnchorValue = 0;
	_anchorAdress = 0;
	_totalDnaSize = 0;
	//_totalDnaCompressedSize = 0;
	_compressedSize = 0;
	
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
	
	if(! _isFasta)
	{
		endQualCompression();
	}
}


void Leon::endDnaCompression(){
	
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes.size());
	for(int i=0; i<_blockSizes.size(); i++){
		//cout << "block size: " << _blockSizes[i] << endl;
		CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _blockSizes[i]);
	}
	_blockSizes.clear();
	
	writeBloom();
	writeAnchorDict();
	
	_dnaCompRate = ((double)_compressedSize / _totalDnaSize);
	
	//#ifdef PRINT_DEBUG
	//	cout << endl;
	//	cout << endl;
	cout << "\tEnd reads compression" << endl;
	cout << "\t\tReads count: " << _readCount << endl;
	cout << "\t\tReads size: " << _totalDnaSize << endl;
	std::cout.precision(4);
	cout << "\t\tCompression rate: " << (float)_dnaCompRate << "  (" << _compressedSize << ")"<< endl;
	std::cout.precision(2);
	cout << "\t\t\tBloom: " << ((_bloom->getSize()*100) / (double)_compressedSize) << "  (" << _bloom->getSize() << ")"<< endl;
	cout << "\t\t\tAnchors dict: " << ((_anchorDictSize*100) / (double)_compressedSize) << "  (" << _anchorDictSize << ")"<< endl;
	u_int64_t readsSize = _anchorAdressSize+_anchorPosSize+_readSizeSize+_bifurcationSize;
	cout << "\t\t\tReads: " << ((readsSize*100) / (double)_compressedSize) << "  (" << readsSize<< ")"<< endl;
	cout << "\t\t\t\tAnchor adress: " << ((_anchorAdressSize*100) / (double)_compressedSize) << "  (" << _anchorAdressSize << ")" << endl;
	cout << "\t\t\t\tAnchor pos: " << ((_anchorPosSize*100) / (double)_compressedSize) << "  (" << _anchorPosSize << ") "<< endl;
	cout << "\t\t\t\tRead size: " << ((_readSizeSize*100) / (double)_compressedSize) << "  (" << _readSizeSize << ")"<< endl;
	cout << "\t\t\t\tBifurcation: " << ((_bifurcationSize*100) / (double)_compressedSize) << "  (" << _bifurcationSize << ")"<< endl;
	cout << "\t\t\tRead without anchor: " << ((_noAnchorSize*100) / (double)_compressedSize) << "  (" << _noAnchorSize << ")"<< endl;
		
		
		
	//cout << "\t\tReads stats" << endl;
	//cout << System::file().getBaseName(getInput()->getStr(STR_URI_FILE)) << endl;
	if(_anchorAdress!=0)
		cout << "\t\tReads per anchor: " << _readCount / _anchorAdress << endl;
	//cout << "\tBit per anchor: " << log2(_anchorKmerCount) << endl;
	//cout << "\tAnchor count: " << _anchorKmerCount << endl;
	//cout << endl;
	cout << "\t\tRead without anchor: " << ((double)_readWithoutAnchorCount*100) / _readCount << endl;
	cout << "\t\tDe Bruijn graph" << endl;
	
	//cout << "\t\t\t\tTotal encoded nt: " << _MCtotal << endl;
	//cout << "\t\t\t\tTotal _MCmultipleNoSolid nt: " << _MCmultipleNoSolid << endl;

	cout << "\t\t\tSimple path: " << ((_MCuniqSolid*100)/ (double)_MCtotal) << endl;
	cout << "\t\t\tBifurcation: " << ((_MCmultipleSolid*100)/(double)_MCtotal) << endl;
	cout << "\t\t\tBreak: " << ((_MCnoAternative*100)/(double)_MCtotal) << endl;
	cout << "\t\t\tError: " << ((_MCuniqNoSolid*100)/(double)_MCtotal) << endl;
	//cout << "\t\t\tOther: " << ((_MCmultipleNoSolid*100)/(double)_MCtotal) << endl;
	
	//cout << "\t\tWith N: " << (_noAnchor_with_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	//cout << "\t\tFull N: " << (_noAnchor_full_N_kmer_count*100) / _readWithoutAnchorCount << "%" << endl;
	//cout << "\tAnchor hash size: " << _anchorKmerCount*(sizeof(kmer_type)+sizeof(int)) << endl;
	//cout << "Total kmer: " <<  _leon->_total_kmer << endl;
	//cout << "Kmer in bloom:  " << (_leon->_total_kmer_indexed*100) / _leon->_total_kmer << "%" << endl;
	//cout << "Uniq mutated kmer: " << (_leon->_uniq_mutated_kmer*100) / _leon->_total_kmer << "%" << endl;
	//cout << "anchor Hash:   memory_usage: " << _anchorKmers.memory_usage() << "    memory_needed: " << _anchorKmers.size_entry ()*_anchorKmerCount << " (" << (_anchorKmers.size_entry ()*_leon->_anchorKmerCount*100) / _anchorKmers.memory_usage() << "%)" << endl;
	//cout << endl;
	//#endif
	
	//u_int64_t readWithAnchorCount = _readCount - _readWithoutAnchorCount;
	
	//@ anchor kmer
	//_readWithAnchorSize += ((double)readWithAnchorCount * log2(_anchorKmerCount)) / 8;
	//_readWithAnchorSize += readWithAnchorCount*2;
	//_readWithAnchorSize += _readWithAnchorMutationChoicesSize;
	
	//anchor dict
	//_totalDnaCompressedSize += _anchorKmerCount*sizeof(kmer_type);
	//bloom
    //_totalDnaCompressedSize += _bloom->getSize();
	//read without anchor (read_size * 0.375) 3 bit/nt
	//_totalDnaCompressedSize += _readWithoutAnchorSize;
	//_totalDnaCompressedSize += _readWithAnchorSize;
	
	//#ifdef PRINT_DEBUG
		//cout << "\tCompression rate: " << ((double)_totalDnaCompressedSize / _totalDnaSize) << "    " << (float)_dnaCompRate << endl;
		//cout << "\t\tAnchor kmer dictionnary: " << ((System::file().getSize(_outputFilename + ".adtemp")*100) / (double)_totalDnaCompressedSize) << "%" << endl;
		//cout << "\t\tRead with anchor: " << ((_readWithAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
		//cout << "\t\t\t@anchor: " << (((double)readWithAnchorCount * log2(_anchorKmerCount)*100 / 8) / (double)_readWithAnchorSize) << "%" << endl;
		//cout << "\t\t\tMutations choices: " << ((_readWithAnchorMutationChoicesSize*100) / (double)_readWithAnchorSize)  << endl;
		//cout << "\t\t\tOther: " << ((readWithAnchorCount*2*100)/(double)_readWithAnchorSize) << "%" << endl;
		//cout << "\t\tRead without anchor: " << ((_readWithoutAnchorSize*100) / (double)_totalDnaCompressedSize) << "%" << endl;
		//cout << endl;
	//#endif
	//_anchorKmers.clear();
	//_anchorKmers->clear();
	
	delete _anchorKmers;
	System::file().remove(_dskOutputFilename);


}

void Leon::writeBloom(){
	//_bloom->save(_outputFilename + ".bloom.temp");
	//cout << _bloom->getBitSize() << endl;
	_compressedSize += _bloom->getSize();
	
	//_outputFile->fwrite(_anchorRangeEncoder.getBuffer(), size, 1);

	//u_int64_t size = _bloom->getSize();
	//CompressionUtils::encodeNumeric(_rangeEncoder, _numericSizeModel, _numericModel, size);
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _bloom->getBitSize());
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _bloom->getNbHash());
	
	//u_int8_t*& bloomData = _bloom->getArray();
	_outputFile->fwrite(_bloom->getArray(), _bloom->getSize(), 1);
	//cout << "Bloom size: " << _bloom->getSize() << endl;
	

	//_outputFile->fwrite(_rangeEncoder.getBuffer(), _rangeEncoder.getBufferSize(), 1);
}

void Leon::writeAnchorDict(){
	_anchorRangeEncoder.flush();
	
	_dictAnchorFile->write( (const char*) _anchorRangeEncoder.getBuffer(), _anchorRangeEncoder.getBufferSize());
	_dictAnchorFile->flush();
	_dictAnchorFile->close();
	_anchorRangeEncoder.clear();
	
	//cout << "lololol: " << _outputFile->tell() << endl;
	
	u_int64_t size = System::file().getSize(_outputFilename + ".adtemp");
	_anchorDictSize = size;
	//u_int64_t size = _anchorRangeEncoder.getBufferSize();
	_compressedSize += size;
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, size);
	
	//Encode anchors count
	CompressionUtils::encodeNumeric(_rangeEncoder, _numericModel, _anchorAdress);
	
	//printf("should encode %u anchors \n",_anchorAdress);
	//cout << "Anchor dict size: " << System::file().getSize(_outputFilename + ".adtemp") << endl;
	//cout << "\t pos: " << _outputFile->tell() << endl;
	//cout << "count: " << _anchorKmerCount << endl;
	
	//_dictAnchorFile->seeko(0, SEEK_SET);
	//_outputFile->fwrite(_dictAnchorFile, size, 1);
	ifstream tempFile((_outputFilename + ".adtemp").c_str(), ios::in|ios::binary);
	
	
	int bufsize = 4096*8;
	char * buffer = new char [bufsize];
	

    while (tempFile.good()) {
		tempFile.read(buffer, bufsize);
        _outputFile->fwrite(buffer, tempFile.gcount(), 1);
    }
    
	tempFile.close();
	remove((_outputFilename + ".adtemp").c_str());
	delete [] buffer;
}

bool Leon::anchorExist(const kmer_type& kmer, u_int32_t* anchorAdress){
	
	if (_anchorKmers->get(kmer,anchorAdress)) //avec Hash16
	{
		return true;
	}
	
//	if(_anchorKmers.find( kmer ) != _anchorKmers.end()){ //avec std map
//		*anchorAdress = _anchorKmers[kmer];
//		return true;
//	}
	
	return false;
	//return _anchorKmers.get(kmer, (int*)anchorAdress); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!changer OAHash in to u_int32_t
	//u_int32_t anchorAdress;
	//bool exist = _anchorKmers.get(kmer, &anchorAdress);
	//return _anchorKmers.has_key(kmer);
}

int Leon::findAndInsertAnchor(const vector<kmer_type>& kmers, u_int32_t* anchorAdress){
	
	pthread_mutex_lock(&findAndInsert_mutex);

		
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
		
		
		/*
		if(_bloom->contains(kmerMin)){
			maxAbundance = 0;
			bestPos = i;
			bestKmer = kmerMin;
			break;
		}
		*/
		
		int abundance;
		if(_kmerAbundance->get(kmerMin, &abundance)){
			if(abundance > maxAbundance){
				maxAbundance = abundance;// + ((kmers.size()-i)*2);
				bestKmer = kmerMin;
				bestPos = i;
				//cout << maxAbundance << endl;
				//cout << bestPos << " " << abundance << " " << kmer.toString(_kmerSize) << " " << revcomp(kmer, _kmerSize).toString(_kmerSize) << endl;
			}
			//cout << abundance << endl;
		}
		else if(maxAbundance == -1 && _bloom->contains(kmerMin)){
			maxAbundance = _nks;
			bestKmer = kmerMin;
			bestPos = i;
			//cout << maxAbundance << endl;
		}
	}
	
	if(maxAbundance == -1)
	{
		pthread_mutex_unlock(&findAndInsert_mutex);
		return -1;
	}

	encodeInsertedAnchor(bestKmer);
	
	_anchorKmers->insert(bestKmer,_anchorAdress); //with Hash16
	//_anchorKmers[bestKmer] = _anchorAdress;
	//_anchorKmers.insert(bestKmer, _anchorAdress);
		
	
	*anchorAdress = _anchorAdress;
	//_anchorKmerCount += 1;
	_anchorAdress += 1;
	
	pthread_mutex_unlock(&findAndInsert_mutex);
	return bestPos;
}

void Leon::encodeInsertedAnchor(const kmer_type& kmer){

	/*
	u_int64_t kmerValue = kmer.getVal();
	u_int64_t deltaValue = 0;

	int deltaType = CompressionUtils::getDeltaValue(kmerValue, _lastAnchorValue, &deltaValue);
	_anchorRangeEncoder.encode(_anchorDictModel, deltaType);
	CompressionUtils::encodeNumeric(_anchorRangeEncoder, _numericModel, deltaValue);
	_lastAnchorValue = kmerValue;

	cout << deltaValue << endl;
	*/

	//static int i = 0;
	
	string kmerStr = kmer.toString(_kmerSize);

	for(int i=0; i<kmerStr.size(); i++){
		_anchorRangeEncoder.encode(_anchorDictModel, Leon::nt2bin(kmerStr[i]));
	}
	//i+= 1;
	//cout << i << endl;
	if(_anchorRangeEncoder.getBufferSize() >= 4096){
		_dictAnchorFile->write((const char*) _anchorRangeEncoder.getBuffer(), _anchorRangeEncoder.getBufferSize());
		_anchorRangeEncoder.clearBuffer();
	}
}







typedef struct
{
	QualDecoder * qual_decoder;
	HeaderDecoder * header_decoder;
	DnaDecoder * dna_decoder;
} thread_arg_decoder;




void * decoder_all_thread(void * args)
{
	
	thread_arg_decoder * targ = (thread_arg_decoder*) args;
	QualDecoder * qual_decoder = targ->qual_decoder;
	DnaDecoder * dna_decoder = targ->dna_decoder;
	HeaderDecoder * header_decoder = targ->header_decoder;
	
	if(qual_decoder!=NULL)
		qual_decoder->execute();
	
	if(header_decoder!=NULL)
		header_decoder->execute();

	dna_decoder->execute();
	
	pthread_exit(0);
}





void * decoder_header_thread(void * args)
{
	HeaderDecoder * header_decoder = (HeaderDecoder*) args;
	//header_decoder->execute();
 	pthread_exit(0);
}

void * decoder_dna_thread(void * args)
{
	DnaDecoder * dna_decoder = (DnaDecoder*) args;
	dna_decoder->execute();
 	pthread_exit(0);
}

void * decoder_qual_thread(void * args)
{
	QualDecoder * qual_decoder = (QualDecoder*) args;
	qual_decoder->execute();
	pthread_exit(0);
}

void Leon::executeDecompression(){

	_filePos = 0;
	
	cout << "Start decompression" << endl;
    _inputFilename = getInput()->getStr(STR_URI_FILE);
    //string inputFilename = prefix + ".txt"; //".leon"
	//_outputFile = System::file().newFile(outputFilename, "wb");
	cout << "\tInput filename: " << _inputFilename << endl;
	
	
	string dir = System::file().getDirectory(_inputFilename);

	_descInputFile = new ifstream(_inputFilename.c_str(), ios::in|ios::binary);
	_inputFile = new ifstream(_inputFilename.c_str(), ios::in|ios::binary);
	

if ( (_inputFile->rdstate() & std::ifstream::failbit ) != 0 )	
{
fprintf(stderr,"cannot open file %s\n",_inputFilename.c_str());
exit( EXIT_FAILURE);
}
	//Go to the end of the file to decode blocks informations, data are read in reversed order (from right to left in the file)
	//The first number is the number of data blocks
	_descInputFile->seekg(0, _descInputFile->end);
	_rangeDecoder.setInputFile(_descInputFile, true);
	
	//_rangeDecoder.setInputFile(_descInputFile);
	
	//Output file
	string prefix = System::file().getBaseName(_inputFilename);;
	while(prefix.find('.') != string::npos){
		int lastindex = prefix.find_last_of(".");
		prefix = prefix.substr(0, lastindex);
	}

    //string prefix = System::file().getBaseName(_inputFilename);
	_outputFilename = dir + "/" + prefix;
	
	//Decode the first byte of the compressed file which is an info byte
	u_int8_t infoByte = _rangeDecoder.nextByte(_generalModel);
	
	//the first bit holds the file format. 0: fastq, 1: fasta
	_isFasta = ((infoByte & 0x01) == 0x01);
	
	

	//Second bit : option no header
	_noHeader = ((infoByte & 0x02) == 0x02);
	
	//printf("info byte %i  _noHeader %i  _isFasta %i \n",infoByte,_noHeader,_isFasta);
	
	if(! _isFasta)
	{
		
		_FileQualname =    dir + "/" +  System::file().getBaseName(_inputFilename) + ".qual";
		_inputFileQual = new ifstream(_FileQualname.c_str(), ios::in|ios::binary);
		cout << "\tQual filename: " << _FileQualname << endl;
	}
	
	
if(_noHeader)
{
	cout << "Headers were not stored, will number reads \n" << endl;
}
	
	if(_isFasta){
		cout << "\tOutput format: Fasta" << endl;
		_outputFilename += "_d.fasta";
	}
	else{
		cout << "\tOutput format: Fastq" << endl;
		_outputFilename += "_d.fastq";
	}
	
	_outputFile = System::file().newFile(_outputFilename, "wb"); 

	//Get kmer size
	_kmerSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	cout << "\tKmer size: " << _kmerSize << endl;
	cout << endl;
	
	startDecompressionAllStreams();
	
	/*
	startHeaderDecompression();
	startDnaDecompression();

	if(! _isFasta){
		startQualDecompression();
	}
	*/
	
	endDecompression();
}


void Leon::startDecompressionAllStreams(){
	
	
	
	_filePosHeader = 0;
	_filePosDna = 0;
	
	if(! _noHeader)
	{
		
	///////// header setup  /////////
	//Decode the first header
	u_int16_t firstHeaderSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(int i=0; i<firstHeaderSize; i++){
		_firstHeader += _rangeDecoder.nextByte(_generalModel);
	}
	
	setupNextComponent(_headerBlockSizes);
	
	}
	
	
	/////// dna setup ////////////
	
	//need to init _filePosDna here
	for(int ii=0; ii<_headerBlockSizes.size(); ii+=2 )
	{
		_filePosDna += _headerBlockSizes[ii];
	}
	
	setupNextComponent(_dnaBlockSizes);
	
	_kmerModel = new KmerModel(_kmerSize);

	
	decodeBloom();
	decodeAnchorDict();
	
	
	/////////// qualities setup //////////
	if(! _isFasta)
	{
	_filePosQual =0;
	
	//read block sizes and _blockCount
	_qualBlockSizes.clear();
	_inputFileQual->seekg(- sizeof(u_int64_t),_inputFileQual->end);
	
	_inputFileQual->read((char *)&_blockCount,sizeof(u_int64_t));
	//cout << "\tBlock count: " << _blockCount/2 << endl;
	
	_qualBlockSizes.resize(_blockCount,0);
	char * databuff = (char * )& _qualBlockSizes[0];
	
	_inputFileQual->seekg(- (sizeof(u_int64_t)*(_blockCount+1)),_inputFileQual->end);
	_inputFileQual->read( databuff ,sizeof(u_int64_t) *  _blockCount);
	
	}
	
	///////////////
	
	switch (getInput()->getInt(STR_VERBOSE))
	{
		case 0: default:    _progress_decode =  new IteratorListener ();break;
		case 1:             _progress_decode = new ProgressSynchro ( new ProgressTimer ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );break;
			
		case 2:             _progress_decode = new ProgressSynchro ( new Progress ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );break;
	}
	
	
	cout << "\tBlock count: " << _blockCount/2 << endl;


	//delete _progress_decode;
	//_progress_decode = new ProgressSynchro ( new ProgressTimer ( _blockCount/2, "Decompressing all streams"), System::thread().newSynchronizer()   );
	_progress_decode->init();
	
	
	
	vector<QualDecoder*> qualdecoders;
	vector<DnaDecoder*> dnadecoders;
	vector<HeaderDecoder*> headerdecoders;
	
	for(int i=0; i<_nb_cores; i++){
		
		if(! _isFasta)
		{
			QualDecoder* qd = new QualDecoder(this, _FileQualname);
			qualdecoders.push_back(qd);
		}
		
		DnaDecoder* dd = new DnaDecoder(this, _inputFilename);
		dnadecoders.push_back(dd);
		
		if(! _noHeader)
		{
		HeaderDecoder* hd = new HeaderDecoder(this, _inputFilename);
		headerdecoders.push_back(hd);
		}
	}
	
	pthread_t * tab_threads = new pthread_t [_nb_cores];
	
	thread_arg_decoder *  targ = new thread_arg_decoder [_nb_cores];

	
	int i = 0;
	int livingThreadCount = 0;
	
	

	
	while(i < _dnaBlockSizes.size()){
		
		for(int j=0; j<_nb_cores; j++){
			

			if(i >= _dnaBlockSizes.size()) break;
			
			
			u_int64_t blockSize;
			int sequenceCount;
			
			livingThreadCount = j+1;
			
			QualDecoder* qdecoder;
			HeaderDecoder* hdecoder;
			DnaDecoder* ddecoder;
			
			//header decoder
			if(! _noHeader)
			{
				blockSize = _headerBlockSizes[i];
				sequenceCount = _headerBlockSizes[i+1];
				hdecoder = headerdecoders[j];
				hdecoder->setup(_filePosHeader, blockSize, sequenceCount);
				_filePosHeader += blockSize;
			}
			else
			{
				hdecoder= NULL;
			}
			
			//dna decoder
			blockSize = _dnaBlockSizes[i];
			sequenceCount = _dnaBlockSizes[i+1];
			ddecoder = dnadecoders[j];
			ddecoder->setup(_filePosDna, blockSize, sequenceCount);
			_filePosDna += blockSize;

			//qual decoder setup
			//here test if in fastq mode, put null pointer otherwise
			if(! _isFasta)
			{
				blockSize = _qualBlockSizes[i];
				sequenceCount = _qualBlockSizes[i+1];
				qdecoder = qualdecoders[j];
				qdecoder->setup(_filePosQual, blockSize, sequenceCount);
				_filePosQual += blockSize;
			}
			else
			{
				qdecoder= NULL;
			}
			
			
			targ[j].qual_decoder = qdecoder;
			targ[j].dna_decoder = ddecoder;
			targ[j].header_decoder = hdecoder;
			
			pthread_create(&tab_threads[j], NULL, decoder_all_thread, targ + j);
			

			i += 2;
			
			this->_progress_decode->inc(1);

		}

		
		
		
		for(int j=0; j < livingThreadCount; j++){
			
			pthread_join(tab_threads[j], NULL);
			
			HeaderDecoder* hdecoder = NULL;
			QualDecoder* qdecoder = NULL;
			DnaDecoder* ddecoder = dnadecoders[j];
			

			std::istringstream  * stream_qual;
			std::istringstream  * stream_header;

			if(! _isFasta)
			{
				qdecoder = qualdecoders[j];
				stream_qual = new std::istringstream (qdecoder->_buffer);
				qdecoder->_buffer.clear();

			}
			
			if(! _noHeader)
			{
				hdecoder = headerdecoders[j];
				stream_header = new std::istringstream (hdecoder->_buffer);
				hdecoder->_buffer.clear();

			}
			
			std::istringstream stream_dna (ddecoder->_buffer);
			
			ddecoder->_buffer.clear();

			
			std::string line;
			std::string output_buff;

			
			
			output_buff.reserve(READ_PER_BLOCK * 500);
			
			bool reading = true;
			
			
			
			u_int64_t readid=0;
			while(reading){
				
				stringstream sint;
				sint << readid;
				
				if( ! _noHeader)
				{
					if(getline(*stream_header, line)){
						if(_isFasta)
							output_buff += ">";
						else
							output_buff += "@";
						
						output_buff +=  line + '\n';
					}
					else
						reading = false;
				}
				else
				{
					if(_isFasta)
						output_buff += "> " + sint.str() + '\n';
					else
						output_buff += "@ " + sint.str() + '\n';
					
					readid++;
				}
				 
				
				
				if(getline(stream_dna, line)){
					output_buff +=  line + '\n';
				}
				else
					reading = false;
				
				
				if( ! _isFasta)
				{
					if(getline(*stream_qual, line)){
						output_buff += "+\n";
						output_buff +=  line + '\n';
					}
					else
						reading = false;
				}
				 
			}
			
			 
			_outputFile->fwrite(output_buff.c_str(), output_buff.size(), 1);


			
		}
		
		livingThreadCount = 0;
	}
	
	
	
	_outputFile->flush();

	for(int i=0; i<dnadecoders.size(); i++){
		delete dnadecoders[i];
	}
	dnadecoders.clear();
	
	for(int i=0; i<headerdecoders.size(); i++){
		delete headerdecoders[i];
	}
	headerdecoders.clear();
	
	for(int i=0; i<qualdecoders.size(); i++){
		delete qualdecoders[i];
	}
	qualdecoders.clear();
	
	
	
	delete [] tab_threads;
	delete [] targ;
	cout << endl;
	
	
	_progress_decode->finish();
	
	delete _kmerModel;

}



void Leon::setupNextComponent( 		vector<u_int64_t>   & blockSizes    ){
	//Go to the data block position (position 0 for headers, position |headers data| for reads)
	_inputFile->seekg(_filePos, _inputFile->beg);
	
	blockSizes.clear();
	//u_int64_t size = 0;
	
	_blockCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	for(int i=0; i<_blockCount; i++){
		u_int64_t blockSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
		blockSizes.push_back(blockSize);
		//size += blockSize;
	}
	
	

	
	//cout << "\tBlock count: " << _blockCount/2 << endl;
	/*
	for(int i=0; i<_blockSizes.size(); i++){
		cout << _blockSizes[i] << " ";
	}
	cout << endl;*/
	
}



void Leon::decodeBloom(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\tDecode bloom filter" << endl;
	#endif
	
	//pos = tous les block sizes des header
	u_int64_t total_header_block_size = 0 ;
	
	for(int ii=0; ii<_headerBlockSizes.size(); ii+=2 )
	{
		total_header_block_size  += _headerBlockSizes[ii];
	}
	
	u_int64_t bloomPos =  total_header_block_size ; // was _inputFile->tellg();

	//printf("\npos at beginning of decodebloom %llu  total_header_block_size %llu nb vecblock size %i\n",bloomPos,total_header_block_size,_headerBlockSizes.size());
	for(int i=0; i<_dnaBlockSizes.size(); i++){
		bloomPos += _dnaBlockSizes[i];
		i += 1;
	}
	//printf(".. pos after setup %llu \n",bloomPos);

	//cout << "Anchor dict pos: " << dictPos << endl;
	_inputFile->seekg(bloomPos, _inputFile->beg);
	
	//u_int64_t bloomSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericSizeModel, _numericModel);
	
	//char* buffer = new char[bloomSize];
    //_inputFile->read(buffer, bloomSize);
    
    //cout << bloomSize << endl;
    
	u_int64_t bloomBitSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	u_int64_t bloomHashCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	//u_int64_t tai = (bloomSize) * 8LL;
	
	
	//cout << tai << endl;
	
	//_bloom->contains("lala");
	//nchar  = (1+tai/8LL);
	//_bloom = new BloomCacheCoherent<kmer_type>(bloomBitSize, bloomHashCount);
	_bloom   =  new BloomNeighborCoherent<kmer_type> (bloomBitSize,_kmerSize,bloomHashCount);
	

	_inputFile->read((char*)_bloom->getArray(), _bloom->getSize());
	//fread(_bloom->getArray(), sizeof(unsigned char), result->getSize(), file);
	
	//_dskOutputFilename = "/local/gbenoit/leontest/SRR747737.h5";
	//createBloom();
	//cout << _bloom->getSize() << endl;
	//memcpy(_bloom, buffer, _bloom->getSize());
	//u_int8_t*& bloomData = _bloom->getArray();
	//_outputFile->fwrite(bloomData, _bloom->getSize(), 1);
	
	//_inputFile->seekg(_bloom->getSize(), _inputFile->cur);
	#ifdef PRINT_DEBUG_DECODER
		cout << "Bloom size: " << _bloom->getSize() << endl;
		cout << "Anchor dict pos: " << _inputFile->tellg() << endl;
	#endif
	

}

void Leon::decodeAnchorDict(){
	#ifdef PRINT_DEBUG_DECODER
		cout << "\tDecode anchor dict" << endl;
	#endif
	
	//_anchorDictFilename = _outputFilename + ".temp.dict";
	//ofstream anchorDictFile(_anchorDictFilename.c_str(), ios::out);
	
	u_int64_t anchorDictSize = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);
	//cout << anchorDictSize << endl;
	
	u_int64_t anchorCount = CompressionUtils::decodeNumeric(_rangeDecoder, _numericModel);

	_anchorRangeDecoder.setInputFile(_inputFile);
	string anchorKmer = "";
	//int lala =0;
	
	//_vecAnchorKmers.push_back(kmer_type()); //Insert a random kmer at first pos because adresses start at 1 (OAHash)
	u_int64_t dictPos = _inputFile->tellg();
	
	//KmerModel model(_kmerSize, KMER_DIRECT);
	//int i=0;
	//cout << _inputFile->tellg()  << " " << dictPos+anchorDictSize << endl;
	//while(_inputFile->tellg() < dictPos+anchorDictSize){
	u_int64_t currentAnchorCount = 0;
	
	while(currentAnchorCount < anchorCount){


		u_int8_t c = _anchorRangeDecoder.nextByte(_anchorDictModel);
		anchorKmer += Leon::bin2nt(c); //convert to char
		if(anchorKmer.size() == _kmerSize){
			

			//cout << anchorKmer << endl;
			//if(i<=10) cout << anchorKmer << endl;
			//cout << "1: " << anchorKmer << endl;
			
			kmer_type kmer = _kmerModel->codeSeed(anchorKmer.c_str(), Data::ASCII).value() ; //then convert to bin
			

			//could be optimized if needed
			//cout << "2: " << model.toString(kmer) << endl;
			//lala += 1;
			_vecAnchorKmers.push_back(kmer);
			

			//anchorDictFile << anchorKmer;// + '\n'; //A remettre pour le mode anchor dict sur disk
			//cout << anchorKmer << endl;
			anchorKmer.clear();
			//i++;
			//if(i > 50) break;
			currentAnchorCount += 1;

		}
	}
	
	#ifdef PRINT_DEBUG_DECODER
		cout << "\t\tAnchor count: " << _vecAnchorKmers.size() << endl;
	#endif
	
	//anchorDictFile.flush();
	//anchorDictFile.close();
}



kmer_type Leon::getAnchor(ifstream* anchorDictFile, u_int32_t adress){
	
	return _vecAnchorKmers[adress];  //icii
	
	//anchorDictFile->seekg(_kmerSize*adress);
	
	//char buffer[_kmerSize];
	
	//anchorDictFile->read(buffer, _kmerSize);
	//kmer_type kmer = model.codeSeed(anchorKmer.c_str(), Data::ASCII);
	//return _vecAnchorKmers[adress];
	//return _kmerModel->codeSeed(buffer, Data::ASCII);
}

void Leon::endDecompression(){
	
	cout << "\tOutput filename: " << _outputFile->getPath() << endl;

	gettimeofday(&_tim, NULL);
	_wfin_leon  = _tim.tv_sec +(_tim.tv_usec/1000000.0);
	
	printf("\tTime: %.2fs\n", (  _wfin_leon - _wdebut_leon) );
	printf("\tSpeed: %.2f mo/s\n", (System::file().getSize(_outputFilename)/1000000.0) / (  _wfin_leon - _wdebut_leon) );
	
	
	//Test decompressed file against original reads file (decompressed and original read file must be in the same dir)
	if(getParser()->saw (Leon::STR_TEST_DECOMPRESSED_FILE)){
		
		cout << endl << "\tChecking decompressed file" << endl;
		
		string dir = System::file().getDirectory(_inputFilename);

		string prefix = System::file().getBaseName(_inputFilename);;
		while(prefix.find('.') != string::npos){
			int lastindex = prefix.find_last_of(".");
			prefix = prefix.substr(0, lastindex);
		}
		//string prefix = System::file().getBaseName(_inputFilename);
		
		string originalFilename;
		IBank* originalBank;
		IBank* newBank;
		Iterator<Sequence>* originalBankIt;
		Iterator<Sequence>* newBankIt;
		
		if(_isFasta)
			originalFilename = dir + "/" + prefix + ".fasta";
		else
			originalFilename = dir + "/" + prefix + ".fastq";
		
		
		cout << "\t\tOriginal file: " << originalFilename << endl;
		cout << "\t\tNew file: " << _outputFile->getPath() << endl;
		
		originalBank = Bank::open(originalFilename);
		originalBankIt = originalBank->iterator();
		originalBankIt->first();
		newBank = Bank::open(_outputFile->getPath());
		newBankIt = newBank->iterator();
		newBankIt->first();
		
		//int i=0;
		
		while(true){
			if(newBankIt->isDone()){
				if(originalBankIt->isDone())
					cout << "\tOK" << endl;
				else
					cout << "\tDecompressed file end but not the original file" << endl;
				break;
			}
			if(originalBankIt->isDone()){
				if(newBankIt->isDone())
					cout << "\tOK" << endl;
				else
					cout << "\tOriginal file end but not the decomrpessed file" << endl;
				break;
			}
			
			string originalHeader = (*originalBankIt)->getComment();
			string originalDna = (string((*originalBankIt)->getDataBuffer())).substr(0, (*originalBankIt)->getDataSize());
			
			
			string newHeader = (*newBankIt)->getComment();
			string newDna = (string((*newBankIt)->getDataBuffer())).substr(0, (*newBankIt)->getDataSize());
			
			if(originalHeader != newHeader){
				cout << "\tSeq " << (*newBankIt)->getIndex() << "    Header different" << endl;
				cout << "\t\t" << originalHeader << endl;
				cout << "\t\t" << newHeader << endl;
				break;
			}
			
			if(originalDna != newDna){
				cout << "\tSeq " << (*newBankIt)->getIndex() << "    Dna different" << endl;
				cout << "\t\t" << originalDna << endl;
				cout << "\t\t" << newDna << endl;
				break;
			}
			
			originalBankIt->next();
			newBankIt->next();
			
			//i ++;
			//cout << i << endl;
			//if(i > 20) return;
		}
	}
}






