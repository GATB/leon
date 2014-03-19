//
//  Leon.cpp
//  leon
//
//  Created by Guillaume Rizk on 17/01/2014.
//  Copyright (c) 2014 rizk. All rights reserved.
//

#include "Leon.hpp"
#include <DSK.hpp>

using namespace std;




//const char* Leon::STR_SMOOTHING_THRESHOLD         = "-smooth";




class CompressReads
{
public:
    void operator() ( Sequence& current_seq) //iterator over sequences, this operator analyzes one sequence
    {
        
        int nb_quals_smoothed = 0;

        int readlen = current_seq.getDataSize();


        
        char * readseq = current_seq.getDataBuffer(); // the nucleotide sequence of the read
        

        
        size_t sizeKmer = _leon->_kmerSize;
        
        
//        printf("%.*s\n",readlen,readseq);
        
        
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
    
    
    
    
    
    Bloom<kmer_type> * _bloom; // the bloom containing the solid kmers
    Leon *         _leon; // the parent Leon object
    
    int *  _nb_living;
    int _thread_id;
    
	
    
    ISynchronizer* _synchro;
    ISynchronizer* getSynchro ()  { return _synchro; }
    
    u_int64_t *  _total_nb_solid_kmers_in_reads;
    u_int64_t   _local_nb_solid_kmers_in_reads;
    

    
    
    
    //default constructor
    CompressReads ()
    : _bloom(NULL),  _leon(NULL),
    _total_nb_solid_kmers_in_reads (NULL), _local_nb_solid_kmers_in_reads(0)
    {
        
    }
    
    
    //main constructor
    CompressReads (Bloom<kmer_type>* bloom, Leon * leon, u_int64_t* nb_solids_kmers, int nb_cores, int * nbliving)
    : _bloom(bloom), _leon(leon),
    _total_nb_solid_kmers_in_reads (nb_solids_kmers), _local_nb_solid_kmers_in_reads(0),
    _synchro(System::thread().newSynchronizer()), _nb_living(nbliving)
    {
        _thread_id = __sync_fetch_and_add (_nb_living, 1);

        
    }
    
    ~CompressReads ()
    {
        /** We increase the global number of corrected errors. */
        _thread_id = __sync_fetch_and_add (_total_nb_solid_kmers_in_reads, _local_nb_solid_kmers_in_reads);
        
    }

    
    
    //copy construct
    CompressReads(const CompressReads& cr) //called by dispatcher iterate to create N functors
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
    
};




Leon::Leon () : Tool("leon"), _kmerSize(27) //, _inputBank (0)
{
    
    /** We get an OptionsParser for DSK. */
    OptionsParser parserDSK = DSK::getOptionsParser();
    getParser()->add (parserDSK);
    

    /** We add options specific to this tool. */

    
}


void Leon::execute ()
{
    

    /*************************************************/
    // We set some attributes (shortcuts).
    /*************************************************/
    _kmerSize           = getInput()->getInt (STR_KMER_SIZE);
    _solidFile          = getInput()->getStr (STR_KMER_SOLID);
    

    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    _bloom = createBloom ();
    LOCAL (_bloom);

    //iterate over initial file
    //BankFasta inbank (getInput()->getStr(STR_URI_FILE));
    IBank* inbank = BankRegistery::singleton().getFactory()->createBank(getInput()->getStr(STR_URI_FILE));


    /*************************************************/
    // We create a sequence iterator for the bank
    /*************************************************/
    Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          inbank->iterator(),
                                                          inbank->estimateNbSequences(),
                                                          "Iterating sequences"
                                                          );
    LOCAL (itSeq);
    
    /*************************************************/
    // We create the modified file
    /*************************************************/
    
    /** We get the basename from the provided URI (ie remove directory path and suffix). */
    string prefix = System::file().getBaseName (getInput()->getStr(STR_URI_FILE));
    

    /** We set the fileonme as the base name + a specific suffix. */
    string fileName;
    
    

    fileName = prefix + string(".leon");


        
    

    u_int64_t total_nb_solid_kmers_in_reads = 0;
    int nb_threads_living;
    

    /*************************************************/
    // We iterate over sequences and correct them
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "sequences correction");

        nb_threads_living = 0 ;
        
#ifdef SERIAL
        setDispatcher (new SerialDispatcher());
#else
        setDispatcher (  new Dispatcher (getInput()->getInt(STR_NB_CORES)) );
#endif
        

    
        getDispatcher()->iterate (itSeq,  CompressReads (_bloom , this, &total_nb_solid_kmers_in_reads,getInput()->getInt(STR_NB_CORES),&nb_threads_living),10000); // , 10000
        
    }
    

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    getInfo()->add (1, "result");
    getInfo()->add (2, "nb solid kmers in reads", "%ld", total_nb_solid_kmers_in_reads);
    
    
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
    BloomBuilder<> builder (estimatedBloomSize, 7,tools::collections::impl::BloomFactory::CACHE,getInput()->getInt(STR_NB_CORES));
    Bloom<kmer_type>* bloom = builder.build (itKmers);

    /** We return the created bloom filter. */
    return bloom;
}



