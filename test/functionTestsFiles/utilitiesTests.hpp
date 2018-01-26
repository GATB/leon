#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <string>
#include <vector>

#include <gatb/gatb_core.hpp>


/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;


typedef kmer::impl::Kmer<>::ModelDirect KmerModel;
typedef kmer::impl::Kmer<>::Type        kmer_type;
typedef kmer::impl::Kmer<>::Count       kmer_count;

using namespace std;



kmer_type getKmerType(const char* kmer_chars, KmerModel* _kmerModel);
void getKmerChars(kmer_type kmer, char* kmer_chars);
string getKmerString(kmer_type kmer);
Node getKmerNode(char* kmer_chars);
Node getKmerNode(string kmer_string, KmerModel* _kmerModel);
Node getKmerNode(kmer_type kmer);


//return reverse complement of a string
string getReverseComplement(string kmer);
//fill the vector "kmers" with the kmers of size "kmerSize" from the string "read"
void getKmers(string read, vector<string>* kmers, int kmerSize);