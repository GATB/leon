#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <cstdio>
#include <string>
#include <sstream>
#include <list>
#include <stdio.h>
#include <sys/wait.h>
#include <cstdlib>
#include <ctime>
#include <ratio>
#include <chrono>

#include <gatb/gatb_core.hpp>
#include <sys/time.h>

#include "utilitiesTests.hpp"

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



//ask the graph if it contains the kmers in the vector "kmers", and
//fill the informations in "graphStats"
void graphContains(Graph graph, vector<string>* kmers, struct GraphStats* graphStats, size_t kmerSize, KmerModel* kmerModel);
void testGraphContains(string inputFileName, Graph graph, size_t kmerSize, KmerModel* kmerModel);