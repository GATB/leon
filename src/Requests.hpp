#ifndef __REQUESTS__HPP__
#define __REQUESTS__HPP__



#include <iostream>
#include <gatb/gatb_core.hpp>
#include "../thirdparty/gatb-core/gatb-core/src/gatb/gatb_core.hpp"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/kmer/impl/Model.hpp"
#include "../thirdparty/gatb-core/gatb-core/src/gatb/gatb_core.hpp"
#include <sys/time.h>

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

#include <string>
#include <vector>
#include <ctype.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

//#include "RangeCoder.hpp"
#include "Leon.hpp"
#include <bitset>
//#include "CompressionUtils.hpp"

using namespace std;
class Leon;
//====================================================================================
// ** AbstractHeaderCoder
//====================================================================================
class Requests
{
	public:
		//AbstractHeaderCoder(Leon* leon);
		Requests(IBank* inputBank, string outputFilename, Graph graph, Partition<kmer_count> & solidCollection);
		

		void printSignatures();

		IBank* _inputBank;
		Iterator<Sequence>* _itBank;
		std::vector<Iterator<Sequence>*> _itBanks;
		int _nbBanks;

		string _outputFilename;

		//Kmer<>::ModelCanonical::Iterator _itKmer;

		Graph _graph;

		u_int64_t _solidFileSize;
		u_int64_t _nb_kmers_infile;

		unsigned char* _signature_array;
		unsigned char* _color_array;
		
	protected:
		
};

#endif __REQUESTS__HPP__ 

