//
//  Leon.hpp
//  leon
//
//  Created by Guillaume Rizk on 17/01/2014.
//  Copyright (c) 2014 G.Rizk, R.Uricaru. All rights reserved.
//

#ifndef __leon__Leon__
#define __leon__Leon__

#include <iostream>
#include <gatb/gatb_core.hpp>

#include <string>
#include <sstream>


/** NOTE: we should not include namespaces here => only to make user life easier... */
using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;


typedef kmer::impl::Kmer<>::Model KmerModel;
typedef kmer::impl::Kmer<>::Type  kmer_type;
typedef kmer::impl::Kmer<>::Count kmer_count;

char char2phred(char c);
double char2proba(char c);



class Leon : public misc::impl::Tool
{
private:
 
    static const char* STR_GZ;

    void execute ();

    collections::impl::Bloom<kmer_type>* _bloom;

    virtual collections::impl::Bloom<kmer_type>* createBloom ();
    
public:
    
    size_t          _kmerSize;
    std::string     _solidFile;
    Leon();

};
#endif /* defined(__leon__Leon__) */
