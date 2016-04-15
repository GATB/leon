# Leon 

| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

#What is Leon?
Leon is a software to compress Next Generation Sequencing data. It can compress Fasta or Fastq format.
The method does not require any reference genome, instead a reference is built de novo from the set of reads as a probabilist de Bruijn Graph. It uses the disk streaming k-mer counting algorithm contained in the GATB library, and inserts solid k-mers in a bloom-filter. Each read is then encoded as a path in this graph, storing only an anchoring kmer and a list of bifurcations indicating which path to follow in the graph if several are possible.

G. Benoit, C. Lemaitre, D. Lavenier, E. Drezen, T. Dayris, R. Uricaru, G. Rizk. (2015) [Reference-free compression of high throughput sequencing data with a probabilistic de Bruijn graph](http://www.biomedcentral.com/1471-2105/16/288). BMC Bioinformatics, 2015, 16:288.
								
# Getting the latest source code

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Instructions

    # get a local copy of source code
    git clone --recursive https://github.com/GATB/leon.git
    
    # compile the code an run a simple test on your computer
    cd leon
    sh INSTALL

#MANUAL	 
								
## Usage

Mandatory options:

    -file  (1 arg) :    input file (e.g. FASTA/FASTQ for compression or .leon file for decompression)  
    -c             :    compress  
    -d             :    decompress  


Optional parameters:
 
    -nb-cores: number of cores 
    (default is the maximum available number of cores)  


Optional Compression parameters:

    -kmer-size <int>  :    kmer size  (default 31), currently should be <32 (or need to recompile, see below)
    -abundance <int>  :    minimal abundance threshold for solid kmers  (default: automatic)
    -lossless         :    switch to lossless compression for qualities (default is lossy. lossy has much higher compression rate, and the loss is in fact a gain. lossy is better!)
    -seq-only         :    store dna sequence only, header and qualitiess are discarded, will decompress to fasta (same as -noheader -noqual)
    -noheader         :    discard header
    -noqual           :    discard quality scores


Examples : 

    leon -file data/toy.fasta -c 
       -> generates the file toy.leon
 
    leon -file toy.leon  -d 
       -> "restores" the file toy.fasta

 

Note:
 
In order to use k values larger than 31, recompilation is necessary.

In the sequence of commands given in the INSTALL file, change the command: 

    cmake ..

by 

    cmake -Dk1=64 -Dk2=96 -Dk3=128 -Dk4=162 ..

this will allow to use k<63. For larger k, change the values such that they are multiple of 32 and k1<k2<k3<k4

#CHANGELOG

* version 1.0.0: April 16, 2015:
 bug fixes

* version 0.3:  March 19, 2015
 added quality compression, and other optimizations

* version 0.2.1:  Dec 18, 2014:
 bug fixes

* version 0.2: Oct 31, 2014:
 major performance improvement (about ~ 2 times faster)

* version 0.1.2  Aug 10, 2014:
 initial public release

#Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/