# Leon 

| **Linux** | **Mac OSX** | **Functional tests** |
|-----------|-------------|----------------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-build-macos-10.9.5-gcc-4.2.1/) | [![Build Status](https://ci.inria.fr/gatb-core/view/Leon/job/tool-leon-functional-tests/badge/icon)](https://github.com/GATB/gatb-core/tree/master/gatb-core/test/jenkins/leon) |

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is Leon?

Leon is a software to compress Next Generation Sequencing data. It can compress Fasta or Fastq format; plain text and gzipped files are supported.

The method does not require any reference genome, instead a reference is built de novo from the set of reads as a probabilist de Bruijn Graph. It uses the disk streaming k-mer counting algorithm contained in the [GATB-Core library](https://github.com/GATB/gatb-core), and inserts solid k-mers in a bloom-filter. Each read is then encoded as a path in this graph, storing only an anchoring kmer and a list of bifurcations indicating which path to follow in the graph if several are possible.

G. Benoit, C. Lemaitre, D. Lavenier, E. Drezen, T. Dayris, R. Uricaru, G. Rizk. (2015) [Reference-free compression of high throughput sequencing data with a probabilistic de Bruijn graph](http://www.biomedcentral.com/1471-2105/16/288). BMC Bioinformatics, 2015, 16:288.
								
# Getting the latest source code

## Requirements

CMake 3.1+; see http://www.cmake.org/cmake/resources/software.html

c++/11 compiler; compilation was tested with gcc and g++ version>=4.8 (Linux) and clang version>=3.6 (Mac OSX).

## Instructions

    # get a local copy of source code
    git clone --recursive https://github.com/GATB/leon.git
    
    # compile the code an run a simple test on your computer
    cd leon
    sh INSTALL

# MANUAL	 
								
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

this will allow to use k<63. For larger k, change the values such that they are multiple of 32 and k1&lt;k2&lt;k3&lt;k4.

# CHANGELOG

* version 1.1.0: July 13, 2017:

    + Leon compressor integrated into GATB-Core 1.4.0
    + Leon compressed data is now stored into a single HDF5 binary file: '.leon'
    + bug fixes on compression/decompression on small reads
    + Extensive test suite added; [read more](https://github.com/GATB/gatb-core/tree/master/gatb-core/test/jenkins/leon)
    + Leon now successfully works on data sets reported in [Y. Zhang et al, 2017](https://github.com/GATB/gatb-core/tree/master/gatb-core/test/jenkins/leon) as breaking Leon
 
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

# Contact

To contact a developer, request help, etc: https://gatb.inria.fr/contact/