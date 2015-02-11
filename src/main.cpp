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

//#include <gatb/gatb_core.hpp>

#include <Leon.hpp>

using namespace std;


/********************************************************************************/


int main (int argc, char* argv[])
{
    // We dump some information about the gatb core library
    std::cout << LibraryInfo::getInfo();

	bool compress =false;
	bool decompress = false;
	
    // We define a try/catch block in case some method fails
    try
    {
		if(argc < 2)
		{
			printf("Usage for leon\n");
			printf("    ./leon [-c|-d] -file filename [options] \n");
			printf("Options  \n");
			printf("    -file               :    input file (e.g. FASTA/FASTQ for compress or .leon file for decompress)  \n");
			printf("    -c                  :    compression  \n");
			printf("    -d                  :    decompression  \n");
			printf("    -nb-cores           :    number of cores (default is the available number of cores)  \n");
			//printf("    -verbose                  :    verbosity level  \n");
			//printf("    -help                  :    display help about possible options  \n");
			printf("Compression options  \n");
			printf("    -kmer-size          :    size of a kmer  (default 31)\n");
			printf("    -abundance          :    abundance threshold for solid kmers  (default inferred)\n");
			printf("    -lossless           :    switch to lossless compression for qualities (default is lossy. lossy has much higher compression rate, and the loss is in fact a gain. lossy is better)\n");
			
			//printf("    -max-disk                  :    display help about possible options  \n");
			//printf("    -max-memory                  :    max memory in MBytes (default 1000)  \n");
			//printf("    -out                  :    output file (if not set basename of the input file)  \n");
			printf("Examples : \n");
			printf("    ./leon -file read.fasta -c \n");
			printf("    ./leon -file read.leon  -d \n");

    
    
			return EXIT_FAILURE;
		}
		
		for (int i=1; i< argc; i++)
		{
			if (strncmp("-c",argv[i],2)  == 0)  {  compress   = true;  }
            if (strncmp("-d",argv[i],2)  == 0)  {  decompress = true;  }
		}
		
        Leon (compress, decompress).run (argc, argv);
    }
    
    catch (misc::impl::OptionFailure& e)
    {
        return e.displayErrors (cout);
    }
    
    catch (gatb::core::system::Exception& e)
    {

        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
    
}

