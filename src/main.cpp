/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  R.Chikhi, G.Rizk, E.Drezen
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
#include <DSK.hpp>

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
        misc::impl::ToolComposite tool;
		
		
		for (int i=1; i< argc; i++) {
			if( strncmp("-c",argv[i],2)  == 0)
			{
				tool.add (new DSK  () );
				compress = true;
			}
		}
		
		for (int i=1; i< argc; i++) {
			if( strncmp("-d",argv[i],2)  == 0)
			{
				//tool.add (new DSK  () );
				decompress = true;
			}
		}
		
        tool.add (new Leon    (compress, decompress));
        
        tool.run (argc, argv);
    }
    
    catch (misc::impl::OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }
    
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
    
}

