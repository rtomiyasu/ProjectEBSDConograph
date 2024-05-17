/*
 * The MIT License

   EBSD-Conograph (Conograph software for EBSD ab-initio indexing)

Copyright (c) <2020> <Ryoko Oishi-Tomiyasu, Kyushu University>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 *
 */

#include <cmath>
#include <fstream>
#include <sstream>

#include "ProfileData.hh"
#include "conograph_codes/utility_func/zstring.hh"
#include "zlog/zlog.hh"

static string printErrorLog(const string& fname, const int& lnumber, const string& funcname, const string& err_message)
{
   	stringstream os;
	os << fname + "(line " << lnumber << ") : " + funcname + " fails.\n"
		+ err_message + ".\n";
	return os.str();
}

ProfileData::ProfileData() 
{
}

ProfileData::~ProfileData()
{
}

string ProfileData::readFile(const string& filename)
{
	m_data.clear();

    ifstream ifs( filename.c_str() );
    if (!ifs){
    	return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot open file: "  + filename );
    }


	string s;

    //  Reads the first header.
    while( getnewline(ifs, s) )
	{
		if(is_blank(s)) continue;
		if(s.at(0) == '#') continue;

		istringstream iss(s);
		iss >> get_unitcell_scale;
		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the 1st flag: "  + filename );

		iss >> wavelength;
		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the wave length: "  + filename );
		break;
	}

    //  Reads the second header.
	DataSet data;
	data.sigma_begin = 0.0;
	data.sigma_end = 0.0;

	static const size_t MaxNumColumns = 4;
	size_t NumColumns = 0;
	double tray[MaxNumColumns];
	while( getnewline(ifs, s) )
	{
    	if(is_blank(s)) continue;
    	if(s.at(0) == '#') continue;

    	istringstream iss(s);

    	for(NumColumns=0; NumColumns<MaxNumColumns; NumColumns++)
    	{
        	iss >> tray[NumColumns];
       		if( iss.fail() ) break;
    	}
    	if( NumColumns < 2 || NumColumns > 4 )
    	{
   	   		return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "The number of columns must be 2--4: "  + filename );
    	}

    	int index=0;
    	data.phi = tray[index++];

   		if( NumColumns == 2 || NumColumns == 4 )
   		{
   	   		data.sigma = tray[index++];
   	   		if( data.sigma < 0.0 || data.sigma >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma is not in [0, PI/2): "  + filename );
   		}
   		if( NumColumns == 3 || NumColumns == 4 )
   		{
	   		data.sigma_begin = tray[index++];
	   		if( data.sigma_begin >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma_begin is not [0, PI/2): "  + filename );

	   		data.sigma_end = tray[index++];
	   		if( data.sigma_end <-data.sigma_begin ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "sigma_begin + sigma_end < 0: "  + filename );
	   		if( data.sigma_end < data.sigma_begin ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "sigma_begin > sigma_end: "  + filename );
	   		if( data.sigma_end >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma_end is not in [0, PI/2): "  + filename );

	   		if( NumColumns == 3 ) data.sigma = ( data.sigma_begin + data.sigma_end )*0.5;
   		}

   		m_data.push_back( data );
   		break;
	}

ZLOG_INFO( "Number of columns = " + num2str(NumColumns) + ".\n\n" );

	while( getnewline(ifs, s) )
	{
    	if(is_blank(s)) continue;
    	if(s.at(0) == '#') continue;

    	istringstream iss(s);

   		iss >> data.phi;
   		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the 1st column (phi): "  + filename );

   		if( NumColumns == 2 || NumColumns == 4 )
   		{
   	   		iss >> data.sigma;
   	   		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the 2nd column (sigma): "  + filename );
   	   		if( data.sigma < 0.0 || data.sigma >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma is not in [0, PI/2): "  + filename );
   		}
   		if( NumColumns == 3 || NumColumns == 4 )
   		{
	   		iss >> data.sigma_begin;
	   		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the 3rd column (sigma_begin): "  + filename );
	   		if( data.sigma_begin >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma_begin is not [0, PI/2): "  + filename );

	   		iss >> data.sigma_end;
	   		if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read the 4th column (sigma_end): "  + filename );
	   		if( data.sigma_end <-data.sigma_begin ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "sigma_begin + sigma_end < 0: "  + filename );
	   		if( data.sigma_end < data.sigma_begin ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "sigma_begin > sigma_end: "  + filename );
	   		if( data.sigma_end >= 90.0 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Some sigma_end is not in [0, PI/2): "  + filename );

	   		if( NumColumns == 3 )
	   		{
	   			data.sigma = ( data.sigma_begin + data.sigma_end )*0.5;
	   		}
   		}

   		m_data.push_back( data );
	}
    if( m_data.size() < 5 ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Number of input data is too small: "  + num2str(m_data.size()) );

    static const Double RadDeg = M_PI /180.0;
    for(vector<DataSet>::iterator it=m_data.begin(); it<m_data.end(); it++)
    {
   		it->phi *= RadDeg;
   		it->sigma *= RadDeg;
   		it->sigma_begin *= RadDeg;
   		it->sigma_end *= RadDeg;
    }
    return "";
}
