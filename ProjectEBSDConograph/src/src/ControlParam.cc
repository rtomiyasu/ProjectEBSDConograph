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

#include "ControlParam.hh"
#include "conograph_codes/utility_func/zstring.hh"
#include "zlog/zlog.hh"

static string printErrorLog(const string& fname, const int& lnumber, const string& funcname, const string& err_message)
{
   	stringstream os;
	os << fname + "(line " << lnumber << ") : " + funcname + " fails.\n"
		+ err_message + ".\n";
	return os.str();
}

ControlParam::ControlParam()
{
	ABC_AXIS = B_Axis;
	RH_AXIS = Rho_Axis;
	output_tex_format = true;
}

ControlParam::~ControlParam()
{
}

string ControlParam::readFile(const string& filename)
{
    ifstream ifs( filename.c_str() );
    if (!ifs){
    	return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot open file: "  + filename );
    }

	string s, str;

	ZLOG_INFO("***                     ***\n");
	ZLOG_INFO("*** Input parameter set ***\n");
	ZLOG_INFO("***                     ***\n");
    //  Reads the first header.
	int count = 0;
    while( count < 12 && getnewline(ifs, s) )
	{
		if(is_blank(s)) continue;
		if(s.at(0) == '#') continue;

		istringstream iss(s);
		if(count==0){
ZLOG_INFO( " Quick search/exhaustive search                                    : " );
			iss >> search_level;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <SearchLevel>: "  + filename );
		}
		else if(count==1){
ZLOG_INFO( " Upper bound on errors in phi, sigma, sigma_begin, sigma_end: " );
			iss >> ErrorForAngles_degree;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <CriticalValueForLinearSum>: "  + filename );
		}
		else if(count==2){
ZLOG_INFO( " Tolerance level for errors in the unit-cell scales         : " );
			iss >> toleranceOnBandWidth;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <ToleranceLevelOnBandWidth>: "  + filename );
		}
		else if(count==3){
ZLOG_INFO( " Resolution for Bravais-type determination                  : " );
			iss >> resolutionForBLD;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <Resolution>: "  + filename );
		}
		else if(count==4){
ZLOG_INFO( " Resolution for selecting output solutions                  : " );
			iss >> resolutionForOutput;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <Resolution>: "  + filename );
		}
		else if(count==5){
ZLOG_INFO( " Number of hkl generated  for indexing                      : " );
			iss >> num_hkl;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <MaxNumberOfPeaksForFOM>: "  + filename );
		}
		else if(count==6){
ZLOG_INFO( " Max |h|,|k|,|l| used for indexing                          : " );
			iss >> max_hkl;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <MaxAbsoluteHKLForFOM>: "  + filename );
		}
		else if(count==7){
ZLOG_INFO( " Flags to fit the projection center shifts                  : " );
			for(int i=0; i<3; i++){
				iss >> fitProjectionCenterShift[i];
				if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <FixProjectionShift>: "  + filename );
			}
		}
		else if(count==8){
ZLOG_INFO( " Lower threshold of FOM:                                    : " );
			iss >> lower_threshold_for_FOM;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <MinFigureOfMerit>: "  + filename );
		}
		else if(count==9){
ZLOG_INFO( " Axis for rhombohedral symmetry                             : " );
			iss >> str;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <AxisForRhombohedralSymmetry>: "  + filename );
			if( str == "Rhombohedral" ){
				RH_AXIS = Rho_Axis;
			}
			else if( str == "Hexagonal" ){
				RH_AXIS = Hex_Axis;
			}
			else{
				return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <AxisForRhombohedralSymmetry>: "  + filename );
			}
		}
		else if(count==10){
ZLOG_INFO( " Axis for monoclinic symmetry                               : " );
			iss >> str;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <AxisForMonoclinicSymmetry>: "  + filename );
			if( str == "A" ){
				ABC_AXIS = A_Axis;
			}
			else if( str == "B" ){
				ABC_AXIS = B_Axis;
			}
			else if( str == "C" ){
				ABC_AXIS = C_Axis;
			}
			else{
				return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <AxisForMonoclinicSymmetry>: "  + filename );
			}
		}
		else if(count==11){
ZLOG_INFO( " Output in latex style                                      : " );
			iss >> output_tex_format;
			if( iss.fail() ) return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read <OutputInTexStyle>: "  + filename );
		}
		count++;

ZLOG_INFO( s + "\n" );
	}
ZLOG_INFO( "\n" );

    if( count < 11 )
    {
		return printErrorLog(__FILE__, __LINE__, __FUNCTION__, "Cannot read " + num2str(count+1) + "'th parameter: "  + filename );
    }

    return "";
}
