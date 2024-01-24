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
#include <ctime>
#include "Indexing.hh"
#include "ControlParam.hh"
#include "ProfileData.hh"
#include "conograph_codes/utility_func/zstring.hh"
#include "marquardt_codes/zerror_type/error_out.hh"
#include "zlog/zlog.hh"

using namespace std;

int main()
{
	clock_t start;
	start = clock();    /* Record the starting time. */
	try{
    CRLog::append(new CCoutListner());
    CRLog::append(new FileoutListner("LOG_CONOGRAPH.txt", zListnerID(1)));

	ControlParam cdata;
	{
ZLOG_INFO( "Reading input.txt ...\n" );
		const string& str = cdata.readFile("input.txt");
		if( str != "" )
		{
			// Some error happened.
			ZLOG_ERROR( str + "\n" );
			return 0;
		}
	}

	ProfileData pdata;
	{
ZLOG_INFO( "Reading data.txt ...\n" );
		const string& str = pdata.readFile("data.txt");
		if( str != "" )
		{
			// Some error happened.
			ZLOG_ERROR( str + "\n" );
			return 0;
		}
	}

	Indexing::conditions cons;
	{
		cons.search_level = cdata.search_level;
		cons.fitProjectionCenterShift[0] = cdata.fitProjectionCenterShift[0]; // false:fix Delta_X, true:fit Delta_X.
		cons.fitProjectionCenterShift[1] = cdata.fitProjectionCenterShift[1]; // false:fix Delta_Y, true:fit Delta_Y.
		cons.fitProjectionCenterShift[2] = cdata.fitProjectionCenterShift[2]; // false:fix Delta_Z, true:fit Delta_Z. (低対称性でバンド幅を使わない場合、z方向はfixする必要がある.)

		cons.num_hkl = cdata.num_hkl; // 理論値を求めるとき生成するミラー指数の数.
		cons.max_hkl = cdata.max_hkl; // 理論値を求めるとき生成するミラー指数h,k,lの許容される絶対値の最大数.

		cons.toleranceOnBandWidth = cdata.toleranceOnBandWidth; // when s1, s2 are the unit-cell scale from distinct bands, (1/toleranceOnBandWidth)*s2 <= s1 <= toleranceOnBandWidth*s2 is allowed.
		cons.resolutionForBLD = cdata.resolutionForBLD;
		cons.resolutionForOutput = cdata.resolutionForOutput;

		cons.ErrorForAngles_radian = cdata.ErrorForAngles_degree*M_PI/180.0;

		cons.lower_threshold_for_FOM = cdata.lower_threshold_for_FOM;
		cons.ABC_AXIS = cdata.ABC_AXIS;
		cons.RH_AXIS = cdata.RH_AXIS;
		cons.output_tex_format = cdata.output_tex_format;
	}

	Indexing ind;
	ind.set(cons, pdata);
	ind.run(cons);
	ind.print("out.txt", cons);

ZLOG_INFO( "The program has finished in CPU time : " + num2str( ( (clock() - start)*10 / CLOCKS_PER_SEC )*0.1 ) + " [sec.]\n\n" );
	}
	catch(const ZErrorMessage& etype)
	{
		ZLOG_ERROR( etype.printErrorLog() );
		return 0;
	}

	return 1;
}
