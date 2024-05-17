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

#ifndef SRC_CONTROLPARAM_HH_
#define SRC_CONTROLPARAM_HH_

#include <vector>
#include <string>
#include "conograph_codes/bravais_type/enumAxis.hh"

using namespace std;

class ControlParam
{
public:
	eABCaxis ABC_AXIS;
	eRHaxis RH_AXIS;
	bool output_tex_format;

	// 0: quick search, 1: exhaustive search.
    int search_level;

	// [0]: false:fix Delta_X, true:fit Delta_X.
	// [1]: false:fix Delta_Y, true:fit Delta_Y.
	// [2]: false:fix Delta_Z, true:fit Delta_Z. (低対称性でバンド幅を使わない場合、z方向はfixする必要がある.)
    bool fitProjectionCenterShift[3];

    double ErrorForAngles_degree;
    size_t num_hkl; // 理論値を求めるとき生成するミラー指数の数.
    int max_hkl;    // 理論値を求めるとき生成するミラー指数h,k,lの許容される絶対値の最大数.
    double resolutionForBLD;    // 格子ベクトルそれぞれの長さ|l|に対し、|l|*resolの差異を許す。Used for Bravais-lattice determination.
    double resolutionForOutput; // 格子ベクトルそれぞれの長さ|l|に対し、|l|*resolの差異を許す。Used for selection output solutions.
    double toleranceOnBandWidth; // when s1, s2 are the unit-cell scale from distinct bands, (1/toleranceOnBandWidth)*s2 <= s1 <= toleranceOnBandWidth*s2 is allowed.
    double lower_threshold_for_FOM;
    bool output_txt_format;

    ControlParam();
    ~ControlParam();

    string readFile(const string& fname);
};



#endif /* SRC_CONTROLPARAM_HH_ */
