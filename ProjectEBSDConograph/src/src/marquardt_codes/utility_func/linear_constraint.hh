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

#ifndef LINEAR_CONSTRAINT_HH_
#define LINEAR_CONSTRAINT_HH_

#include<string>
#include "../../conograph_codes/utility_data_structure/index_set.hh"
#include "../zparam/etype_ID.hh"

// Vec_DP_save are used as a container of vectors for the saving memory.
typedef vector< index_set<Double> > Vec_DP_save;

struct constr_DP{
	etype_ID ID;
	Vec_DP_save constr;
};

// Mat_DP_save are used as a container of flags and constraints for the saving memory.
typedef vector< constr_DP > Mat_DP_constr;


// Returns the inner product.(lhs, rhs are vectors.)
inline Double product_constr(const Vec_DP_save& lhs, const Double* rhs)
{
	Double t = 0.0;
	for(UInt4 l=0; l<lhs.size(); l++) t += rhs[lhs[l].index] * lhs[l].element;
	return t;
}

// Change the type of a vector from Vec_DP_save to Vec_DP.
inline void chToVector(const Vec_DP_save& lhs, const Int4 pn, Vec_DP& rhs)
{
 	rhs.clear();
 	rhs.resize(pn, 0.0);
 	for(UInt4 j=0; j<lhs.size(); j++) rhs[lhs[j].index] += lhs[j].element;
}

static const double thred0 = 1.0e-8;

void gauss_elimination_row(vector< vector<double> >&, const double& thred = thred0);
void putExplicitConstraints(const size_t&, vector< vector<double> >&, Mat_DP_constr&, const double& thred = thred0);

// First call the function gauss_elimination_row before using the following function.
bool checkIfFirstConstraintsSatisfySecond(const vector< vector<double> >&, const vector< vector<double> >&, const double& thred = thred0);

string chEquationToString(const int index, const Vec_DP_save& cont, const vector<string>& par_label);

#endif /*LINEAR_CONSTRAINT_HH_*/
