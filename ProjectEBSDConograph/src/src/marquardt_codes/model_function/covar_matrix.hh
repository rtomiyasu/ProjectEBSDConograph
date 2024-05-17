/*
 * The MIT License

   Conograph (powder auto-indexing program)

Copyright (c) <2012> <Ryoko Oishi-Tomiyasu, KEK>

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
#ifndef _COVAR_MATRIX_HH_
#define _COVAR_MATRIX_HH_

#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../../conograph_codes/utility_data_structure/SymMat.hh"
#include "../utility_func/linear_constraint.hh"

// Return the square value of the argument.
// ( If the argument is  negative, return 0.0.)
inline Double sqrt_d(Double r)
{
	if(r<=0.0) return 0.0;
	else return sqrt(r);
};

// Calculate the covariant matrix of the parameters from start'th to (end-1)'th
// using the covariant matrix of the parameters independently fit.
void calPartOfCovariantMatrix(const Mat_DP_constr& Cmat,
								const SymMat<Double>& covar,
								const Int4& start,
								const Int4& end,
								SymMat<Double>& covar_all
							);

void calPartOfCovariantMatrix(const Mat_DP_constr& Cmat,
								const SymMat<Double>& covar,
								const Vec_INT& index_tray,
								SymMat<Double>& covar_all
							);

// Copy the covariant matrix of the parameters from start'th to (end-1)'th.
void copyPartOfCovariantMatrix(const SymMat<Double>&, const Int4& start, const Int4& end, SymMat<Double>&);

// Copy the covariant matrix from index_tray[i]'th entries.
void copyPartOfCovariantMatrix(const SymMat<Double>&, const Vec_INT& index_tray, SymMat<Double>&);

#endif /*_COVAR_MATRIX_HH_*/
