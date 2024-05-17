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
#ifndef _SVdcmp_h_
#define _SVdcmp_h_
// SVdcmp.hh

#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../../conograph_codes/utility_data_structure/nrutil_nr.hh"
#include "../../conograph_codes/utility_data_structure/SymMat.hh"

using namespace std;
 
/*
	Decompose a symmetric matrix alpha into transpose(V) * T * V using matrices as followings.
	T is a diagonal matrix. V is a matrix such that	
	transpose(V) * V is a diagonal matrix whose j'th diagonal element equals to the input alpha[j][j].  
    On Output, alpha is replaced by V^-1.
    adiag_in[j] equal the j'th diagonal element of T.
*/   
void SVdcmp_precondition(NRMat<Double>&, NRVec<Double>&, const Int4& ibegin, const Int4& iend);

// Calculate the inverse matrix of transpose(V^-1) * adiag * V^-1. (adiag is a diagonal matrix.)
void getInverseMatrix(const NRMat<Double>& V, const NRVec<Double>& adiag,
						const Int4& ibegin, const Int4& iend,
						const Double& thred, SymMat<Double>& InvMat);

// Singlular value decomposition.
void SVdcmp(NRMat<Double>&, NRVec<Double>&, const Int4& ibegin, const Int4& iend);
void tred2(NRMat<Double>&, NRVec<Double>&, NRVec<Double>&, const Int4& ibegin, const Int4& iend);
void tqli(NRMat<Double>&, NRVec<Double>&, NRVec<Double>&, const Int4& ibegin, const Int4& iend);

template<class T>
inline const T SQR(const T a){return a*a;};

template<class T>
inline const T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);};

inline Double pythag(const Double a, const Double b)
{
	Double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
};

void putMatrixDiagonalAscendingOrder(const NRMat<Double>& alpha, Vec_INT& row_order);

void permSquareMatrixElement(const Vec_INT& row_order, NRMat<Double>& alpha);

void permSquareMatrixRow(const Vec_INT& row_order, NRMat<Double>& alpha);

#endif
