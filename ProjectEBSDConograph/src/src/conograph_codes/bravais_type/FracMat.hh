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
#ifndef FRACMAT_HH_
#define FRACMAT_HH_
#ifdef DEBUG
	#include <iostream>
#endif
#include "../HC_algorithm/GCD_1_3_1.hh"
#include "../utility_data_structure/nrutil_nr.hh"

template<class Integer>
class FracMat
{
public:
	NRMat<Integer> mat;
	Integer denom;

	FracMat(){};
	FracMat(const NRMat<Integer>& data) : mat(data), denom(1){};
	FracMat(const NRMat<Integer>& data, const Integer& i) : mat(data), denom(i){};
};


// Returns the inverse of rhs/denom
template<class Integer>
inline FracMat<Integer> FInverse3(const FracMat<Integer>& rhs)
{
	const NRMat<Integer>& imat = rhs.mat;
	assert( rhs.mat.nrows() == 3 && rhs.mat.ncols() == 3 );

	const Integer& denom = rhs.denom;

	const Integer det01 = imat[0][0]*imat[1][1]-imat[0][1]*imat[1][0];
	const Integer det02 = imat[0][0]*imat[2][2]-imat[0][2]*imat[2][0];
	const Integer det12 = imat[1][1]*imat[2][2]-imat[1][2]*imat[2][1];
	const Integer det01_02 = imat[0][0]*imat[1][2]-imat[0][2]*imat[1][0];
	const Integer det02_01 = imat[0][0]*imat[2][1]-imat[0][1]*imat[2][0];
	const Integer det01_12 = imat[0][1]*imat[1][2]-imat[0][2]*imat[1][1];
	const Integer det02_12 = imat[0][1]*imat[2][2]-imat[0][2]*imat[2][1];
	const Integer det12_02 = imat[1][0]*imat[2][2]-imat[1][2]*imat[2][0];
	const Integer det12_01 = imat[1][0]*imat[2][1]-imat[1][1]*imat[2][0];

	const Integer det = imat[0][0]*det12 - imat[0][1]*det12_02 + imat[0][2]*det12_01;
	assert( det != 0 );
	const Integer cdiv = GCD_Euclid(det, denom);

	NRMat<Integer> ans(3,3);
	ans[0][0] = det12;
	ans[0][1] = -det02_12;
	ans[0][2] = det01_12;
	ans[1][0] = -det12_02;
	ans[1][1] = det02;
	ans[1][2] = -det01_02;
	ans[2][0] = det12_01;
	ans[2][1] = -det02_01;
	ans[2][2] = det01;

	if( det > 0 ) return FracMat<Integer>(ans * (denom / cdiv), det/cdiv);
	else return FracMat<Integer>(ans * (-denom / cdiv), -det/cdiv);
}


template<class T>
static NRMat<T> operator/(const NRMat<T>& lhs, const T& rhs)
{
	const int nn = lhs.nrows();
	const int mm = lhs.ncols();
	NRMat<T> ans = lhs;
	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++)
		{
			assert( ans[k][j] % rhs == 0 );
			ans[k][j] /= rhs;
		}
	return ans;
}

template<class Integer>
inline NRMat<Integer> FIInverse3(const FracMat<Integer>& rhs)
{
	FracMat<Integer> ans = FInverse3(rhs);

	assert( ans.mat[0][0] % ans.denom == 0
			&& ans.mat[0][1] % ans.denom == 0
			&& ans.mat[0][2] % ans.denom == 0
			&& ans.mat[1][0] % ans.denom == 0
			&& ans.mat[1][1] % ans.denom == 0
			&& ans.mat[1][2] % ans.denom == 0
			&& ans.mat[2][0] % ans.denom == 0
			&& ans.mat[2][1] % ans.denom == 0
			&& ans.mat[2][2] % ans.denom == 0 );

	return ans.mat / ans.denom;
}


#endif /* FRACMAT_HH_ */
