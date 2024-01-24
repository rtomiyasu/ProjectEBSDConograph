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

#ifndef _MatrixInverse_2_2_7_HH_
#define _MatrixInverse_2_2_7_HH_

#include "../utility_data_structure/nrutil_nr.hh"

// Algorithm 2.2.2 page# 50
template<class T>
bool Inverse_Matrix(NRMat<T> m,  NRMat<T>& x)
{
    const int n = m.nrows();
    assert( n == m.ncols() );
    T c[n];
    NRMat<T> B(n, n, T(0));
	for (int j=0;j<n;j++) B[j][j] = T(1);

	int i;
	for (int j=0;j<n;j++)
	{
		for (i=j;i<n;i++)
			if ( m[i][j] != 0 ) break;
		if( i >= n )
		{
			return false;
		}

		if( i > j )
		{
			for (int l=j;l<n;l++)
			{
				swap(m[i][l], m[j][l]);
			}
			for (int l=0;l<n;l++)
			{
				swap(B[i][l], B[j][l]);
			}
		}

		T d = T(1) / m[j][j];

		for(int k=j+1; k<n; k++)
		{
			c[k] = m[k][j] * d;
		}

		for(int k=j+1; k<n; k++)
		{
			for (int l=j+1; l<n; l++)
				m[k][l] -= c[k] * m[j][l];

			for (int l=0;l<n;l++)
				B[k][l] -= c[k] * B[j][l];
		}

	}

	NRVec<T> sigma(n);
	x = NRMat<T>(n, n);
	for (int i=n-1; i>=0; i--)
	{
		sigma = T(0);
		for (int j=i+1;j<n;j++)
		{
			for (int ii=0; ii<n; ii++)
			{
				sigma[ii] += m[i][j] * x[j][ii];
			}
		}

		for (int ii=0; ii<n; ii++)
		{
			x[i][ii] = ( B[i][ii] - sigma[ii] ) / m[i][i];
		}
	}

	return true;


}

//int main ()
//{
////	NRMat<double> m(2, 2);
////	NRMat<double> x(2, 2);
////
////	m[0][0] = 1.;
////	m[0][1] = 2.;
////	m[1][0] = 3.;
////	m[1][1] = 4.;
//
//	NRMat<double> m(3, 3);
//	NRMat<double> x(3, 3);
//	const int n = m.nrows();
//
//	m[0][0] =  2.;
//	m[0][1] =  1. ;
//	m[0][2] = -1. ;
//	m[1][0] = -3. ;
//	m[1][1] = -1.;
//	m[1][2] =  2.;
//	m[2][0] = -2. ;
//	m[2][1] =  1.;
//	m[2][2] =  2.;
//
//	if( !Inverse_Matrix(m,x) )
//	{
//		cout<<" M is not invertible. "<<endl;
//	}
//
//	for (int j=0; j<n; j++)
//			for (int i=0;i<n;i++)
//				std::cout<<"x ["<<j<<"]["<<i<<"]= "<<x[j][i]<<std::endl;
//
//
//	return 0;
//}

#endif
