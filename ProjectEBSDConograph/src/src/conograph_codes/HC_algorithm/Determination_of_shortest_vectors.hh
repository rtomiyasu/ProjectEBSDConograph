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

#ifndef DETERMINATION_OF_SHORTEST_VECTORS_CC_
#define DETERMINATION_OF_SHORTEST_VECTORS_CC_

#include <algorithm>
//#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "../utility_data_structure/nrutil_nr.hh"
#include "Choleskydcmp_pos_def.hh"
using namespace std;


// On output, A satisfies the following equality for any (x1, ..., xm):
// | x1*b1 + ... + xm*bm |^2 = \sum_{i=1}^m A[i][i]*(xi + \sum_{j=i+1}^m aij*xj)^2.
template<class T>
void make_aij(SymMat<T> Qmat, NRMat<T>& A)
{
	const int m = Qmat.size();
	assert( m == A.nrows() );
	assert( m == A.ncols() );

	//For  i = 1,...,n-1  set  aii = bii, aik = bik/aii  (k  =  i+l .... ,n)
	for(int i=0; i<m; i++){
		A[i][i] = Qmat(i,i);
		for(int k=i+1; k<m; k++){
			A[i][k] = Qmat(i,k) / A[i][i];
		}
		for(int j=i+1; j<m; j++){
			for(int k=j; k<m; k++){
				Qmat(j,k) -= A[i][j]*Qmat(i,k);
			}
		}
	}
}


// Input:
// Qmat(*, *): Gram matrix.
// beta_i (1<i<m): real number.
// On output, in shortest_vec_a, the vectors (n_1, ..., n_m)
// of length | (beta_1+n_1)*b_1 + ... + (beta_1+n_m)*b_m |^2 < Max_length are saved.
template<class T>
bool get_short_vectors(const SymMat<T>& Qmat,
		const NRVec<T>& beta, const T& Max_length, vector< pair<NRVec<int>, T> >& shortest_vec_a)
{
	assert( Qmat.size() > 0 );

	const int m = Qmat.size();
//	const int n = b[0].size();

	NRMat<double> Amat(m,m);
	make_aij(Qmat, Amat);

	shortest_vec_a.clear();

	// SM(i) = beta_i + \sum_{j=i+1}^m aij * (beta_j + m_j)
	// SL(i) = beta_i + m_i + \sum_{j=i+1}^m aij * (beta_j + m_j)
	NRVec<T> SM(m), SZ(m);

	// step1
//	NRVec<T> vec_b(n, 0);
	for(int j = 0; j < m; j++)
	{
		SM[j] = beta[j];
//		vec_b += b[j]*beta[j];
	}

	SZ[m-1] = Max_length;

//	NRVec<T> SL(m);
	NRVec<T> lower_bound(m), upper_bound(m);
	NRVec<int> vec_m(m);
	for(int i=m-1; i>=0; i--)
	{
		// step2
		const T t = sqrt(SZ[i] / Amat[i][i]);
		lower_bound[i] = T(ceil(-t-SM[i]));
		upper_bound[i] = T(floor(t-SM[i]));

		for(vec_m[i]=lower_bound[i]; ; vec_m[i]++)
		{
			// Step3
			if( vec_m[i] > upper_bound[i] )
			{
				i++;
				if( i < m ) continue;
				else return true;
			}

			// step4
			if(i != 0) break;

			// step7
			double length = 0.0;
			for(int j=0; j<m; j++)
			{
				length += Qmat(j,j) * (beta[j] + vec_m[j]) * (beta[j] + vec_m[j]);
				for(int j2=j+1; j2<m; j2++)
				{
					length += Qmat(j,j2) * (beta[j] + vec_m[j]) * (beta[j2] + vec_m[j2]) * 2.0;
				}
			}
			if(length > Max_length) continue;
			shortest_vec_a.push_back(make_pair(vec_m, length));

			// Shortest vector�����߂�Ƃ��͈ȉ����s���B
//			for(int j=0; j<m; j++)
//			{
//				SZ[j] -= SZ[m-1]-inner_product(vec_c,vec_c);
//			}
//			assert( SZ[m-1] >= T(0) );

			i = 0;
			// step8
			while(SZ[i] < T(0))
			{
				i++;
				assert( i < m );
			}
		}

		// step5
		const T SLi = vec_m[i] + SM[i];
		SZ[i-1] = SZ[i] - SLi * SLi * Amat[i][i];
		SM[i-1] = beta[i-1];
		for(int j=i; j<m; j++){
			SM[i-1] += Amat[i-1][j] * (vec_m[j] + beta[j]);
		}
	}
	return false;
}

template<class T>
inline bool operator<(const pair<NRVec<int>, T>& lhs, const pair<NRVec<int>, T>& rhs)
{
	return lhs.second < rhs.second;
}

template<class T>
bool get_short_vectors(const SymMat<T>& Qmat,
		const T& Max_length, vector< pair<NRVec<int>, T> >& shortest_vec_a)
{
	return get_short_vectors(Qmat, NRVec<T>(Qmat.size(), 0.0), Max_length, shortest_vec_a);
}

/*
int main ()
{
	vector< NRVec<double> > b(3, NRVec<double>(3));
    const int m = b.size();
    const int n = b[0].size();
    b[0][0] = 0;
	b[0][1] = 1;
	b[0][2] = 1;
//	b[0][3] = 0;
//	b[0][4] = 0;
//	b[0][5] = 0;
//	b[0][6] = 0;
	b[1][0] = 1;
	b[1][1] = 0;
	b[1][2] = 1;
//	b[1][3] = 0;
//	b[1][4] = 0;
//	b[1][5] = 0;
//	b[1][6] = 0;
	b[2][0] = 1;
	b[2][1] = 1;
	b[2][2] = 0;
//	b[2][3] = 0;
//	b[2][4] = 0;
//	b[2][5] = 0;
//	b[2][6] = 0;
//	b[3][0] = 0;
//	b[3][1] = 0;
//	b[3][2] = 0;
//	b[3][3] = 1;
//	b[3][4] = 0;
//	b[3][5] = 0;
//	b[3][6] = 0;
//	b[4][0] = 0;
//	b[4][1] = 0;
//	b[4][2] = 0;
//	b[4][3] = 0;
//	b[4][4] = 1;
//	b[4][5] = 0;
//	b[4][6] = 0;
//	b[5][0] = 0;
//	b[5][1] = 0;
//	b[5][2] = 0;
//	b[5][3] = 0;
//	b[5][4] = 0;
//	b[5][5] = 1;
//	b[5][6] = 0;
//	b[6][0] = 0;
//	b[6][1] = 0;
//	b[6][2] = 0;
//	b[6][3] = 0;
//	b[6][4] = 0;
//	b[6][5] = 0;
//	b[6][6] = 1;
//
	NRVec<double> beta(m,0.);
	beta[0] = 0.0;
	beta[1] = 0.0;
	beta[2] = 0.0;
//	beta[3] = 0.5;
//	beta[4] = 0.5;
//	beta[5] = 0.5;
//	beta[6] = 0.5;

	vector< pair<NRVec<double>, double> > ans;

	if( !Determination_of_shortest_vectors(b, beta, 6., ans) )
	{
		cout << "The rows of b are linearly dependent.\n";
	}
	sort(ans.begin(), ans.end());

	 cout<<"shortest vector a: "<<endl;
	 for (size_t i=0; i<ans.size(); i++)
	 {
		 cout << "(" << ans[i].first[0];
		 for (int k=1; k<m; k++)
		 {
			 cout << ", " << ans[i].first[k];
		 }
		 cout << "), " << ans[i].second << "\n";
	 }
	 cout<<endl;

	return 0;
}
*/
#endif /* DETERMINATION_OF_SHORTEST_VECTORS_CC_ */
