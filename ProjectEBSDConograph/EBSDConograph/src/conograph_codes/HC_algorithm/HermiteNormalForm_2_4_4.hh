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

#ifndef _HermiteNormalForm_2_4_4_HH_
#define _HermiteNormalForm_2_4_4_HH_

#include "../utility_data_structure/nrutil_nr.hh"

using namespace std;
// Algorithm Hermite Normal Form,  2.4.4 page# 68
/*Algorithm 2.4.4 (Hermite Normal Form). Given an m x n matrix A with
integer coefficients (a;,j) this algorithm finds the Hermite normal form W of A.
As usual, we write w;,j for the coefficients of W, A; (resp. Wi) for the columns
of A (resp. W).*/

template<class Fraction, class Integer>
void Hermite_Normal_Form(vector< NRVec<Fraction> > A, vector< NRVec<Fraction> >& W)
{
	if( A.empty() )
	{
		W.clear();
		return;
	}

    const int n = A.size();
    const int m = A[0].size();

    int j0, j;
    Fraction b, bb;
    Integer q;
    const int l = (m > n?m - n:0);
    
    for(int i = m - 1, k = n - 1; ; )
    {
      	 for (j=0; j<k; j++)
      		 if ( A[j][i] != 0) break;

      	 if (j>=k)
      	 {
      		 if (A[k][i] < 0)
      		 {
        		 for (int jj= 0; jj<m; jj++) A[k][jj] *= -1;
      		 }

      			// step5:
      	   	   	 b = A[k][i];
      	   	   	 if ( b == 0 ){
      	   	   		 k++;
      	   	   	 }
      	   	   	 else
      	   	   	 {
      	   	   		 for (j=k+1; j<n; j++)
      	   	   		 {
      	   	   			 q = (A[j][i] / b).floor();
      	   	   			 for (int jj=0; jj<m; jj++)
      	   	   			      A[j][jj] -= A[k][jj] * q;
      	   	   		 }
      	   	   	 }

      	   	 //step6:
      	   	 if (i == l)
      	   	 {
      	   		 W.resize(n-k, NRVec<Fraction>(m));
      	   		 for (int j=0; j<n-k; j++)
      	   			 for (int jj= 0; jj<m; jj++)
      	   				 W[j][jj] = A[j+k][jj];
      	   		 return;
      	   	 }
      	   	 else {i--;k--;continue;}
      	 }

         // step3:
      	 if (A[j][i] != 0)
      	 {
         	 bb = abs(A[j][i]);
         	 j0 = j;
      	 }
      	 for (j=0; j<=k; j++)
      	 {
      		 if (A[j][i] != 0)
      		 {
      			 if ( abs(A[j][i]) <= bb ) { bb = abs(A[j][i]); j0 = j; }
      		 }
      	 }
      	 if ( j0 < k ){
      		 for (int jj= 0; jj<m; jj++)
      		 {
      		    swap( A[k][jj], A[j0][jj]);
      		 }
      	 }
      	 if ( A[k][i] < 0 ){
      		 for (int jj= 0; jj<m; jj++)
      			 A[k][jj] *= -1;
      	 }
      	 b = A[k][i];

      	 //step4:
      	 for (j=0; j<=k-1; j++)
      	 {
     		 q = ( A[j][i] / b ).lrint();
     		 for (int jj= 0; jj<m; jj++)
      		      A[j][jj] -= A[k][jj] * q;
      	 }
    }
}

#endif
