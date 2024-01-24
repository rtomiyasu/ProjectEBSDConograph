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
#ifndef _TRANSFORM_SYM_MATRIX_HH_
#define _TRANSFORM_SYM_MATRIX_HH_

#include"../RietveldAnalysisTypes.hh"
#include"../zerror_type/error_out.hh"
#include"../utility_data_structure/nrutil_nr.hh"
#include"../utility_data_structure/SymMat.hh"

// Returns lhs * rhs * transpose(lhs), where rhs is a symmetric matrix.
template <class U, class T>
inline SymMat<T> transform_sym_matrix(const NRMat<U>& lhs, const SymMat<T>& rhs)
{
	const Int4 icol = lhs.ncols();
	const Int4 isize = lhs.nrows();

	assert( rhs.size() == icol );
	
	SymMat<T> ans(isize, 0);
	for(int i=0; i<isize; i++)
	{
		for(int i2=i; i2<isize; i2++)
		{
			for(int j=0; j<icol; j++)
			{
				T& ans_r = ans(i,i2);
				ans_r += rhs(j,j)*( lhs[i][j]*lhs[i2][j] );
				for(int j2=j+1; j2<icol; j2++)
				{
					ans_r += rhs(j,j2)*( lhs[i][j]*lhs[i2][j2] + lhs[i][j2]*lhs[i2][j] );
				}
			}
		}
	}
	return ans;
}


// Returns lhs * rhs * transpose(lhs), where rhs is a symmetric matrix.
template <class U, class T>
inline SymMat<T> transform_sym_covar(const NRMat<U>& lhs, //const SymMat<T>& rhs,
		const SymMat<T>& rhs_var)
{
//	const SymMat<T> ans = transform(lhs, rhs);

	const Int4 icol = lhs.ncols();
	const Int4 isize = lhs.nrows();
	
	const Int4 icol2 = icol*(icol+1)/2;
	const Int4 isize2 = isize*(isize+1)/2;
	
	assert( rhs_var.size() == icol2 );
	SymMat<T> ans_var( isize2 );
	
	NRMat<Int4> trans_mat_var(isize2, icol2);
	
	for(Int4 i=0, index=0; i<isize; i++)
	{
		for(Int4 i2=0, index2=0; i2<icol; i2++)
		{
			trans_mat_var[index][index2] = lhs[i][i2]*lhs[i][i2];
			index2++;
			for(Int4 j2=i2+1; j2<icol; j2++, index2++)
				trans_mat_var[index][index2] = (lhs[i][i2]*lhs[i][j2]+lhs[i][j2]*lhs[i][i2]);
		}
		index++;
		for(Int4 j=i+1; j<isize; j++, index++)
		{
			for(Int4 i2=0, index2=0; i2<icol; i2++)
			{
				trans_mat_var[index][index2] = lhs[i][i2]*lhs[j][i2];
				index2++;
				for(Int4 j2=i2+1; j2<icol; j2++, index2++)
					trans_mat_var[index][index2] = (lhs[i][i2]*lhs[j][j2]+lhs[i][j2]*lhs[j][i2]);
			}
		}
	}

	ans_var = transform_sym_matrix(trans_mat_var, rhs_var);
	
	return ans_var;
}

#endif /*SUPER_BASIS3_HH_*/
