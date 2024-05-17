/*
 * The MIT License

   BLDConograph (Bravais lattice determination module in Conograph)

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
#ifndef _MATRIX_NBYN_HH_
#define _MATRIX_NBYN_HH_

#include"../ZAnalysisTypes.hh"
#include"../utility_data_structure/SymMat.hh"
#include"../utility_data_structure/nrutil_nr.hh"

inline void put_complement_set4(const Int4& i, const Int4& j, Int4& k, Int4& l)
{
	static const Int4 dim = 4;

	assert( i != j );
	assert( 0 <= i && i < dim );
	assert( 0 <= j && j < dim );

	bool flag_tray[dim] = {false, false, false, false};
	flag_tray[i]=true;
	flag_tray[j]=true;

	Int4 tray[2];
	for(Int4 n=0, index=0; index<2; n++)
	{
		if( !flag_tray[n] ) tray[index++] = n;
	}
	k = tray[0];
	l = tray[1];
}


// 0<=i,j<ISIZE, i!=j.
inline NRMat<Int4> put_perm_mat(const Int4& ISIZE, const Int4& i, const Int4& j)
{
	NRMat<Int4> transmat(ISIZE,ISIZE,0);
	for(Int4 k=0; k<ISIZE; k++) transmat[k][k] = 1;
	transmat[i][i] = 0;
	transmat[j][j] = 0;
	transmat[i][j] = 1;
	transmat[j][i] = 1;

	return transmat;
}
// 0<=i<4, 0<=j<4, i!=j.
inline const NRMat<Int4>& put_permutation_matrix_dim_4(const Int4& i, const Int4& j)
{
	assert( i < j );

	static const NRMat<Int4> trans_mat4[6]
	  = {
		  put_perm_mat(4, 0, 1),
		  put_perm_mat(4, 0, 2),
		  put_perm_mat(4, 0, 3),
		  put_perm_mat(4, 1, 2),
		  put_perm_mat(4, 1, 3),
		  put_perm_mat(4, 2, 3)
	  };


	return trans_mat4[ i*(5-i)/2 + j-1 ];
}


// 0<=i<4, 0<=j<4, i!=j.
static NRMat<Int4> put_reduct_mat_dim_4(const Int4& i, const Int4& j)
{
	Int4 k, l;
	put_complement_set4(i, j, k, l);
	
	NRMat<Int4> transmat(4,4,0);
	transmat[i][i] = -1;
	transmat[j][j] = 1;
	transmat[k][k] = 1; transmat[k][i] = 1;
	transmat[l][l] = 1; transmat[l][i] = 1;

	return transmat;
}


// 0<=i<4, 0<=j<4, i!=j.
inline const NRMat<Int4>& put_reduction_matrix_dim_4(const Int4& i, const Int4& j)
{
	assert( i < j );
	static const NRMat<Int4> trans_mat4[6]
	  = {
		  put_reduct_mat_dim_4(0, 1),
		  put_reduct_mat_dim_4(0, 2),
		  put_reduct_mat_dim_4(0, 3),
		  put_reduct_mat_dim_4(1, 2),
		  put_reduct_mat_dim_4(1, 3),
		  put_reduct_mat_dim_4(2, 3)
	  };


	return trans_mat4[ i*(5-i)/2 + j-1 ];
}


static NRMat<Int4> put_trans_mat43()
{
	static const Int4 isize = 3;
	NRMat<Int4> TransMat(isize+1,isize,0);
	for(Int4 k=0; k<isize; k++) TransMat[k][k] = 1;
	for(Int4 k=0; k<isize; k++) TransMat[isize][k] = -1;

	return TransMat;
}


inline const NRMat<Int4>& put_transform_matrix_row3to4()
{
	static const NRMat<Int4> mat43 = put_trans_mat43();
	return mat43;
}


inline NRMat<Int4> put_transform_matrix_row3to4(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 3 );
	
	const Int4 icol = rhs.ncols();
	
	NRMat<Int4> TransMat(4, icol);
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 k2=0; k2<icol; k2++)
		{
			TransMat[k][k2] = rhs[k][k2];
		}
	}
	for(Int4 k2=0; k2<icol; k2++)
	{
		TransMat[3][k2] = -( rhs[0][k2] + rhs[1][k2] + rhs[2][k2] );
	}

	return TransMat;
}


inline NRMat<Int4> put_transform_matrix_row4to3(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 4 );
	
	const Int4 icol = rhs.ncols();
	NRMat<Int4> TransMat(3, icol);
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 k2=0; k2<icol; k2++)
		{
			TransMat[k][k2] = rhs[k][k2];
		}
	}
	return TransMat;
}


inline NRMat<Int4> put_transform_matrix_size4to3(const NRMat<Int4>& rhs)
{
	assert( rhs.nrows() == 4 && rhs.ncols() == 4 );

	NRMat<Int4> TransMat(3, 3);
	for(Int4 k=0; k<3; k++)
	{
		for(Int4 k2=0; k2<3; k2++)
		{
			TransMat[k][k2] = rhs[k][k2] - rhs[k][3];
		}
	}
	return TransMat;
}


template<class T>
inline SymMat<T> put_sym_matrix_sizeNplus1toN(const SymMat<T>& rhs)
{
	assert( rhs.size() > 0 );

	static const Int4 isize = rhs.size() - 1;
	SymMat<T> ans(isize);
	
	for(Int4 k=0; k<isize; k++)
	{
		for(Int4 k2=k; k2<isize; k2++)
		{
			ans(k,k2) = rhs(k,k2);
		}
	}
	return ans;
}


template<class T>
inline SymMat<T> put_sym_matrix_sizeNtoNplus1(const SymMat<T>& rhs)
{
	static const Int4 isize = rhs.size();
	SymMat<T> ans(isize+1);
	
	// Copy.
	for(Int4 k=0; k<isize; k++)
	{
		for(Int4 k2=k; k2<isize; k2++)
		{
			ans(k,k2) = rhs(k,k2);
		}
	}

	ans(isize,isize) = 0.0;
	for(Int4 k=0; k<isize; k++)
	{
		ans(k, isize) = 0.0;
		for(Int4 k2=0; k2<isize; k2++)
		{
			ans(k, isize) -= rhs(k, k2);
		}
		ans(isize,isize) -= ans(k, isize);
	}

	return ans;
}

#endif
