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
#ifndef MATRIX_3BY3_HH_
#define MATRIX_3BY3_HH_

#include "../utility_data_structure/nrutil_nr.hh"

inline Int4 put_complement_set3(const Int4& i, const Int4& j)
{
	assert( i != j );
	assert( 0 <= i && i < 3 );
	assert( 0 <= j && j < 3 );
	bool flag[3] = {true,true,true};
	flag[i] = false;
	flag[j] = false;
	if( flag[0] ) return 0;
	else if( flag[1] ) return 1;
	else return 2;
}


static NRMat<Int4> put_sign_mat3(const Int4& i1)
{
	NRMat<Int4> ans(3,3,0);
	ans[0][0] = 1;
	ans[1][1] = 1;
	ans[2][2] = 1;
	ans[i1][i1] = -1;
	return ans;
}


inline const NRMat<Int4>& put_sign_matrix3(const Int4& axis)
{
	static const NRMat<Int4> mat_tray[3] = { put_sign_mat3(0), put_sign_mat3(1), put_sign_mat3(2) };
	return mat_tray[axis];
}


//             -1  0  1
// Returns T =  0  1 -1
//             -1 -1 -1
static NRMat<Int4> put_perm_mat3(const Int4& i1, const Int4& i2, const Int4& i3)
{
	NRMat<Int4> ans(3,3,0);
	ans[0][i1] = 1;
	ans[1][i2] = 1;
	ans[2][i3] = 1;
	return ans;
}


inline const NRMat<Int4>& put_matrix_XYZ()
{
	static const NRMat<Int4> ans = put_perm_mat3(0,1,2);
	return ans;
}

inline const NRMat<Int4>& put_matrix_XZY()
{
	static const NRMat<Int4> ans = put_perm_mat3(0,2,1);
	return ans;
}

inline const NRMat<Int4>& put_matrix_YXZ()
{
	static const NRMat<Int4> ans = put_perm_mat3(1,0,2);
	return ans;
}

inline const NRMat<Int4>& put_matrix_YZX()
{
	static const NRMat<Int4> ans = put_perm_mat3(1,2,0);
	return ans;
}

inline const NRMat<Int4>& put_matrix_ZXY()
{
	static const NRMat<Int4> ans = put_perm_mat3(2,0,1);
	return ans;
}

inline const NRMat<Int4>& put_matrix_ZYX()
{
	static const NRMat<Int4> ans = put_perm_mat3(2,1,0);
	return ans;
}


#endif
