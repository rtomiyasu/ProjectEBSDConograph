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
#include "../bravais_type/FracMat.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_func/zmath.hh"
#include "check_equiv.hh"
#include "ReducedLatticeToCheckBravais.hh"


static void put_transform_matrix_from_sell_to_neighbor_base(vector< NRMat<Int4> >& arg,
		const bool& does_prudent_search)
{
	// { hB * g : g^{-1} is one of the following matrices } are the matrices in Table 6, Oishi-Tomiyasu, Acta. Cryst. (2012).
	static const Int4 ISIZE = 69;
	static const Int4 mat_tray[ISIZE][3][3]
	= {
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { 0, 1, 1 } },
			{ { 1, 0, 0 },
			  { -1, 1, 0 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { -1, 1, 0 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { -1, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { -1, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { -1, 0, -1 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { -1, -1, -1 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, 1, 0 },
			  { 1, 0, 0 },
			  { 0, 0, 1 } },
			{ { -1, 1, 0 },
			  { 1, 0, 0 },
			  { 0, -1, -1 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { -1, 1, 0 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { 0, -1, -1 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { -1, -1, -1 } },
			{ { -1, 1, 0 },
			  { 0, 0, 1 },
			  { 1, 0, 0 } },
			{ { 0, 0, 1 },
			  { -1, 1, 0 },
			  { 1, 0, 0 } },
			{ { 0, 0, 1 },
			  { -1, -1, -1 },
			  { 1, 0, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { -1, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { -1, 0, 1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 1, -1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { -1, 1, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { 0, 0, -1 } },
			{ { 1, 0, 0 },
			  { 0, 0, -1 },
			  { 0, 1, 1 } },
			{ { 1, 0, 0 },
			  { 0, 0, -1 },
			  { -1, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { -1, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, -1, 0 },
			  { 0, 0, -1 } },
			{ { 1, 0, 0 },
			  { -1, -1, 0 },
			  { 0, 1, 1 } },
			{ { 1, 0, 0 },
			  { 0, 0, 1 },
			  { -1, 1, -1 } },
			{ { 1, 0, 0 },
			  { -1, 1, -1 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { -1, 0, 1 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { -1, 0, 1 } },
			{ { -1, 0, 1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { 0, -1, -1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { -1, 1, -1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { 0, 1, 1 },
			  { 1, 0, 0 },
			  { 0, 0, -1 } },
			{ { 0, 0, -1 },
			  { 1, 0, 0 },
			  { 0, 1, 1 } },
			{ { 0, 0, -1 },
			  { 1, 0, 0 },
			  { -1, -1, 0 } },
			{ { 0, 1, 1 },
			  { 1, 0, 0 },
			  { -1, -1, 0 } },
			{ { -1, -1, 0 },
			  { 1, 0, 0 },
			  { 0, 0, -1 } },
			{ { -1, -1, 0 },
			  { 1, 0, 0 },
			  { 0, 1, 1 } },
			{ { 0, 0, 1 },
			  { 1, 0, 0 },
			  { -1, 1, -1 } },
			{ { -1, 1, -1 },
			  { 1, 0, 0 },
			  { 0, 0, 1 } },
			{ { -1, 0, 1 },
			  { 1, 0, 0 },
			  { 0, -1, -1 } },
			{ { 0, -1, -1 },
			  { 1, 0, 0 },
			  { -1, 0, 1 } },
			{ { 0, 1, 1 },
			  { 0, 0, -1 },
			  { 1, 0, 0 } },
			{ { 0, 0, -1 },
			  { 0, 1, 1 },
			  { 1, 0, 0 } },
			{ { 0, 0, -1 },
			  { -1, -1, 0 },
			  { 1, 0, 0 } },
			{ { 0, 1, 1 },
			  { -1, -1, 0 },
			  { 1, 0, 0 } },
			{ { -1, -1, 0 },
			  { 0, 0, -1 },
			  { 1, 0, 0 } },
			{ { -1, -1, 0 },
			  { 0, 1, 1 },
			  { 1, 0, 0 } },
			{ { 0, 0, 1 },
			  { -1, 1, -1 },
			  { 1, 0, 0 } },
			{ { -1, 1, -1 },
			  { 0, 0, 1 },
			  { 1, 0, 0 } },
			{ { -1, 0, 1 },
			  { 0, -1, -1 },
			  { 1, 0, 0 } },
			{ { 0, -1, -1 },
			  { -1, 0, 1 },
			  { 1, 0, 0 } },
			{ { 1, 1, 0 },
			  { 0, 0, 1 },
			  { 0, -1, -1 } },
			{ { 1, 1, 0 },
			  { 0, -1, -1 },
			  { 0, 0, 1 } },
			{ { 0, 0, 1 },
			  { 1, 1, 0 },
			  { 0, -1, -1 } },
			{ { 0, -1, -1 },
			  { 1, 1, 0 },
			  { 0, 0, 1 } },
			{ { 0, 0, 1 },
			  { 0, -1, -1 },
			  { 1, 1, 0 } },
			{ { 0, -1, -1 },
			  { 0, 0, 1 },
			  { 1, 1, 0 } }
	};

	const Int4 ISIZE2 = (does_prudent_search?ISIZE:21);
	arg.clear();
	arg.resize(ISIZE2, NRMat<Int4>(3,3));
	for(Int4 i=0; i<ISIZE2; i++)
	{
		NRMat<Int4>& arg_ref = arg[i];
		const Int4 (*mat)[3] = mat_tray[i];
		for(Int4 i2=0; i2<3; i2++)
		{
			for(Int4 j2=0; j2<3; j2++)
			{
				arg_ref[i2][j2] = mat[i2][j2];
			}
		}
	}
}


static void put_transform_matrix_from_sell_to_neighbor_rhom(vector< NRMat<Int4> >& arg,
		const bool& does_prudent_search)
{

	// The inverses of the matrices in Table 6, Oishi-Tomiyasu, Acta. Cryst. (2012).
	static const Int4 ISIZE = 64;
	static const Int4 mat_tray[ISIZE][3][3]
	= {
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { -1, -1, -1 } },
			{ { 1, 0, 0 },
			  { -1, -1, -1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 0, 1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, -1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { -1, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { 0, 1, -1 } },
			{ { 1, 0, 0 },
			  { -1, 0, 1 },
			  { 0, 1, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, -1 },
			  { -1, 0, 1 } },
			{ { -1, -1, -1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, 0, 1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { 0, 1, -1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { -1, 0, 1 },
			  { 1, 0, 0 },
			  { 0, 1, -1 } },
			{ { 0, 1, -1 },
			  { 1, 0, 0 },
			  { -1, 0, 1 } },
			{ { -1, 0, 1 },
			  { 0, 1, -1 },
			  { 1, 0, 0 } },
			{ { 0, 1, -1 },
			  { -1, 0, 1 },
			  { 1, 0, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 0, -1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { 0, 1, 1 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { -1, 0, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 1 },
			  { -1, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, -1, 0 },
			  { 0, 1, 1 } },
			{ { 0, 1, 1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { -1, 0, -1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { 0, 1, 1 },
			  { 1, 0, 0 },
			  { -1, -1, 0 } },
			{ { -1, -1, 0 },
			  { 1, 0, 0 },
			  { 0, 1, 1 } },
			{ { 0, 1, 1 },
			  { -1, -1, 0 },
			  { 1, 0, 0 } },
			{ { -1, -1, 0 },
			  { 0, 1, 1 },
			  { 1, 0, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, 0, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { -1, 0, 1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { -1, -1, 1 } },
			{ { 1, 0, 0 },
			  { 0, 0, -1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 0, 1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { -1, -1, 1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 1, -1 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { -1, 1, -1 } },
			{ { 1, 0, 0 },
			  { -1, 1, 0 },
			  { 0, -1, -1 } },
			{ { 1, 0, 0 },
			  { 0, -1, -1 },
			  { -1, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { 0, -1, 1 } },
			{ { 1, 0, 0 },
			  { 0, 1, 0 },
			  { -1, 0, -1 } },
			{ { 1, 0, 0 },
			  { 0, 0, -1 },
			  { 0, -1, 0 } },
			{ { -1, 0, 1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { 0, -1, -1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, -1, 1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, 1, -1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 0 },
			  { -1, 1, 1 } },
			{ { -1, 1, 0 },
			  { 1, 0, 0 },
			  { 0, -1, -1 } },
			{ { 0, -1, -1 },
			  { 1, 0, 0 },
			  { -1, 1, 0 } },
			{ { 1, 0, 0 },
			  { 0, -1, 1 },
			  { 0, 1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 0, -1 },
			  { 0, 1, 0 } },
			{ { 0, -1, 1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, 0, -1 },
			  { 1, 0, 0 },
			  { 0, 1, 0 } },
			{ { -1, 1, 1 },
			  { 1, 0, 0 },
			  { 0, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, 1, 1 },
			  { 0, -1, 0 } },
			{ { -1, 1, 0 },
			  { 0, -1, -1 },
			  { 1, 0, 0 } },
			{ { 0, -1, -1 },
			  { -1, 1, 0 },
			  { 1, 0, 0 } },
			{ { -1, -1, 0 },
			  { 0, 1, -1 },
			  { 1, 0, 0 } },
			{ { 0, 1, -1 },
			  { -1, -1, 0 },
			  { 1, 0, 0 } },
			{ { -1, -1, 0 },
			  { 1, 0, 0 },
			  { 0, 1, -1 } },
			{ { 0, 1, -1 },
			  { 1, 0, 0 },
			  { -1, -1, 0 } },
			{ { 1, 0, 0 },
			  { -1, -1, 0 },
			  { 0, 1, -1 } },
			{ { 1, 0, 0 },
			  { 0, 1, -1 },
			  { -1, -1, 0 } }
	};

	const Int4 ISIZE2 = (does_prudent_search?ISIZE:16);
	arg.clear();
	arg.resize(ISIZE2, NRMat<Int4>(3,3));
	for(Int4 i=0; i<ISIZE2; i++)
	{
		NRMat<Int4>& arg_ref = arg[i];
		const Int4 (*mat)[3] = mat_tray[i];
		for(Int4 i2=0; i2<3; i2++)
		{
			for(Int4 j2=0; j2<3; j2++)
			{
				arg_ref[i2][j2] = mat[i2][j2];
			}
		}
	}
}




// The second variable is the inverse matrix of the first variable.
static vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > put_Transform_Matrix_base()
{
	static const NRMat<Int4> tmat_prim_to_Acell1 = transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseA_Axis) );
	static const NRMat<Int4> tmat_prim_to_Bcell1 = transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseB_Axis) );
	static const NRMat<Int4> tmat_prim_to_Ccell1 = transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseC_Axis) );

	vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > S_min_to_sell(6);
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minA_to_sell_qck = S_min_to_sell[(size_t)BaseA_Axis*2];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minA_to_sell_prd = S_min_to_sell[(size_t)BaseA_Axis*2+1];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minB_to_sell_qck = S_min_to_sell[(size_t)BaseB_Axis*2];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minB_to_sell_prd = S_min_to_sell[(size_t)BaseB_Axis*2+1];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minC_to_sell_qck = S_min_to_sell[(size_t)BaseC_Axis*2];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_minC_to_sell_prd = S_min_to_sell[(size_t)BaseC_Axis*2+1];

	vector< NRMat<Int4> > mat_tray;
	NRMat<Int4> mat(3,3);
	put_transform_matrix_from_sell_to_neighbor_base(mat_tray, false);

	for(vector< NRMat<Int4> >::const_iterator it=mat_tray.begin(); it!=mat_tray.end(); it++)
	{
		mat = mprod(*it, tmat_prim_to_Acell1);
		S_minA_to_sell_qck.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

		mat = mprod(*it, tmat_prim_to_Bcell1);
		S_minB_to_sell_qck.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

		mat = mprod(*it, tmat_prim_to_Ccell1);
		S_minC_to_sell_qck.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
	}

	put_transform_matrix_from_sell_to_neighbor_base(mat_tray, true);

	for(vector< NRMat<Int4> >::const_iterator it=mat_tray.begin(); it!=mat_tray.end(); it++)
	{
		mat = mprod(*it, tmat_prim_to_Acell1);
		S_minA_to_sell_prd.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
		mat = mprod(*it, tmat_prim_to_Bcell1);
		S_minB_to_sell_prd.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
		mat = mprod(*it, tmat_prim_to_Ccell1);
		S_minC_to_sell_prd.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
	}

	return S_min_to_sell;
}


// The second variable is the inverse matrix of the first variable.
static vector< pair< NRMat<Int4>, FracMat<Int4> > > put_Transform_Matrix_face()
{
	// \tr{hF^{-1}}(2 3), hF is the matrix in Table 4, Oishi-Tomiyasu, Acta Cryst. (2012).
	static const NRMat<Int4> tmat_prim_to_face = transpose( BravaisType::putTransformMatrixFromPrimitiveToFace() );

	vector< pair< NRMat<Int4>, FracMat<Int4> > > S_min_to_sell;

	// < \tr{hF^{-1}}(2 3), (2 3) \tr{hF} >
	S_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( tmat_prim_to_face, FInverse3<Int4>( tmat_prim_to_face ) ) );

	// < (1 3 2) \tr{hF^{-1}}(2 3), (2 3) \tr{hF} (1 2 3) >
	NRMat<Int4> mat = mprod(put_matrix_ZXY(), tmat_prim_to_face);
	S_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

	// < (1 2 3) \tr{hF^{-1}}(2 3), (2 3) \tr{hF} (1 3 2) >
	mat = mprod(put_matrix_YZX(), tmat_prim_to_face);
	S_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

	return S_min_to_sell;
}


// The second variable is the inverse matrix of the first variable.
static vector< pair< NRMat<Int4>, FracMat<Int4> > > put_Transform_Matrix_body()
{
	// \tr{hI^{-1}}(2 3), hI is the matrix in Table 5, Oishi-Tomiyasu, Acta Cryst. (2012).
	static const NRMat<Int4> tmat_prim_to_body = BravaisType::putTransformMatrixFromBodyToPrimitive();

	vector< pair< NRMat<Int4>, FracMat<Int4> > > InvS_min_to_sell;

	// < \tr{hI^{-1}}(2 3), (2 3) \tr{hI} >
	InvS_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( tmat_prim_to_body, FInverse3<Int4>( tmat_prim_to_body ) ) );

	// < (1 3 2) \tr{hI^{-1}}(2 3), (2 3) \tr{hI} (1 2 3) >
	NRMat<Int4> mat = mprod(put_matrix_ZXY(), tmat_prim_to_body);
	InvS_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

	// < (1 2 3) \tr{hI^{-1}}(2 3), (2 3) \tr{hI} (1 3 2) >
	mat = mprod(put_matrix_YZX(), tmat_prim_to_body);
	InvS_min_to_sell.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );

	return InvS_min_to_sell;
}


// The second variable is the inverse matrix of the first variable.
static vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > put_Transform_Matrix_rhom()
{
	// Definition of \tr{hR}.
	static const NRMat<Int4> tmat_prim_to_rhomhex = transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() );

	vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > S_min_to_sell(4);
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_min_rho_to_sell_qck = S_min_to_sell[(size_t)Rho_Axis*2];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_min_rho_to_sell_prd = S_min_to_sell[(size_t)Rho_Axis*2+1];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_min_hex_to_sell_qck = S_min_to_sell[(size_t)Hex_Axis*2];
	vector< pair< NRMat<Int4>, FracMat<Int4> > >& S_min_hex_to_sell_prd = S_min_to_sell[(size_t)Hex_Axis*2+1];

	vector< NRMat<Int4> > mat_tray;
	NRMat<Int4> mat(3,3);
	put_transform_matrix_from_sell_to_neighbor_rhom(mat_tray, false);

	for(vector< NRMat<Int4> >::const_iterator it=mat_tray.begin(); it!=mat_tray.end(); it++)
	{
		// < g^{-1}, g >
		S_min_rho_to_sell_qck.push_back( pair< NRMat<Int4>, FracMat<Int4> >( *it, FInverse3<Int4>( *it ) ) );

		// < g^{-1}*\tr{hR}, \tr{hR}^{-1}*g >
		mat = mprod(*it, tmat_prim_to_rhomhex);
		S_min_hex_to_sell_qck.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
	}

	put_transform_matrix_from_sell_to_neighbor_rhom(mat_tray, true);

	for(vector< NRMat<Int4> >::const_iterator it=mat_tray.begin(); it!=mat_tray.end(); it++)
	{
		// < g^{-1}, g >
		S_min_rho_to_sell_prd.push_back( pair< NRMat<Int4>, FracMat<Int4> >( *it, FInverse3<Int4>( *it ) ) );

		// < g^{-1}*\tr{hR}, \tr{hR}^{-1}*g >
		mat = mprod(*it, tmat_prim_to_rhomhex);
		S_min_hex_to_sell_prd.push_back( pair< NRMat<Int4>, FracMat<Int4> >( mat, FInverse3<Int4>( mat ) ) );
	}

	return S_min_to_sell;
}


const vector< pair< NRMat<Int4>, FracMat<Int4> > > ReducedLatticeToCheckBravais::m_trans_mat_red_F = put_Transform_Matrix_face();
const vector< pair< NRMat<Int4>, FracMat<Int4> > > ReducedLatticeToCheckBravais::m_trans_mat_red_I = put_Transform_Matrix_body();
const vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > ReducedLatticeToCheckBravais::m_trans_mat_red_rhom = put_Transform_Matrix_rhom();
const vector< vector< pair< NRMat<Int4>, FracMat<Int4> > > > ReducedLatticeToCheckBravais::m_trans_mat_red_base = put_Transform_Matrix_base();

ReducedLatticeToCheckBravais::ReducedLatticeToCheckBravais(
		const eABCaxis& axis1,
		const eRHaxis& axis2,
		const bool& does_prudent_sym_search,
		const Double& resol, const SymMat43_Double & S_red)
 : m_monoclinic_b_type(put_monoclinic_b_type(axis1)),
   m_rhombohedral_type(put_rhombohedral_type(axis2)),
   m_S_super_obtuse( transform_sym_matrix(S_red.second, S_red.first) )
{
	put_S_Buerger_reduced_IF(resol, m_S_super_obtuse, m_S_red_body, false);
	put_S_Buerger_reduced_base(m_monoclinic_b_type, does_prudent_sym_search, resol, m_S_super_obtuse, m_S_red_base);
	put_S_Buerger_reduced_rhom(m_rhombohedral_type, does_prudent_sym_search, resol, m_S_super_obtuse, m_S_red_rhom);

	const SymMat<Double> S_super_obtuse3( put_sym_matrix_sizeNplus1toN(m_S_super_obtuse) );
	const SymMat<Double> inv_S( Inverse3( S_super_obtuse3 ) );

	// Calculate the inverse of m_S_red.
	SymMat<Double> inv_S_super_obtuse(4);
	NRMat<Int4> tmat_inv_S_super_obtuse(4,3);

	// inv_S_super_obtuse = transpose( tmat_inv_S_super_obtuse)  * inverse(S_super_obtuse3) * tmat_inv_S_super_obtuse.
	put_Selling_reduced_dim_3(inv_S, inv_S_super_obtuse, tmat_inv_S_super_obtuse);
	moveSmallerDiagonalLeftUpper(inv_S_super_obtuse, tmat_inv_S_super_obtuse);
	tmat_inv_S_super_obtuse = put_transform_matrix_row4to3(tmat_inv_S_super_obtuse);
	transpose_square_matrix(tmat_inv_S_super_obtuse);

	const SymMat<Double> S_inv_super_obtuse
			= put_sym_matrix_sizeNtoNplus1( transform_sym_matrix( Inverse3(tmat_inv_S_super_obtuse), S_super_obtuse3 ) );

	put_S_Buerger_reduced_IF(resol, S_inv_super_obtuse, m_S_red_face, true);

	for(map< SymMat<Double>, NRMat<Int4> >::iterator it=m_S_red_face.begin(); it!=m_S_red_face.end(); it++)
	{
		it->second = put_transform_matrix_row3to4( mprod(tmat_inv_S_super_obtuse, put_transform_matrix_row4to3(it->second) ) );
	}
}


ReducedLatticeToCheckBravais::~ReducedLatticeToCheckBravais()
{
}


// On input, inv_flag = false indicates that S_super_obtuse_equiv is Selling-reduced,
// and inv_flag = true indicates that Inverse(S_super_obtuse_equiv) is Selling-reduced.
// In the former case, on output, S_red_body are symmetric matrices having a body-centered and Buerger-reduced inverse.
// In the latter case, on output, S_red_IF are symmetric matrices having a face-centered and Buerger-reduced inverse.
void ReducedLatticeToCheckBravais::put_S_Buerger_reduced_IF(
		const Double& resol, const SymMat<Double>& S_super_obtuse,
		map< SymMat<Double>, NRMat<Int4> >& S_red_IF,
		const bool& inv_flag)
{
	const vector< pair< NRMat<Int4>, FracMat<Int4> > >& tmat_red_IF = (inv_flag?m_trans_mat_red_F:m_trans_mat_red_I);
	S_red_IF.clear();

	NRMat<Int4> tmat;
	SymMat<Double> S2_red0(3), S2_red(3);

	for(vector< pair< NRMat<Int4>, FracMat<Int4> > >::const_iterator it=tmat_red_IF.begin(); it!=tmat_red_IF.end(); it++)
	{
		const FracMat<Int4>& inv_mat = it->second;
		S2_red0 = transform_sym_matrix(inv_mat.mat, put_sym_matrix_sizeNplus1toN(S_super_obtuse) ) / (inv_mat.denom*inv_mat.denom);
		S2_red = S2_red0;

		cal_average_crystal_system(D2h, S2_red);
		if( !check_equiv_m(S2_red0, S2_red, resol) ) continue;

		tmat = identity_matrix<Int4>(3);
		moveLargerDiagonalLeftUpper(S2_red, tmat);
		tmat = mprod( put_transform_matrix_row3to4(it->first), transpose(tmat) );	// inverse(tmat) = transpose(tmat).

		S_red_IF.insert( SymMat43_Double(S2_red, tmat) );
	}
}

void ReducedLatticeToCheckBravais::put_S_Buerger_reduced_rhom(
		const BravaisType& rhombohedral_type,
		const bool& does_prudent_sym_search,
		const Double& resol, const SymMat<Double>& S_super_obtuse,
		map< SymMat<Double>, NRMat<Int4> >& S_red_rhomhex)
{
	const vector< pair< NRMat<Int4>, FracMat<Int4> > >& tmat_red_rhom = m_trans_mat_red_rhom[(size_t)rhombohedral_type.enumRHaxis()*2+(does_prudent_sym_search?1:0)];
	S_red_rhomhex.clear();

	NRMat<Int4> tmat;
	SymMat<Double> S2_red0(3), S2_red(3);

	for(vector< pair< NRMat<Int4>, FracMat<Int4> > >::const_iterator it=tmat_red_rhom.begin(); it!=tmat_red_rhom.end(); it++)
	{
		const FracMat<Int4>& inv_mat = it->second;
		S2_red0 = transform_sym_matrix(inv_mat.mat, put_sym_matrix_sizeNplus1toN(S_super_obtuse) ) / (inv_mat.denom*inv_mat.denom);
		S2_red = S2_red0;

		cal_average_crystal_system(rhombohedral_type.enumLatticeSymmetry(), S2_red);
		if( !check_equiv_m(S2_red0, S2_red, resol) ) continue;

		tmat = put_transform_matrix_row3to4(it->first);

		S_red_rhomhex.insert( SymMat43_Double(S2_red, tmat) );
	}
}


// On input, S_red is Buerger-reduced and S_super_obtuse_equiv is Selling-reduced.
// On output, S_red_base are symmetric matrices having a base-centered and Buerger-reduced inverse.
void ReducedLatticeToCheckBravais::put_S_Buerger_reduced_base(
		const BravaisType& monoclinic_b_type,
		const bool& does_prudent_sym_search,
		const Double& resol, const SymMat<Double>& S_super_obtuse,
		map< SymMat<Double>, NRMat<Int4> >& S_red_base)
{
	const eBASEaxis ibase_axis = monoclinic_b_type.enumBASEaxis();
	const vector< pair< NRMat<Int4>, FracMat<Int4> > >& tmat_red_base = m_trans_mat_red_base[(size_t)ibase_axis*2+(does_prudent_sym_search?1:0)];

	S_red_base.clear();

	NRMat<Int4> tmat;
	SymMat<Double> S2_red0(3), S2_red(3);

	for(vector< pair< NRMat<Int4>, FracMat<Int4> > >::const_iterator it=tmat_red_base.begin(); it!=tmat_red_base.end(); it++)
	{
		const FracMat<Int4>& inv_mat = it->second;
		S2_red0 = transform_sym_matrix(inv_mat.mat, put_sym_matrix_sizeNplus1toN(S_super_obtuse) ) / (inv_mat.denom*inv_mat.denom);
		S2_red = S2_red0;

		cal_average_crystal_system(monoclinic_b_type.enumLatticeSymmetry(), S2_red);

		if( !check_equiv_m(S2_red0, S2_red, resol) ) continue;

		tmat = put_transform_matrix_row3to4(it->first);
		putBuergerReducedMonoclinicB(monoclinic_b_type, S2_red, tmat);

		S_red_base.insert( SymMat43_Double(S2_red, tmat) );
	}
}


const map< SymMat<Double>, NRMat<Int4> >& ReducedLatticeToCheckBravais::checkCentringType(const BravaisType& brat) const
{
	if( brat == m_monoclinic_b_type )
	{
		return m_S_red_base;
	}
	else if( brat.enumCenteringType() == CenteringType::Face )
	{
		return m_S_red_face;
	}
	else if( brat.enumCenteringType() == CenteringType::Inner )
	{
		return m_S_red_body;
	}
	else if( brat == m_rhombohedral_type )
	{
		return m_S_red_rhom;
	}
	else
	{
		assert(false);
		return m_S_red_body;
	}
}

