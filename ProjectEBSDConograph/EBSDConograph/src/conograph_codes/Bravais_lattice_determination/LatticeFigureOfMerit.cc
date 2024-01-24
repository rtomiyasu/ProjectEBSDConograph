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
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"

#include "ReducedLatticeToCheckBravais.hh"
#include "LatticeFigureOfMerit.hh"

const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_face = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToFace() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_body = put_transform_matrix_row3to4( BravaisType::putTransformMatrixFromBodyToPrimitive() );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_rhomhex = put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToRhomHex() ) );
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_base[3] =
		{
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseA_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseB_Axis) ) ),
				put_transform_matrix_row3to4( transpose( BravaisType::putTransformMatrixFromPrimitiveToBase(BaseC_Axis) ) )
		};
const NRMat<Int4> LatticeFigureOfMerit::m_tmat_prim_to_prim = put_transform_matrix_row3to4();

LatticeFigureOfMerit::LatticeFigureOfMerit()
	: m_S_red(3)
{
}


LatticeFigureOfMerit::LatticeFigureOfMerit(const BravaisType& brat, const SymMat43_Double& S) : m_S_red(3)
{
	m_brat = brat;

	NRMat<Int4> trans_mat;
	setInverseOfBuergerReducedForm(trans_mat, S);	// Set m_S_red from S.
}

// Set m_S_red from S_red.
// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) equals the original S.
LatticeFigureOfMerit::LatticeFigureOfMerit(const BravaisType& brat, const SymMat43_Double& S, NRMat<Int4>& trans_mat) : m_S_red(3)
{
	m_brat = brat;

	setInverseOfBuergerReducedForm(trans_mat, S);	// Set m_S_red from S.
}

#ifdef DEBUG
static bool checkInitialLatticeParameters(
		const BravaisType& brat,
		const SymMat<Double>& S_red)
{
	const SymMat<Double> inv_S_red( Inverse3(S_red) );

	if( brat.enumLatticeSymmetry() == C2h_Y && brat.enumCenteringType() == CenteringType::Prim )
	{
		assert( inv_S_red(0,2) <= 0.0 &&
				inv_S_red(0,0) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumLatticeSymmetry() == C2h_Z && brat.enumCenteringType() == CenteringType::Prim )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLatticeSymmetry() == C2h_X && brat.enumCenteringType() == CenteringType::Prim )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumLatticeSymmetry() == C2h_Y && brat.enumCenteringType() == CenteringType::BaseC )
	{
		assert( inv_S_red(0,2) <= 0.0
				&& fabs( inv_S_red(0,2) ) * 0.9999 < inv_S_red(2,2)
				&& fabs( inv_S_red(0,2) ) * 1.9999 < inv_S_red(0,0) );
	}
	else if( brat.enumLatticeSymmetry() == C2h_Z && brat.enumCenteringType() == CenteringType::BaseA )
	{
		assert( inv_S_red(0,1) <= 0.0
				&& fabs( inv_S_red(0,1) ) * 0.9999 < inv_S_red(0,0)
				&& fabs( inv_S_red(0,1) ) * 1.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLatticeSymmetry() == C2h_X && brat.enumCenteringType() == CenteringType::BaseB )
	{
		assert( inv_S_red(1,2) <= 0.0
				&& fabs( inv_S_red(1,2) ) * 0.9999 < inv_S_red(1,1)
				&& fabs( inv_S_red(1,2) ) * 1.9999 < inv_S_red(2,2) );
	}
	else if( brat.enumBravaisType() == BravaisType::Orthorhombic_C )
	{
		assert( brat.enumCenteringType() == CenteringType::BaseC );
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1) );
	}
	else if( brat.enumLatticeSymmetry() == D2h && brat.enumCenteringType() == CenteringType::Prim )
	{
		assert( inv_S_red(0,0) * 0.9999 < inv_S_red(1,1)
				&& inv_S_red(1,1) * 0.9999 < inv_S_red(2,2) );
	}
	return true;
}
#endif

static void putTransformMatrixToBuergerReduced(
		const SymMat<Double>& S, NRMat<Int4>& trans_mat)
{
	assert( S.size() == 3 );

	SymMat<Double> S_super_obtuse(4);
	put_Selling_reduced_dim_3(S, S_super_obtuse, trans_mat);
	moveSmallerDiagonalLeftUpper(S_super_obtuse, trans_mat);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	NRMat<Int4> trans_mat2;
	putBuergerReducedMatrix(S_super_obtuse, S_red, trans_mat2);
	trans_mat = mprod( trans_mat2, put_transform_matrix_row4to3(trans_mat) );
}


// Set m_S_red from S.
// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) equals the original S.
void LatticeFigureOfMerit::setInverseOfBuergerReducedForm(NRMat<Int4>& trans_mat, const SymMat43_Double& S_optimized)
{
	if( m_brat.enumBravaisType() == BravaisType::Triclinic )
	{
		// trans_mat * Inverse(S_optimized.first) * transpose(trans_mat) is Buerger-reduced
		// <=> Inverse of transpose(Inverse(trans_mat)) * S_optimized.first * Inverse(trans_mat) is Buerger-reduced.
		putTransformMatrixToBuergerReduced(Inverse3(S_optimized.first), trans_mat);
		transpose_square_matrix(trans_mat);
		m_S_red = transform_sym_matrix(Inverse3(trans_mat), S_optimized.first);
	}
	else
	{
		m_S_red = S_optimized.first;
		trans_mat = identity_matrix<Int4>(3);
		if( m_brat.enumBravaisType() == BravaisType::Monoclinic_P )
		{
			if( m_brat.enumLatticeSymmetry() == C2h_X )
			{
				putBuergerReducedMonoclinicP(1, 2, m_S_red, trans_mat);
			}
			else if( m_brat.enumLatticeSymmetry() == C2h_Y )
			{
				putBuergerReducedMonoclinicP(0, 2, m_S_red, trans_mat);
			}
			else //if( m_brat.enumPointGroup() == C2h_Z )
			{
				putBuergerReducedMonoclinicP(0, 1, m_S_red, trans_mat);
			}
		}
		else if( m_brat.enumBravaisType() == BravaisType::Monoclinic_B )
		{
			putBuergerReducedMonoclinicB(m_brat, m_S_red, trans_mat);
		}
		else if( m_brat.enumLatticeSymmetry() == D2h )
		{
			putBuergerReducedOrthorhombic(m_brat.enumCenteringType(), m_S_red, trans_mat);
		}
	}

	assert( checkInitialLatticeParameters(m_brat, m_S_red) );
}




inline bool checkIfFirstEntryIsPositive(const VecDat3<Int4>& rhs)
{
	for(Int4 i=0; i<3; i++)
	{
		if( rhs[i] == 0 ) continue;
		if( rhs[i] > 0 ) return true;
		else return false;
	}
	return false;
}
