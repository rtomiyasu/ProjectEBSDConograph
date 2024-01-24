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
#ifdef DEBUG
	#include <iostream>
#endif
#include "../utility_data_structure/index_set.hh"
#include "../utility_func/zstring.hh"
#include "../utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "../utility_lattice_reduction/put_Buerger_reduced_lattice.hh"
#include "BravaisLattice.hh"
#include "LatticeFigureOfMeritToCheckSymmetry.hh"
#include "ReducedLatticeToCheckBravais.hh"


BravaisLattice::BravaisLattice()
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		JudgeSymmetry[i] = true;
	}

	m_DoesPrudentSymSearch = false;
   	m_resol = 1.0e-8;
}


BravaisLattice::~BravaisLattice()
{
}


// Set the member variables.
void BravaisLattice::setParam(const bool& does_exhaustive_search, const Double& resolution)
{
	for(Int4 i=0; i<NUM_LS; i++)
	{
		JudgeSymmetry[i] = true;
	}

	m_DoesPrudentSymSearch = does_exhaustive_search;
	m_resol = resolution;
}




void BravaisLattice::putCentringTypes(const ReducedLatticeToCheckBravais& RLCB,
		const BravaisType& brat,
		vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result)
{
	lattice_result.clear();
	
	const map< SymMat<Double>, NRMat<Int4> >& S_red_tray = RLCB.checkCentringType(brat);
	if( S_red_tray.empty() ) return;

	// The lattice of RLCB has at least the symmetry given by eblat.
	SymMat<Double> S_super(4);
	NRMat<Int4> trans_mat(4,3);

	for(map< SymMat<Double>, NRMat<Int4> >::const_iterator it=S_red_tray.begin(); it!=S_red_tray.end(); it++)
	{
		S_super = transform_sym_matrix(it->second, it->first);
		trans_mat = identity_matrix<Int4>(4);

		// S_super = trans_mat * it->second * it->first * Transpose(trans_mat * it->second) is Delone reduced.
		if( !put_Selling_reduced_dim_3(S_super, trans_mat) )
		{
			assert( false );
		}
		moveSmallerDiagonalLeftUpper(S_super, trans_mat);

		lattice_result.push_back( LatticeFigureOfMeritToCheckSymmetry( brat, SymMat43_Double(it->first, mprod(trans_mat, it->second) ) ) );

		lattice_result.back().setTransformToOriginalLattice( Inverse3( put_transform_matrix_size4to3(trans_mat) ) );
	}
}



static SymMat43_Double putTransformMatrixFromSellingReducedToBuergerReduced(
		const SymMat<Double>& S_super_obtuse)
{
	NRMat<Int4> trans_mat(4,3);

	// S_red = trans_mat * S_super_obtuse * transpose(trans_mat).
	SymMat<Double> S_red(3);
	putBuergerReducedMatrix(S_super_obtuse, false, S_red, trans_mat);

	return SymMat43_Double( S_red, put_transform_matrix_row3to4( Inverse3( trans_mat ) ) );
}


void BravaisLattice::putLatticeCandidatesForEachBravaisTypes(const eABCaxis& abc_axis, const eRHaxis& rh_axis)
{
	for(Int4 i=1; i<NUM_LS; i++)
	{
		m_result[i].clear();
	}
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = m_result[(size_t)BravaisType::Triclinic];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_P = m_result[(size_t)BravaisType::Monoclinic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_mono_B = m_result[(size_t)BravaisType::Monoclinic_B];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_P = m_result[(size_t)BravaisType::Orthorhombic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_B = m_result[(size_t)BravaisType::Orthorhombic_C];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_I = m_result[(size_t)BravaisType::Orthorhombic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_ortho_F = m_result[(size_t)BravaisType::Orthorhombic_F];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_P = m_result[(size_t)BravaisType::Tetragonal_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tetra_I = m_result[(size_t)BravaisType::Tetragonal_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_rhom = m_result[(size_t)BravaisType::Rhombohedral];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_hex = m_result[(size_t)BravaisType::Hexagonal];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_P = m_result[(size_t)BravaisType::Cubic_P];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_I = m_result[(size_t)BravaisType::Cubic_I];
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_cubic_F = m_result[(size_t)BravaisType::Cubic_F];

	assert( lattice_result_tri.size() == 1 );

	LatticeFigureOfMeritToCheckSymmetry& latFOM = lattice_result_tri[0];
	const SymMat43_Double& S_red = latFOM.putReducedForm();

	// Calculate figures of merit as triclinic
	const ReducedLatticeToCheckBravais RLCB(abc_axis, rh_axis, m_DoesPrudentSymSearch, m_resol, S_red);

	vector<LatticeFigureOfMeritToCheckSymmetry> latFOM_tray, latFOM_tray2;

	if( JudgeSymmetry[BravaisType::Monoclinic_B] )
	{
		putCentringTypes(RLCB, BravaisType( put_monoclinic_b_type(abc_axis) ), latFOM_tray);
		lattice_result_mono_B.insert(lattice_result_mono_B.end(), latFOM_tray.begin(), latFOM_tray.end());

		if( JudgeSymmetry[BravaisType::Orthorhombic_C] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray2);
				lattice_result_ortho_B.insert(lattice_result_ortho_B.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
	}
	if( JudgeSymmetry[BravaisType::Orthorhombic_I] )
	{
		putCentringTypes(RLCB, BravaisType( pair<CenteringType::Enum, ePointGroup>(CenteringType::Inner, D2h) ), latFOM_tray);
		lattice_result_ortho_I.insert(lattice_result_ortho_I.end(), latFOM_tray.begin(), latFOM_tray.end());

		if( JudgeSymmetry[BravaisType::Tetragonal_I] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray2);
				lattice_result_tetra_I.insert(lattice_result_tetra_I.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
		if( JudgeSymmetry[BravaisType::Cubic_I] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray2);
				lattice_result_cubic_I.insert(lattice_result_cubic_I.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
	}
	if( JudgeSymmetry[BravaisType::Orthorhombic_F] )
	{
		putCentringTypes(RLCB, BravaisType( pair<CenteringType::Enum, ePointGroup>(CenteringType::Face, D2h) ), latFOM_tray);
		lattice_result_ortho_F.insert(lattice_result_ortho_F.end(), latFOM_tray.begin(), latFOM_tray.end());

		if( JudgeSymmetry[BravaisType::Cubic_F] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray2);
				lattice_result_cubic_F.insert(lattice_result_cubic_F.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
	}
	if( JudgeSymmetry[BravaisType::Rhombohedral] )
	{
		putCentringTypes(RLCB, BravaisType( put_rhombohedral_type(rh_axis) ), latFOM_tray);
		lattice_result_rhom.insert(lattice_result_rhom.end(), latFOM_tray.begin(), latFOM_tray.end());
	}
	if( JudgeSymmetry[BravaisType::Monoclinic_P] )
	{
		latFOM.putLatticesOfHigherSymmetry(put_monoclinic_p_type(abc_axis), m_resol, latFOM_tray);
		lattice_result_mono_P.insert(lattice_result_mono_P.end(), latFOM_tray.begin(), latFOM_tray.end());

		if( JudgeSymmetry[BravaisType::Hexagonal] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(D6h, m_resol, latFOM_tray2);
				lattice_result_hex.insert(lattice_result_hex.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
	}
	if( JudgeSymmetry[BravaisType::Orthorhombic_P] )
	{
		latFOM.putLatticesOfHigherSymmetry(D2h, m_resol, latFOM_tray);
		lattice_result_ortho_P.insert(lattice_result_ortho_P.end(), latFOM_tray.begin(), latFOM_tray.end());

		if( JudgeSymmetry[BravaisType::Tetragonal_P] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(D4h_Z, m_resol, latFOM_tray2);
				lattice_result_tetra_P.insert(lattice_result_tetra_P.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
		if( JudgeSymmetry[BravaisType::Cubic_P] )
		{
			for(vector<LatticeFigureOfMeritToCheckSymmetry>::const_iterator it=latFOM_tray.begin(); it!=latFOM_tray.end(); it++)
			{
				it->putLatticesOfHigherSymmetry(Oh, m_resol, latFOM_tray2);
				lattice_result_cubic_P.insert(lattice_result_cubic_P.end(), latFOM_tray2.begin(), latFOM_tray2.end());
			}
		}
	}

//#ifdef DEBUG
//	for(Int4 i=1; i<NUM_LS; i++)
//	{
//		if( !JudgeSymmetry[i] ) continue;
//
//cout << "(" + num2str( i+1 ) + ") The number of candidates for " + BravaisType::ToString(BravaisType::Enum(i), abc_axis)
//			+ " : " + num2str<size_t>( m_result[i].size() ) + "\n";
//	}
//cout << "\n";
//#endif
}


void BravaisLattice::execute(const SymMat<Double>& S_super,
					const eABCaxis& abc_axis,
					const eRHaxis& rh_axis)
{
	vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result_tri = m_result[(size_t)BravaisType::Triclinic];

	const SymMat43_Double S_red = putTransformMatrixFromSellingReducedToBuergerReduced(S_super);
	LatticeFigureOfMeritToCheckSymmetry latFOM(BravaisType( pair<CenteringType::Enum, ePointGroup>(CenteringType::Prim, Ci) ), S_red);

	lattice_result_tri.clear();
	lattice_result_tri.push_back( latFOM );

	putLatticeCandidatesForEachBravaisTypes(abc_axis, rh_axis);
}


const LatticeFigureOfMeritToCheckSymmetry& BravaisLattice::putBravaisLattice() const
{
	for(Int4 i=NUM_LS-1; i>0; i--)
	{
		if( !m_result[i].empty() )
		{
			assert( m_result[i].size() == 1 );
			return *m_result[i].begin();
		}
	}
	assert( m_result[0].size() == 1 );
	return *m_result[0].begin();
}
