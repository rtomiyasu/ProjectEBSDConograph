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
#ifndef LATTICEFIGUREOFMERIT_HH_
#define LATTICEFIGUREOFMERIT_HH_

#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_data_structure/VecDat3.hh"
#include "../utility_data_structure/SymMat.hh"
#include "../utility_func/lattice_constant.hh"
#include "../point_group/enumPointGroup.hh"
#include "../bravais_type/BravaisType.hh"
#include "SymMatNplus1N.hh"

class LatticeFigureOfMeritToCheckSymmetry;

// Class for outputting information about a lattice in index file.
class LatticeFigureOfMerit
{
private:
	static const NRMat<Int4> m_tmat_prim_to_face;
	static const NRMat<Int4> m_tmat_prim_to_body;
	static const NRMat<Int4> m_tmat_prim_to_rhomhex;
	static const NRMat<Int4> m_tmat_prim_to_base[3];
	static const NRMat<Int4> m_tmat_prim_to_prim;

	BravaisType m_brat;

	// The inverse of m_S_red is Buerger-reduced.
	SymMat<Double> m_S_red;

	// Set m_S_red from S.
	// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) equals the original S.
	void setInverseOfBuergerReducedForm(NRMat<Int4>& trans_mat, const SymMat43_Double& S);

public:
	LatticeFigureOfMerit();
	LatticeFigureOfMerit(const BravaisType& brat, const SymMat43_Double& S_red);

	// Set m_S_red from S_red.
	// On output, trans_mat gives the matrix such that trans_mat * m_S_red * transpose(trans_mat) equals the original S.
	LatticeFigureOfMerit(const BravaisType& brat, const SymMat43_Double& S, NRMat<Int4>& trans_mat);

	virtual ~LatticeFigureOfMerit(){};

	// Change the lattice constants to string.
	inline string printOptimizedLatticeConstants(const eABCaxis& axis1,
													const eRHaxis& axis2,
													const Int4& precision) const;

	inline const ePointGroup& enumPointGroup() const { return m_brat.enumLatticeSymmetry(); };
	inline const CenteringType::Enum& enumCentringType() const { return m_brat.enumCenteringType(); };
	inline const BravaisType::Enum& enumBravaisType() const { return m_brat.enumBravaisType(); };
	inline const BravaisType& putBravaisType() const { return m_brat; };

	// The returned matrix is Buerger-reduced matrix equivalent with m_S_red_optimized.
	inline const SymMat<Double>& putInverseOfBuergerReducedForm() const { return m_S_red; };
};

#endif /*LATTICEFIGUREOFMERIT_HH_*/
