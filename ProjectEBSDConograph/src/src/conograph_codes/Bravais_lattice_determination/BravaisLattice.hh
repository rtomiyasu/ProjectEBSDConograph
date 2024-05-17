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
#ifndef _BravaisLattice_hh_
#define _BravaisLattice_hh_
// BravaisLattice.hh

#include <map>

#include "../bravais_type/enumAxis.hh"
#include "../bravais_type/BravaisType.hh"
#include "../utility_data_structure/SymMat.hh"
#include "../ZAnalysisTypes.hh"

class LatticeFigureOfMeritToCheckSymmetry;
class ReducedLatticeToCheckBravais;

using namespace std;


class BravaisLattice
{
private:
	enum{ NUM_LS = 14 };	// 	Triclinic, Monoclinic, Monoclinic(C),
							// Orthorhombic(P), Orthorhombic(C), Orthorhombic(I), Orthorhombic(F),
							// Tetragonal(P), Tetragonal(I), Rhombohedral, Hexiagonl
							// Cubic(P), Cubic(I), Cubic(F).

	bool JudgeSymmetry[NUM_LS];
	bool m_DoesPrudentSymSearch;
    Double m_resol;

    vector<LatticeFigureOfMeritToCheckSymmetry> m_result[NUM_LS];

    static void putCentringTypes(const ReducedLatticeToCheckBravais& RLCB,
										const BravaisType& brat,
										vector<LatticeFigureOfMeritToCheckSymmetry>& lattice_result);

	void putLatticeCandidatesForEachBravaisTypes(
								const eABCaxis& baxis_flag,
								const eRHaxis& rhom_flag);

public:
    BravaisLattice();
    ~BravaisLattice();

    void setParam(const bool& does_exhaustive_search, const Double& resolution);

    // S_super is assumed to be Delone reduced.
    // Returns almost equivalent Bravais lattices after sorting them by the distance from S_super.
	void execute(const SymMat<Double>& S_super,
								const eABCaxis& baxis_flag,
								const eRHaxis& rhom_flag);

	const vector<LatticeFigureOfMeritToCheckSymmetry>& putResult(const BravaisType::Enum& arg) const { return m_result[arg]; };

	const LatticeFigureOfMeritToCheckSymmetry& putBravaisLattice() const;

	static void putLatticeSymmetries(const SymMat<Double>& S_super, vector< NRMat<Int4> >& rotations);
};

#endif
