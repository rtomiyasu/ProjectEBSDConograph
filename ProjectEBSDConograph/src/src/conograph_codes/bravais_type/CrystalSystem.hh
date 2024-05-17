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

#ifndef CRYSTALSYSTEM_HH_
#define CRYSTALSYSTEM_HH_

#include "../point_group/enumPointGroup.hh"
#include "../symmetric_operation/S1.hh"

struct CrystalSystem {
	enum Enum {
		Cubic        = 0,
		Tetragonal   = 1,
		Orthorhombic = 2,
		Hexagonal    = 3,
		Trigonal     = 4,
		Monoclinic   = 5,
		Triclinic    = 6,
	};

//	static Enum From_PhaseData(const vector<Double> &latticeParam, const NRMat<Double>& transf_P, const NRMat<Int4>& transf_P_inv, vector< NRMat<Int4> >& rotations);
	static string ToString(const Enum e);
};

#endif /*CRYSTALSYSTEM_HH_*/
