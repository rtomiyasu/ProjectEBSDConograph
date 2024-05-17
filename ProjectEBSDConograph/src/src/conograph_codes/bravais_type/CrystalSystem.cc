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

#ifdef DEBUG
	#include <iostream>
#endif
#include <cassert>
#include "CrystalSystem.hh"
using namespace std;

//namespace {
///*
//	bool isEqual(const Double a, const Double b) { return abs(a - b) < 1.0e-8; }
//
//	NRMat<Int4> putIDMatrix()
//	{
//		NRMat<Int4> id(3, 3);
//		for(size_t i=0; i<3; i++) {
//			for(size_t j=0; j<3; j++) {
//				id[i][j] = (i == j) ? 1 : 0;
//			}
//		}
//		return id;
//	};
//
//	// 蠑墓焚縺ｮ蜈ｨ縺ｦ縺ｮ蜈ｨ縺ｦ縺ｮ陦悟�励↓蟇ｾ縺励��2荵励＠縺溘→縺阪↓1縺ｫ縺ｪ繧九↑繧液rue繧定ｿ斐☆
//	bool are1WhenSquared(const vector<NRMat<Int4> > &rotations) {
//		static const NRMat<Int4> id = putIDMatrix();
//		for(size_t i=0; i<rotations.size(); i++) {
//			NRMat<Int4> rot = rotations[i];
//			if(!(product(rot, rot)==id)) { return false; }
//		}
//		return true;
//	}
//
//	bool isHexagonal(const vector<NRMat<Int4> > &rotations) {
//		static const NRMat<Int4> id = putIDMatrix();
//		Int4 count = 0;
//		for(size_t i=0; i<rotations.size() && count <= 4; i++) {
//			NRMat<Int4> rot = rotations[i];
//			if(product(rot, rot)==id) { count++; }
//		}
//		if( count > 4 ) return false;
//		return true;
//	}
//*/
//
//	Int4 toInt(const Double d) {
//		Int4 sign = d < 0 ? -1 : 1;
//		Int4 body = Int4(abs(d) + 1.0e-5);
//		return sign * body;
//	}
//}
//
//
//CrystalSystem::Enum CrystalSystem::From_PhaseData(const vector<Double> &latticeParam, const NRMat<Double>& transf_P, const NRMat<Int4>& transf_P_inv, vector< NRMat<Int4> >& rotations) {
//	assert(latticeParam.size() == 6);
//
//	static const Double RadDeg = M_PI /180.0;
//
//	NRMat<Double> InvS(3, 3);
//	InvS[0][0] = latticeParam[0]*latticeParam[0];
//	InvS[1][1] = latticeParam[1]*latticeParam[1];
//	InvS[2][2] = latticeParam[2]*latticeParam[2];
//	InvS[0][1] = latticeParam[0]*latticeParam[1]*cos(latticeParam[5]*RadDeg);
//	InvS[0][2] = latticeParam[0]*latticeParam[2]*cos(latticeParam[4]*RadDeg);
//	InvS[1][2] = latticeParam[1]*latticeParam[2]*cos(latticeParam[3]*RadDeg);
//	InvS[1][0] = InvS[0][1];
//	InvS[2][0] = InvS[0][2];
//	InvS[2][1] = InvS[1][2];
//
//	vector< NRMat<Int4> > sym_opt_tray;
//	putAllSymmetryOperations(sym_opt_tray);
//
//	rotations.clear();
//	NRMat<Double> rot(3,3);
//	bool isSame;
//	for(size_t i=0; i<sym_opt_tray.size(); i++) {
//		const NRMat<Int4>& rot_i = sym_opt_tray[i];
//		rot[0][0] = rot_i[0][0];
//		rot[0][1] = rot_i[0][1];
//		rot[0][2] = rot_i[0][2];
//		rot[1][0] = rot_i[1][0];
//		rot[1][1] = rot_i[1][1];
//		rot[1][2] = rot_i[1][2];
//		rot[2][0] = rot_i[2][0];
//		rot[2][1] = rot_i[2][1];
//		rot[2][2] = rot_i[2][2];
//
//		const NRMat<Double> InvS2 = mprod(mprod(rot, InvS), transpose(rot));
//
//		isSame = true;
//		for(Int4 j=0; j<3 && isSame; j++) {
//			for(Int4 k=j; k<3; k++) {
//				if( fabs(InvS2[j][k] - InvS[j][k]) > max( max( InvS[j][j], InvS[k][k] ), max( InvS2[j][j], InvS2[k][k] ) )*1.0e-10 )
//				{
//					isSame = false;
//					break;
//				}
//			}
//		}
//		if( isSame ) {
//			// primitive cell縺ｮ繧ゅ�ｮ縺ｫ鄂ｮ縺肴鋤縺医ｋ縲�
///*
//for(int i=0; i<3; i++)
//{
//	for(int k=0; k<3; k++)
//	{
//cout << rot_i[i][k] << " ";
//	}
//cout << "\n";
//}
//cout << "\n";
//*/
//			NRMat<Int4> rot2_i( mprod(rot_i, transf_P_inv) );
//			rot[0][0] = rot2_i[0][0]; rot[0][1] = rot2_i[0][1]; rot[0][2] = rot2_i[0][2];
//			rot[1][0] = rot2_i[1][0]; rot[1][1] = rot2_i[1][1]; rot[1][2] = rot2_i[1][2];
//			rot[2][0] = rot2_i[2][0]; rot[2][1] = rot2_i[2][1]; rot[2][2] = rot2_i[2][2];
//			rot = mprod(transf_P, rot);
///*
//for(int i=0; i<3; i++)
//{
//	for(int k=0; k<3; k++)
//	{
//cout << rot[i][k] << " ";
//	}
//	cout << "\n";
//}
//cout << "\n";
//*/
//			rot2_i[0][0] = ::toInt(rot[0][0]); rot2_i[0][1] = ::toInt(rot[0][1]); rot2_i[0][2] = ::toInt(rot[0][2]);
//			rot2_i[1][0] = ::toInt(rot[1][0]); rot2_i[1][1] = ::toInt(rot[1][1]); rot2_i[1][2] = ::toInt(rot[1][2]);
//			rot2_i[2][0] = ::toInt(rot[2][0]); rot2_i[2][1] = ::toInt(rot[2][1]); rot2_i[2][2] = ::toInt(rot[2][2]);
//
//			isSame = true;
//			for(Int4 j=0; j<3 && isSame; j++)
//				for(Int4 k=0; k<3; k++)
//					if( fabs( rot2_i[j][k] - rot[j][k] ) > 1.0e-10 )
//					{
//						isSame = false;
//					}
//			if( isSame ) rotations.push_back(rot2_i);
//		}
//	}
//
//	switch(rotations.size()) {
//		case 2 : return Triclinic;
//		case 4 : return Monoclinic;
//		case 8 : return Orthorhombic;
//		case 12: return Trigonal;
//		case 16: return Tetragonal;
//		case 24: return Hexagonal;
//		case 48: return Cubic;
//		default: assert(false); return Triclinic;
//	}
//}

string CrystalSystem::ToString(const CrystalSystem::Enum e) {
	switch(e) {
		case Cubic       : return "Cubic";
		case Tetragonal  : return "Tetragonal";
		case Orthorhombic: return "Orthorhombic";
		case Hexagonal   : return "Hexagonal";
		case Trigonal    : return "Trigonal";
		case Monoclinic  : return "Monoclinic";
		case Triclinic   : return "Triclinic";
		default: assert(false); return "unknown";
	}
}
