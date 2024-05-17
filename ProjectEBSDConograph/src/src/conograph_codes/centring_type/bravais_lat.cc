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

#include "bravais_lat.hh"
#include "../symmetric_operation/translation_vector.hh"

SimpleLat::SimpleLat()
{
}

const string& SimpleLat::Name() const
{
	static const string name = "P";
	return name;
}

const VecDat3<S1>* SimpleLat::putTranslationVector(Int4& num) const
{
	static const Int4 isize = 1;
	static const VecDat3<S1> vec_array[isize] = { exp2PIi000000() };
	
	num = isize;
	return vec_array;
}

Int4 SimpleLat::putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const
{
	return 1;
}

SimpleLat::~SimpleLat()
{
}


// Class of the base-centered lattice.
const Int4 BaseCenteredLat::m_colattice_basis[3][3][3]={ { {1,0,0},{0,1,1},{0,1,-1} },
														{ {1,0,1},{0,1,0},{1,0,-1} },
														{ {1,1,0},{1,-1,0},{0,0,1} } };

BaseCenteredLat::BaseCenteredLat(const eBASEaxis& num) : m_paxis(num)
{
}

const string& BaseCenteredLat::Name() const
{
	static const string name[] = {"A", "B", "C"};
	return name[(size_t)m_paxis];
}

const VecDat3<S1>* BaseCenteredLat::putTranslationVector(Int4& num) const
{
	static const Int4 isize = 2;
	static const VecDat3<S1> vec_array[][isize] = {	{exp2PIi000000(), exp2PIi001212()},
														{exp2PIi000000(), exp2PIi120012()},
														{exp2PIi000000(), exp2PIi121200()} };
	num = isize;
	return vec_array[(size_t)m_paxis];
}


inline Int4 npaxis(const Int4& b, const eBASEaxis& i)
{
	static const Int4 npaxis[][3] = { {1,2,0}, {2,0,1} };
	return npaxis[b][(size_t)i];
}

Int4 BaseCenteredLat::putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const
{
	if( (hkl[ npaxis(0,m_paxis) ]+hkl[ npaxis(1,m_paxis) ]) % 2 == 0 ) return 2;
	else return 0;
}

BaseCenteredLat::~BaseCenteredLat()
{
}


// Class of the face-centered lattice.
const Int4 FaceCenteredLat::m_colattice_basis[3][3]={ {1,1,1},{-1,1,-1},{-1,-1,1} };

FaceCenteredLat::FaceCenteredLat()
{
}

const string& FaceCenteredLat::Name() const
{
	static const string name = "F";
	return name;
}

const VecDat3<S1>* FaceCenteredLat::putTranslationVector(Int4& num) const
{
	static const Int4 isize = 4;
	static const VecDat3<S1> vec_array[isize] = { exp2PIi000000(), exp2PIi001212(), exp2PIi120012(), exp2PIi121200() };

	num = isize;
	return vec_array;
}

Int4 FaceCenteredLat::putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const
{
	Int4 evenflag = 0, oddflag = 0;
	for(Int4 i=0; i<3; i++){
		if(hkl[i] % 2==0) evenflag++;
		else oddflag++;
	}
	if(oddflag && evenflag) return 0;
	else return 4;
}

FaceCenteredLat::~FaceCenteredLat()
{
}


// Class of the body-centered lattice.
const Int4 BodyCenteredLat::m_colattice_basis[3][3]={ {0,1,1},{1,0,1},{1,1,0} };

BodyCenteredLat::BodyCenteredLat()
{
}

const string& BodyCenteredLat::Name() const
{
	static const string name = "I";
	return name;
}

const VecDat3<S1>* BodyCenteredLat::putTranslationVector(Int4& num) const
{
	static const Int4 isize = 2;
	static const VecDat3<S1> vec_array[isize] = { exp2PIi000000(), exp2PIi121212() };

	num = isize;
	return vec_array;
}

Int4 BodyCenteredLat::putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const
{
	Int4 sum=0;
	for(Int4 i=0; i<3; i++) sum+=hkl[i];
	if(sum % 2==0) return 2;
	else return 0;
}

BodyCenteredLat::~BodyCenteredLat()
{
}


// Class of the rhombohedral lattice.
const Int4 RhombohedralLat::m_colattice_basis[3][3]={ {1,-1,-1},{1,1,0},{1,0,1} };

RhombohedralLat::RhombohedralLat()
{
}

const string& RhombohedralLat::Name() const
{
	static const string name = "R(Hexagonal Axis)";
	return name;
}

const VecDat3<S1>* RhombohedralLat::putTranslationVector(Int4& num) const
{
	static const Int4 isize = 3;
	static const VecDat3<S1> vec_array[isize] = { exp2PIi000000(), exp2PIi231313(), exp2PIi132323() };

	num = isize;
	return vec_array;
}

Int4 RhombohedralLat::putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const
{
	if( (-hkl[0] + hkl[1] + hkl[2]) % 3 == 0 ) return 3;
	else return 0;
}

RhombohedralLat::~RhombohedralLat()
{
}
