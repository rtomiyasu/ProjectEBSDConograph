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

#ifndef BRAVAIS_LAT_HH_
#define BRAVAIS_LAT_HH_

#include "../ZAnalysisTypes.hh"
#include "ICentringType.hh"

// Class of the simple lattice.
class SimpleLat : virtual public ICentringType
{
private:
	static const Int4 m_colattice_basis[3][3];

public:
	SimpleLat();
	const string& Name() const;
	const VecDat3<S1>* putTranslationVector(Int4&) const;
	Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const;
	virtual ~SimpleLat();
};


// Class of the base-centered lattice.
class BaseCenteredLat : virtual public ICentringType
{
private:
	static const Int4 m_colattice_basis[3][3][3];
	const eBASEaxis m_paxis; // m_paxis : principal axis. (0:a-axis, 1:b-axis, 2:c-axis)
	
public:
	BaseCenteredLat(const eBASEaxis&);
	const string& Name() const;
	const VecDat3<S1>* putTranslationVector(Int4&) const;
	Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const;
	virtual ~BaseCenteredLat();
};


// Class of the face-centered lattice.
class FaceCenteredLat : virtual public ICentringType
{
private:
	static const Int4 m_colattice_basis[3][3];

public:
	FaceCenteredLat();
	const string& Name() const;
	const VecDat3<S1>* putTranslationVector(Int4&) const;
	Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const;
	virtual ~FaceCenteredLat();
};


// Class of the body-centered lattice.
class BodyCenteredLat : virtual public ICentringType
{
private:
	static const Int4 m_colattice_basis[3][3];

public:
	BodyCenteredLat();
	const string& Name() const;
	const VecDat3<S1>* putTranslationVector(Int4&) const;
	Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const;
	virtual ~BodyCenteredLat();
};

// Class of the rhombohedral lattice(Hexagonal).
class RhombohedralLat : virtual public ICentringType
{
private:
	static const Int4 m_colattice_basis[3][3];
	static const Int4 m_paxis=2; // principal axis. (c-axis)

public:
	RhombohedralLat();
	const string& Name() const;
	const VecDat3<S1>* putTranslationVector(Int4&) const;
	Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3& hkl) const;
	virtual ~RhombohedralLat();
};

#endif /*BRAVAIS_LAT_HH_*/
