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
#ifndef LAUE_GR_HH_
#define LAUE_GR_HH_

#include "../ZAnalysisTypes.hh"
#include "../bravais_type/enumAxis.hh"
#include "ILaueGroup.hh"

// Class of the triclinic lattice.
class LaueTriclinic : public ILaueGroup
{
public:
	LaueTriclinic();
	string CrystalSystemName() const;
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueTriclinic();
};


// Class of the monoclinic lattice.
class LaueMonoclinic : public ILaueGroup
{
protected:
	const eABCaxis m_paxis; // m_paxis : principal axis. (0:a-axis, 1:b-axis, 2:c-axis)

public:
	LaueMonoclinic(const eABCaxis&);
	string CrystalSystemName() const;
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueMonoclinic();
};


// Class of the orthorhombic lattice.
class LaueOrthorhombic : public ILaueGroup
{
public:
	LaueOrthorhombic();
	string CrystalSystemName() const;
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueOrthorhombic();
};


// Class of the tetragonal(4/m) lattice.
class LaueTetragonal : public ILaueGroup
{
protected:
	const eABCaxis m_paxis; // principal axis.

public:
	LaueTetragonal(const eABCaxis&);
	string CrystalSystemName() const;
	virtual string Name() const;
	virtual Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueTetragonal();
};


// Class of the tetragonal(4/mmm) lattice.
class LaueTetragonal_mmm : public LaueTetragonal
{
public:
	LaueTetragonal_mmm(const eABCaxis&);
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	virtual ~LaueTetragonal_mmm();
};


// Class of the cubic(m-3) lattice.
class LaueCubic : public ILaueGroup
{
public:
	LaueCubic();
	string CrystalSystemName() const;
	virtual string Name() const;
	virtual Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueCubic();
};


// Class of the cubic(m-3m) lattice.
class LaueCubic_m3m : public LaueCubic
{
public:
	LaueCubic_m3m();
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	virtual ~LaueCubic_m3m();
};

// Class of the Hexagonal(6m) lattice.
class LaueHexagonal : public ILaueGroup
{
protected:
	static const eABCaxis m_paxis;

public:
	LaueHexagonal();
	string CrystalSystemName() const;
	virtual string Name() const;
	virtual Int4 LaueMultiplicity(const MillerIndex3& hkil) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueHexagonal();
};


// Class of the Hexagonal(6/mmm) lattice.
class LaueHexagonal_mmm : public LaueHexagonal
{
public:
	LaueHexagonal_mmm();
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkil) const;
	virtual ~LaueHexagonal_mmm();
};


// Class of the rhombohedral(-3) lattice.
class LaueRhombohedral_rho : public ILaueGroup
{
public:
	LaueRhombohedral_rho();
	string CrystalSystemName() const;
	virtual string Name() const;
	virtual Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueRhombohedral_rho();
};

// Class of the rhombohedral(-3) lattice.
class LaueRhombohedral_hex : public ILaueGroup
{
protected:
	static const eABCaxis m_paxis;

public:
	LaueRhombohedral_hex();
	string CrystalSystemName() const;
	virtual string Name() const;
	virtual Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const;
	virtual ~LaueRhombohedral_hex();
};


// Class of the rhombohedral(-3m) lattice.
class LaueRhombohedral_rho3m : public LaueRhombohedral_rho
{
public:
	LaueRhombohedral_rho3m();
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	virtual ~LaueRhombohedral_rho3m();
};


// Class of the rhombohedral(-3m1, -31m) lattice.
class LaueRhombohedral_hex3m : public LaueRhombohedral_hex
{
protected:
	const Int4 m_flag; // 0:-31m,   1:-3m1

public:
	LaueRhombohedral_hex3m(const Int4&);
	string Name() const;
	Int4 LaueMultiplicity(const MillerIndex3& hkl) const;
	virtual ~LaueRhombohedral_hex3m();
};

#endif /*LAUE_GR_HH_*/
