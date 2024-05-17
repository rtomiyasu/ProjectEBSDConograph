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
#ifndef LAUEGROUP_HH_
#define LAUEGROUP_HH_

#include"../ZAnalysisTypes.hh"
#include"../point_group/enumPointGroup.hh"
#include"ILaueGroup.hh"

class LaueTriclinic;
class LaueMonoclinic;
class LaueOrthorhombic;
class LaueTetragonal;
class LaueTetragonal_mmm;
class LaueCubic;
class LaueCubic_m3m;
class LaueHexagonal;
class LaueHexagonal_mmm;
class LaueRhombohedral_rho;
class LaueRhombohedral_hex;
class LaueRhombohedral_rho3m;
class LaueRhombohedral_hex3m;

// Class of the Laue Group.
class LaueGroup
{
private:
	static const LaueTriclinic lg_Triclinic;
	static const LaueMonoclinic lg_Monoclinic[3];
	static const LaueOrthorhombic lg_Orthorhombic;
	static const LaueTetragonal lg_Tetragonal[3];
	static const LaueTetragonal_mmm lg_Tetragonal_mmm[3];
	static const LaueCubic lg_Cubic;
	static const LaueCubic_m3m lg_Cubic_m3m;
	static const LaueHexagonal lg_Hexagonal;
	static const LaueHexagonal_mmm lg_Hexagonal_mmm;
	static const LaueRhombohedral_rho lg_Rhombohedral_rho;
	static const LaueRhombohedral_hex lg_Rhombohedral_hex;
	static const LaueRhombohedral_rho3m lg_Rhombohedral_rho3m;
	static const LaueRhombohedral_hex3m lg_Rhombohedral_hex3m[2];

	const ILaueGroup* m_lg;

public:
	LaueGroup(const ePointGroup&);
	LaueGroup(const LaueGroup&);
	LaueGroup& operator=(const LaueGroup& rhs);
	inline const ILaueGroup* operator->() const;
	~LaueGroup();
};

inline const ILaueGroup* LaueGroup::operator->() const
{
	return m_lg;
}

#endif /*LAUEGROUP_HH_*/
