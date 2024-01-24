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
#include <cassert>
#include "LaueGroup.hh"
#include"laue_gp.hh"

const LaueTriclinic LaueGroup::lg_Triclinic;
const LaueMonoclinic LaueGroup::lg_Monoclinic[3] = {A_Axis, B_Axis, C_Axis};
const LaueOrthorhombic LaueGroup::lg_Orthorhombic;
const LaueTetragonal LaueGroup::lg_Tetragonal[3] = {A_Axis, B_Axis, C_Axis};
const LaueTetragonal_mmm LaueGroup::lg_Tetragonal_mmm[3] = {A_Axis, B_Axis, C_Axis};
const LaueCubic LaueGroup::lg_Cubic;
const LaueCubic_m3m LaueGroup::lg_Cubic_m3m;
const LaueHexagonal LaueGroup::lg_Hexagonal;
const LaueHexagonal_mmm LaueGroup::lg_Hexagonal_mmm;
const LaueRhombohedral_rho LaueGroup::lg_Rhombohedral_rho;
const LaueRhombohedral_hex LaueGroup::lg_Rhombohedral_hex;
const LaueRhombohedral_rho3m LaueGroup::lg_Rhombohedral_rho3m;
const LaueRhombohedral_hex3m LaueGroup::lg_Rhombohedral_hex3m[2] = {false, true};

LaueGroup::LaueGroup(const ePointGroup& num)
{
	m_lg = NULL;
	if(num==Ci) m_lg= &lg_Triclinic;
	else if(num==C2h_X) m_lg= &lg_Monoclinic[(size_t)A_Axis];
	else if(num==C2h_Y) m_lg= &lg_Monoclinic[(size_t)B_Axis];
	else if(num==C2h_Z) m_lg= &lg_Monoclinic[(size_t)C_Axis];
	else if(num==D2h) m_lg= &lg_Orthorhombic;
	else if(num==C4h_X) m_lg= &lg_Tetragonal[(size_t)A_Axis];
	else if(num==C4h_Y) m_lg= &lg_Tetragonal[(size_t)B_Axis];
	else if(num==C4h_Z) m_lg= &lg_Tetragonal[(size_t)C_Axis];
	else if(num==D4h_X) m_lg= &lg_Tetragonal_mmm[(size_t)A_Axis];
	else if(num==D4h_Y) m_lg= &lg_Tetragonal_mmm[(size_t)B_Axis];
	else if(num==D4h_Z) m_lg= &lg_Tetragonal_mmm[(size_t)C_Axis];
	else if(num==Th) m_lg= &lg_Cubic;
	else if(num==Oh) m_lg= &lg_Cubic_m3m;
	else if(num==C6h) m_lg= &lg_Hexagonal;
	else if(num==D6h) m_lg= &lg_Hexagonal_mmm;
	else if(num==C31i_rho) m_lg= &lg_Rhombohedral_rho;
	else if(num==C3i_hex) m_lg= &lg_Rhombohedral_hex;
	else if(num==D31d_rho) m_lg= &lg_Rhombohedral_rho3m;
	else if(num==D3d_0_hex) m_lg= &lg_Rhombohedral_hex3m[0];
	else if(num==D3d_1_hex) m_lg= &lg_Rhombohedral_hex3m[1];
	else assert( false );
}

LaueGroup::LaueGroup(const LaueGroup& rhs)
{
	m_lg=rhs.m_lg;
}

LaueGroup& LaueGroup::operator=(const LaueGroup& rhs)
{
	if(this != & rhs) m_lg=rhs.m_lg;
	return *this;
}

LaueGroup::~LaueGroup()
{
}
