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

#include <assert.h>
#include "CentringType.hh"
#include "bravais_lat.hh"

const SimpleLat CentringType::blat_Primitive;
const BaseCenteredLat CentringType::blat_Base[3] = {BaseA_Axis, BaseB_Axis, BaseC_Axis};
const FaceCenteredLat CentringType::blat_Face;
const BodyCenteredLat CentringType::blat_Body;
const RhombohedralLat CentringType::blat_Rhom;

CentringType::CentringType(const CenteringType::Enum& num)
{
	m_blat = NULL;
	if(num==CenteringType::Prim) m_blat= &blat_Primitive;
	else if(num==CenteringType::BaseA) m_blat= &blat_Base[(size_t)BaseA_Axis];
	else if(num==CenteringType::BaseB) m_blat= &blat_Base[(size_t)BaseB_Axis];
	else if(num==CenteringType::BaseC) m_blat= &blat_Base[(size_t)BaseC_Axis];
	else if(num==CenteringType::Face) m_blat= &blat_Face;
	else if(num==CenteringType::Inner) m_blat= &blat_Body;
	else if(num==CenteringType::Rhom_hex) m_blat= &blat_Rhom;
	else assert( false );
}

CentringType::CentringType(const CentringType& rhs)
{
	m_blat=rhs.m_blat;
}

CentringType& CentringType::operator=(const CentringType& rhs)
{
	if(this != & rhs) m_blat = rhs.m_blat;
	return *this;
}

CentringType::~CentringType()
{
}
