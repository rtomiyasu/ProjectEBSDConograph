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
#ifndef CENTRINGTYPE_HH_
#define CENTRINGTYPE_HH_

#include "../ZAnalysisTypes.hh"
#include "ICentringType.hh"

class SimpleLat;
class BaseCenteredLat;
class FaceCenteredLat;
class BodyCenteredLat;
class RhombohedralLat;

class CentringType
{
private:
	static const SimpleLat blat_Primitive;
	static const BaseCenteredLat blat_Base[3];
	static const FaceCenteredLat blat_Face;
	static const BodyCenteredLat blat_Body;
	static const RhombohedralLat blat_Rhom;

	const ICentringType* m_blat;

public:
	CentringType(const CenteringType::Enum&);
	CentringType(const CentringType&);
	CentringType& operator=(const CentringType& rhs);
	inline const ICentringType* operator->() const;
	virtual ~CentringType();
};

inline const ICentringType* CentringType::operator->() const
{
	return m_blat;
}

#endif /*CentringType_H_*/
