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
#ifndef BRAVAISTYPE_HH_
#define BRAVAISTYPE_HH_

#include "../centring_type/enumCentringType.hh"
#include "../point_group/enumPointGroup.hh"
#include "enumAxis.hh"
#include "matrix_3by3.hh"
#include "FracMat.hh"

class BravaisType
{
public:
	enum Enum{ CSUndefined = -1, Triclinic = 0, Monoclinic_P = 1, Monoclinic_B = 2,
				Orthorhombic_P = 3, Orthorhombic_C = 4, Orthorhombic_I = 5, Orthorhombic_F = 6,
				Tetragonal_P = 7, Tetragonal_I = 8, Rhombohedral = 9, Hexagonal = 10,
				Cubic_P = 11, Cubic_I = 12, Cubic_F = 13 };

	static string ToString(const BravaisType::Enum& e, const eABCaxis& axis);

private:
	CenteringType::Enum m_eblat;
	ePointGroup m_epg;
	BravaisType::Enum m_ecs;

public:

	BravaisType() : m_eblat(CenteringType::Prim), m_epg(Ci), m_ecs(CSUndefined){};
	BravaisType(const pair<CenteringType::Enum, ePointGroup>& arg);
	BravaisType(const BravaisType::Enum& ecs, const eABCaxis& abc_axis, const eRHaxis& rh_axis);

	inline const ePointGroup& enumLatticeSymmetry() const { return m_epg; };
	inline const CenteringType::Enum& enumCenteringType() const { return m_eblat; };
	inline const BravaisType::Enum& enumBravaisType() const { return m_ecs; };
	inline string putBravaisTypeName() const { return ToString(m_ecs, enumABCaxis()); };

	inline eABCaxis enumABCaxis() const;
	inline eBASEaxis enumBASEaxis() const;
	inline eRHaxis enumRHaxis() const;

	inline bool operator==(const BravaisType& rhs) const { return m_eblat == rhs.m_eblat && m_epg == rhs.m_epg; };
	inline bool operator<(const BravaisType& rhs) const { return m_eblat < rhs.m_eblat || ( m_eblat == rhs.m_eblat && m_epg < rhs.m_epg ); };


	// for GUI
	BravaisType(CenteringType::Enum m_eblat, ePointGroup m_epg, BravaisType::Enum m_ecs) : m_eblat(m_eblat), m_epg(m_epg), m_ecs(m_ecs) {};

	//  ** Cell choice 1:
	//  base_axis : 0  0  1 *  a  c  d *  0  1  1  =  b    2d      0
	//         b2 : 1  1  0    c  a  d    0  1 -1    2d 2(a+c)     0
	//      paxis : 1 -1  0    d  d  b    1  0  0     0     0  2(a-c)
	static const NRMat<Int4>& putTransformMatrixFromPrimitiveToBase(const eBASEaxis& baxis);

	//  ** Cell choice 1:
	//  Returns inverse of the method putTransformMatrixFromPrimitiveToBase.
	static const FracMat<Int4>& putTransformMatrixFromBaseToPrimitive(const eBASEaxis& baxis);

	//             1  1  0    0.5   0  0.5 ^{-1}
	// Returns T = 1  1  2 =  0.5   0 -0.5
	//             1 -1  0   -0.5 0.5    0
	// such that
	//	   -2p-q     q     p           -4p       0       0
	// T *     q -2p-q     p * tr(T) =   0 -4(p+r)       0
	//	       p     p -2p-r             0       0 -4(p+q)
	static NRMat<Int4> putTransformMatrixFromPrimitiveToFace();

	//             1  0  1    0.5  0.5  0 ^{-1}
	// Returns T = 1  0 -1  = 0.5  0.5  1
	//            -1 -1  0    0.5 -0.5  0
	// such that
	//	        -2p-q     q     p                -p      0      0
	// T^{-1} *     q -2p-q     p * tr{T^{-1}} =  0 -(p+r)      0
	//	            p     p -2p-r                 0      0 -(p+q)
	static NRMat<Int4> putTransformMatrixFromBodyToPrimitive();

	//             -1  0  1
	// Returns T =  0  1 -1
	//             -1 -1 -1
	// such that
	//	   a d d           2a-2d   -a+d      0
	// T * d a d * tr(T) =  -a+d  2a-2d      0
	//	   d d a               0      0  3a+6d
	static NRMat<Int4> putTransformMatrixFromPrimitiveToRhomHex();
};

inline eABCaxis BravaisType::enumABCaxis() const
{
	if( m_epg == C2h_X ) return A_Axis;
	else if( m_epg == C2h_Y ) return B_Axis;
	else if( m_epg == C2h_Z ) return C_Axis;
	assert(false);
	return A_Axis;
}


inline eBASEaxis BravaisType::enumBASEaxis() const
{
	if( m_eblat == CenteringType::BaseA ) return BaseA_Axis;
	else if( m_eblat == CenteringType::BaseB ) return BaseB_Axis;
	else if( m_eblat == CenteringType::BaseC ) return BaseC_Axis;
	assert(false);
	return BaseA_Axis;
}


inline eRHaxis BravaisType::enumRHaxis() const
{
	if( m_ecs == Rhombohedral )
	{
		if( m_eblat == CenteringType::Rhom_hex ) return Hex_Axis;
		else return Rho_Axis;
	}
	else if( m_ecs == Hexagonal ) return Hex_Axis;
	assert(false);
	return Rho_Axis;
}



// Transform matrix from lhs to rhs. (This function assumes cell choice 1.)
inline NRMat<Int4> put_transform_matrix_monoclinic_b(const eABCaxis& lhs, const eABCaxis& rhs)
{
	if( lhs == rhs ) return put_matrix_XYZ();
	else if( lhs == A_Axis )
	{
		if( rhs == B_Axis ) return put_matrix_ZXY();
		else // if( rhs == C_Axis )
			 return put_matrix_YZX();
	}
	else if( lhs == B_Axis )
	{
		if( rhs == C_Axis ) return put_matrix_ZXY();
		else // if( rhs == A_Axis )
			 return put_matrix_YZX();
	}
	else // if( lhs == C_Axis )
	{
		if( rhs == A_Axis ) return put_matrix_ZXY();
		else // if( rhs == B_Axis )
			 return put_matrix_YZX();
	}
}


inline const ePointGroup& put_monoclinic_p_type(const eABCaxis& axis)
{
	static const ePointGroup brat_tray[3]
	     = { C2h_X, C2h_Y, C2h_Z };

	return brat_tray[(size_t)axis];
}


inline const pair<CenteringType::Enum, ePointGroup>& put_monoclinic_b_type(const eABCaxis& axis)
{
	static const pair<CenteringType::Enum, ePointGroup> brat_tray[3]
	     = { pair<CenteringType::Enum, ePointGroup>(CenteringType::BaseB,C2h_X),
	    	 pair<CenteringType::Enum, ePointGroup>(CenteringType::BaseC,C2h_Y),
	    	 pair<CenteringType::Enum, ePointGroup>(CenteringType::BaseA,C2h_Z) };

	// Cell choice 1.
	return brat_tray[(size_t)axis];
}


inline const pair<CenteringType::Enum, ePointGroup>& put_rhombohedral_type(const eRHaxis& axis)
{
	static const pair<CenteringType::Enum, ePointGroup> brat_tray[2]
	     = { pair<CenteringType::Enum, ePointGroup>(CenteringType::Prim,D31d_rho),
	    		 pair<CenteringType::Enum, ePointGroup>(CenteringType::Rhom_hex,D3d_1_hex) };

	return brat_tray[(size_t)axis];
}

#endif /*enumBravaisType_HH_*/
