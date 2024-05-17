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
#include <map>
#include "BravaisType.hh"
#include "FracMat.hh"

inline NRMat<Int4> put_transform_mat_from_primitive_to_base(const eABCaxis& paxis, const eBASEaxis& baxis)
{
	// Determine b2 such that { axis, b1, b2 } = {0, 1, 2}.
	const Int4 b2 = put_complement_set3((Int4)paxis, (Int4)baxis);

	NRMat<Int4> ans(3,3,0);
	//  base_axis : 0  0  1 = 0 1/2  1/2 ^{-1}
	//         b2 : 1  1  0   0 1/2 -1/2
	//      paxis : 1 -1  0   1   0    0
	ans[(size_t)baxis][2] = 1;
	ans[(size_t)b2][0] = 1;
	ans[(size_t)b2][1] = 1;
	ans[(size_t)paxis][0] = 1;
	ans[(size_t)paxis][1] = -1;
	return ans;
}


BravaisType::BravaisType(const pair<CenteringType::Enum, ePointGroup>& arg)
	: m_eblat(arg.first), m_epg(arg.second)
{
	if( m_eblat == CenteringType::Prim )
	{
		if( m_epg == Ci ) m_ecs = Triclinic;
		else if( m_epg == C2h_X ||  m_epg == C2h_Y ||  m_epg == C2h_Z ) m_ecs = Monoclinic_P;
		else if( m_epg == D2h ) m_ecs = Orthorhombic_P;
		else if( m_epg == D4h_Z ) m_ecs = Tetragonal_P;
		else if( m_epg == D31d_rho ) m_ecs = Rhombohedral;
		else if( m_epg == D6h ) m_ecs = Hexagonal;
		else if( m_epg == Oh ) m_ecs = Cubic_P;
		else assert(false);
	}
	else if( m_eblat == CenteringType::BaseA )
	{
		if( m_epg == C2h_Z ) m_ecs = Monoclinic_B;
//		else if( m_epg == D2h ) m_ecs = Orthorhombic_C;
		else assert(false);
	}
	else if( m_eblat == CenteringType::BaseB )
	{
		if( m_epg == C2h_X ) m_ecs = Monoclinic_B;
//		else if( m_epg == D2h ) m_ecs = Orthorhombic_C;
		else assert(false);
	}
	else if( m_eblat == CenteringType::BaseC )
	{
		if( m_epg == C2h_Y ) m_ecs = Monoclinic_B;
		else if( m_epg == D2h ) m_ecs = Orthorhombic_C;
		else assert(false);
	}
	else if( m_eblat == CenteringType::Inner )
	{
		if( m_epg == C2h_X ||  m_epg == C2h_Y ||  m_epg == C2h_Z ) m_ecs = Monoclinic_B;
		else if( m_epg == D2h ) m_ecs = Orthorhombic_I;
		else if( m_epg == D4h_Z ) m_ecs = Tetragonal_I;
		else if( m_epg == Oh ) m_ecs = Cubic_I;
		else assert(false);
	}
	else if( m_eblat == CenteringType::Face )
	{
		if( m_epg == D2h ) m_ecs = Orthorhombic_F;
		else if( m_epg == Oh ) m_ecs = Cubic_F;
		else assert(false);
	}
	else if( m_eblat == CenteringType::Rhom_hex )
	{
		if( m_epg == D3d_1_hex ) m_ecs = Rhombohedral;
		else assert(false);
	}
	else assert(false);
}


BravaisType::BravaisType(const BravaisType::Enum& ecs, const eABCaxis& abc_axis, const eRHaxis& rh_axis)
{
	m_ecs = ecs;
	if( ecs == Triclinic )
	{
		m_eblat = CenteringType::Prim;
		m_epg = Ci;
	}
	else if( ecs == Monoclinic_P )
	{
		m_eblat = CenteringType::Prim;
		if( abc_axis == A_Axis ) m_epg = C2h_X;
		else if( abc_axis == B_Axis ) m_epg = C2h_Y;
		else // if( abc_axis == C_Axis )
			m_epg = C2h_Z;
	}
	else if( ecs == Monoclinic_B )
	{
		if( abc_axis == A_Axis )
		{
			m_epg = C2h_X;
			m_eblat = CenteringType::BaseB;
		}
		else if( abc_axis == B_Axis )
		{
			m_epg = C2h_Y;
			m_eblat = CenteringType::BaseC;
		}
		else // if( abc_axis == C_Axis )
		{
			m_epg = C2h_Z;
			m_eblat = CenteringType::BaseA;
		}
	}
	else if( ecs == Orthorhombic_P )
	{
		m_epg = D2h;
		m_eblat = CenteringType::Prim;
	}
	else if( ecs == Orthorhombic_C )
	{
		m_epg = D2h;
		m_eblat = CenteringType::BaseC;
	}
	else if( ecs == Orthorhombic_I )
	{
		m_epg = D2h;
		m_eblat = CenteringType::Inner;
	}
	else if( ecs == Orthorhombic_F )
	{
		m_epg = D2h;
		m_eblat = CenteringType::Face;
	}
	else if( ecs == Tetragonal_P )
	{
		m_epg = D4h_Z;
		m_eblat = CenteringType::Prim;
	}
	else if( ecs == Tetragonal_I )
	{
		m_epg = D4h_Z;
		m_eblat = CenteringType::Inner;
	}
	else if( ecs == Rhombohedral )
	{
		if( rh_axis == Hex_Axis )
		{
			m_epg = D3d_1_hex;
			m_eblat = CenteringType::Rhom_hex;
		}
		else // if( rh_axis == Rho_Axis )
		{
			m_epg = D31d_rho;
			m_eblat = CenteringType::Prim;
		}
	}
	else if( ecs == Hexagonal )
	{
		m_epg = D6h;
		m_eblat = CenteringType::Prim;
	}
	else if( ecs == Cubic_P )
	{
		m_epg = Oh;
		m_eblat = CenteringType::Prim;
	}
	else if( ecs == Cubic_I )
	{
		m_epg = Oh;
		m_eblat = CenteringType::Inner;
	}
	else if( ecs == Cubic_F )
	{
		m_epg = Oh;
		m_eblat = CenteringType::Face;
	}
	else assert( false );
};


const NRMat<Int4>& BravaisType::putTransformMatrixFromPrimitiveToBase(const eBASEaxis& baxis)
{
	static const NRMat<Int4> prim_to_base_tray[3]
		                 = {
		                	 put_transform_mat_from_primitive_to_base(C_Axis, BaseA_Axis),
		                	 put_transform_mat_from_primitive_to_base(A_Axis, BaseB_Axis),
		                	 put_transform_mat_from_primitive_to_base(B_Axis, BaseC_Axis)
							};

	return prim_to_base_tray[(size_t)baxis];
}


const FracMat<Int4>& BravaisType::putTransformMatrixFromBaseToPrimitive(const eBASEaxis& baxis)
{
	static const FracMat<Int4> base_to_prim_tray[3]
		                 = { FInverse3<Int4>( putTransformMatrixFromPrimitiveToBase(BaseA_Axis) ),
		                	 FInverse3<Int4>( putTransformMatrixFromPrimitiveToBase(BaseB_Axis) ),
		                	 FInverse3<Int4>( putTransformMatrixFromPrimitiveToBase(BaseC_Axis) ) };

	return base_to_prim_tray[(size_t)baxis];
}


//             1  1  0
// Returns T = 1  1  2
//             1 -1  0
NRMat<Int4> BravaisType::putTransformMatrixFromPrimitiveToFace()
{
	NRMat<Int4> ans(3,3);
	ans[0][0] = 1;
	ans[0][1] = 1;
	ans[0][2] = 0;
	ans[1][0] = 1;
	ans[1][1] = 1;
	ans[1][2] = 2;
	ans[2][0] = 1;
	ans[2][1] = -1;
	ans[2][2] = 0;
	return ans;
}

//             1  0  1    0.5  0.5  0 ^{-1}
// Returns T = 1  0 -1  = 0.5  0.5  1
//            -1  1  0    0.5 -0.5  0
NRMat<Int4> BravaisType::putTransformMatrixFromBodyToPrimitive()
{
	NRMat<Int4> ans(3,3);
	ans[0][0] = 1;
	ans[0][1] = 0;
	ans[0][2] = 1;
	ans[1][0] = 1;
	ans[1][1] = 0;
	ans[1][2] =-1;
	ans[2][0] =-1;
	ans[2][1] = 1;
	ans[2][2] = 0;
	return ans;
}

//             -1  0  1
// Returns T =  0  1 -1
//             -1 -1 -1
NRMat<Int4> BravaisType::putTransformMatrixFromPrimitiveToRhomHex()
{
	NRMat<Int4> ans(3,3);
	ans[0][0] = -1;
	ans[0][1] = 0;
	ans[0][2] = 1;
	ans[1][0] = 0;
	ans[1][1] = 1;
	ans[1][2] = -1;
	ans[2][0] = -1;
	ans[2][1] = -1;
	ans[2][2] = -1;
	return ans;
}

string BravaisType::ToString(const BravaisType::Enum& e, const eABCaxis& axis)
{
	static const size_t NUM_LS = 14;
	static const string name[NUM_LS] = { "Triclinic", "Monoclinic(P)", "Monoclinic",
											"Orthorhombic(P)", "Orthorhombic(C)",
											"Orthorhombic(I)", "Orthorhombic(F)",
											"Tetragonal(P)", "Tetragonal(I)",
											"Trigonal", "Hexagonal",
											"Cubic(P)", "Cubic(I)", "Cubic(F)" };

	if( e == Monoclinic_B )
	{
		if( axis == A_Axis ) return name[(size_t)e] + "(B)";
		else if( axis == B_Axis ) return name[(size_t)e] + "(C)";
		else // if( axis == C_Axis )
			return name[(size_t)e] + "(A)";
	}
	else return name[(size_t)e];
}
