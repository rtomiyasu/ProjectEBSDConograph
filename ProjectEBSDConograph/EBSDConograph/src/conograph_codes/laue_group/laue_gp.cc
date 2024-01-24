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
#include <stdlib.h>
#include "laue_gp.hh"


static const Int4 ldim = 3;
static const Int4 PNLatticeConst = 6;

// Class of triclinic lattice.
LaueTriclinic::LaueTriclinic()
{
}

string LaueTriclinic::CrystalSystemName() const
{
	static const string name = "Triclinic";
	return name;
}

string LaueTriclinic::Name() const
{
	static const string name = "-1";
	return name;
}

Int4 LaueTriclinic::LaueMultiplicity(const MillerIndex3& hkl) const
{
	return 2;
}

bool LaueTriclinic::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	for(Int4 i=0; i<ldim; i++)
	{
		if( length_axis[i] <= 0.0 )
		{
			flag = false;
			length_axis[i] = 1.0;
		}
		if( angle_axis[i] <= 0.0 || angle_axis[i] >= 180.0 )
		{
			flag = false;
			angle_axis[i] = 90.0;
		}
	}
	return flag;
}

LaueTriclinic::~LaueTriclinic()
{
}


// Class of monoclinic lattice.
LaueMonoclinic::LaueMonoclinic(const eABCaxis& num):m_paxis(num)
{
}

string LaueMonoclinic::CrystalSystemName() const
{
	static const string name = "Monoclinic";
	return name + "("+str_axis(m_paxis)+ ")";
}

string LaueMonoclinic::Name() const
{
	static const string name = "2/m";
	return name + "("+str_axis(m_paxis)+ ")";
}

Int4 LaueMonoclinic::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(hkl[m_paxis]==0) return 2;
	if( hkl[npaxis(0, m_paxis)] == 0 && hkl[npaxis(1, m_paxis)] == 0 ) return 2;
	else return 4;
}

bool  LaueMonoclinic::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	for(Int4 i=0; i<ldim; i++)
	{
		if( length_axis[i] <= 0.0 )
		{
			length_axis[i] = 1.0;
			flag = false;
		}
	}

	if( angle_axis[ m_paxis ] <= 0.0 || angle_axis[ m_paxis ] >= 180.0 )
	{
		angle_axis[ m_paxis ] = 90.0;
		flag = false;
	}
	if( fabs( angle_axis[ npaxis(0, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	if( fabs( angle_axis[ npaxis(1, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	angle_axis[ npaxis(0, m_paxis) ] = 90.0;
	angle_axis[ npaxis(1, m_paxis) ] = 90.0;
	return flag;
}

LaueMonoclinic::~LaueMonoclinic()
{
}


// Class of orthorhombic lattice.
LaueOrthorhombic::LaueOrthorhombic()
{
}

string LaueOrthorhombic::CrystalSystemName() const
{
	static const string name = "Orthorhombic";
	return name;
}

string LaueOrthorhombic::Name() const
{
	static const string name = "mmm";
	return name;
}

Int4 LaueOrthorhombic::LaueMultiplicity(const MillerIndex3& hkl) const
{
	Int4 count = 0;
	for(Int4 i=0; i<ldim; i++)
		if(hkl[i]==0) count++;
	if(count == 0) return 8;
	else if(count == 1) return 4;
	else return 2;
}

bool  LaueOrthorhombic::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	for(Int4 i=0; i<ldim; i++)
	{
		if( length_axis[i] <= 0.0 )
		{
			length_axis[i] = 1.0;
			flag = false;
		}
		if( fabs( angle_axis[i] - 90.0 ) > thred * 90.0 ) flag = false;
		angle_axis[i] = 90.0;
	}
	return flag;
}

LaueOrthorhombic::~LaueOrthorhombic()
{
}


// Class of tetragonal lattice.
LaueTetragonal::LaueTetragonal(const eABCaxis& num) : m_paxis(num)
{
}

string LaueTetragonal::CrystalSystemName() const
{
	static const string name = "Tetragonal";
	return name + "("+str_axis(m_paxis)+ ")";
}

string LaueTetragonal::Name() const
{
	static const string name = "4/m";
	return name + "("+str_axis(m_paxis)+ ")";
}

Int4 LaueTetragonal::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(hkl[m_paxis]==0) return 4;
	else if(hkl[npaxis(0, m_paxis)]==0 && hkl[npaxis(1, m_paxis)]==0) return 2;
	else return 8;
}

bool  LaueTetragonal::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	
	if( length_axis[m_paxis] <= 0.0 )
	{
		length_axis[m_paxis] = 1.0;
		flag = false;
	}
	if( length_axis[ npaxis(1, m_paxis) ] <= 0.0 )
	{
		if( length_axis[ npaxis(0, m_paxis) ] > 0.0 ) length_axis[ npaxis(1, m_paxis) ] = length_axis[ npaxis(0, m_paxis) ];
		else length_axis[ npaxis(1, m_paxis) ] = 1.0;
		flag = false;
	}
	if( fabs( length_axis[ npaxis(0, m_paxis) ] - length_axis[ npaxis(1, m_paxis) ] ) > thred * length_axis[ npaxis(1, m_paxis) ] ) flag = false;
	length_axis[ npaxis(0, m_paxis) ] = length_axis[ npaxis(1, m_paxis) ];
	
	for(Int4 i=0; i<ldim; i++)
	{
		if( fabs( angle_axis[i] - 90.0 ) > thred * 90.0 ) flag = false;
		angle_axis[i] = 90.0;
	}
	return flag;
}

LaueTetragonal::~LaueTetragonal()
{
}


// Class of tetragonal(4/mmm) lattice.
LaueTetragonal_mmm::LaueTetragonal_mmm(const eABCaxis& num) : LaueTetragonal(num)
{
}


string LaueTetragonal_mmm::Name() const
{
	static const string name = "4/mmm";
	return name;
}


Int4 LaueTetragonal_mmm::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(hkl[m_paxis]==0)
	{
		if( hkl[npaxis(0, m_paxis)]==0
				|| abs(hkl[npaxis(0, m_paxis)])==abs(hkl[npaxis(1, m_paxis)])
				|| hkl[npaxis(1, m_paxis)]==0 ) return 4;
		else return 8;
	}
	else if( hkl[npaxis(0, m_paxis)]==0 )
	{
		if( hkl[npaxis(1, m_paxis)]==0 ) return 2;
		return 8;
	}
	else if( hkl[npaxis(1, m_paxis)]==0
			|| abs(hkl[npaxis(0, m_paxis)])==abs(hkl[npaxis(1, m_paxis)]) ) return 8;
	else return 16;
}

LaueTetragonal_mmm::~LaueTetragonal_mmm()
{
}


// Class of cubic(m-3) lattice.
LaueCubic::LaueCubic()
{
}

string LaueCubic::CrystalSystemName() const
{
	static const string name = "Cubic";
	return name;
}

string LaueCubic::Name() const
{
	static const string name = "m-3";
	return name;
}

Int4 LaueCubic::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if( abs(hkl[0])==abs(hkl[1]) && abs(hkl[1])==abs(hkl[2]) ) return 8;
	Int4 count = 0;
	for(Int4 i=0; i<ldim; i++)
		if(hkl[i]==0) count++;
	if(count == 0) return 24;
	else if(count == 1) return 12;
	else return 6;
}

bool  LaueCubic::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	
	if( length_axis[2] <= 0.0 )
	{
		if( length_axis[0] > 0.0 ) length_axis[2] = length_axis[0];
		else if	( length_axis[1] > 0.0 ) length_axis[2] = length_axis[1];
		else length_axis[2] = 1.0;
		flag = false;
	}
	if( fabs( length_axis[0] - length_axis[2] ) > thred * length_axis[2] ) flag = false;
	if( fabs( length_axis[1] - length_axis[2] ) > thred * length_axis[2] ) flag = false;
	length_axis[0] = length_axis[2];
	length_axis[1] = length_axis[2];
	
	for(Int4 i=0; i<ldim; i++)
	{
		if( fabs( angle_axis[i] - 90.0 ) > thred * 90.0 ) flag = false;
		angle_axis[i] = 90.0;
	}
	
	return flag;
}

LaueCubic::~LaueCubic()
{
}


// Class of cubic(m-3m) lattice.
LaueCubic_m3m::LaueCubic_m3m()
{
}

string LaueCubic_m3m::Name() const
{
	static const string name = "m-3m";
	return name;
}

Int4 LaueCubic_m3m::LaueMultiplicity(const MillerIndex3& hkl) const
{
	Int4 count = 0;
	for(Int4 i=0; i<ldim; i++)
		if(hkl[i]==0) count++;
	if(count == 2) return 6;
	else if(count == 1){
		if( abs(hkl[0])==abs(hkl[1]) || abs(hkl[1])==abs(hkl[2]) || abs(hkl[0])==abs(hkl[2]) ) return 12;
		else return 24;
	}
	else if( abs(hkl[0])==abs(hkl[1]) ){
		if( abs(hkl[1])==abs(hkl[2]) ) return 8;
		else return 24;
	}
	else if( abs(hkl[0])==abs(hkl[2]) || abs(hkl[1])==abs(hkl[2]) ) return 24;
	else return 48;
}

LaueCubic_m3m::~LaueCubic_m3m()
{
}



const eABCaxis LaueHexagonal::m_paxis = C_Axis;

// Class of hexagonal(6/m) lattice.
LaueHexagonal::LaueHexagonal()
{
}

string LaueHexagonal::CrystalSystemName() const
{
	static const string name = "Hexagonal";
	return name;
}

string LaueHexagonal::Name() const
{
	static const string name = "6/m";
	return name;
}

Int4 LaueHexagonal::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if( hkl[2]==0 ) return 6;
	else if( hkl[0]==0 &&  hkl[1]==0 ) return 2;
	else return 12; 
}

  
bool  LaueHexagonal::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	
	if( length_axis[m_paxis] <= 0.0 )
	{
		length_axis[m_paxis] = 1.0;
		flag = false;
	}
	if( length_axis[ npaxis(1, m_paxis) ] <= 0.0 )
	{
		if( length_axis[ npaxis(0, m_paxis) ] > 0.0 ) length_axis[ npaxis(1, m_paxis) ] = length_axis[ npaxis(0, m_paxis) ];
		else length_axis[ npaxis(1, m_paxis) ] = 1.0;
		flag = false;
	}
	if( fabs( length_axis[ npaxis(0, m_paxis) ] - length_axis[ npaxis(1, m_paxis) ] ) > thred * length_axis[ npaxis(1, m_paxis) ] ) flag = false;
	length_axis[ npaxis(0, m_paxis) ] = length_axis[ npaxis(1, m_paxis) ];
	
	if( fabs( angle_axis[ m_paxis ] - 120.0 ) > thred * 120.0 ) flag = false;
	if( fabs( angle_axis[ npaxis(0, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	if( fabs( angle_axis[ npaxis(1, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	angle_axis[ m_paxis ] = 120.0;
	angle_axis[ npaxis(0, m_paxis) ] = 90.0;
	angle_axis[ npaxis(1, m_paxis) ] = 90.0;
	
	return flag;
}

LaueHexagonal::~LaueHexagonal()
{
}


// Class of Hexagonal(6/mmm) lattice.
LaueHexagonal_mmm::LaueHexagonal_mmm()
{
}

string LaueHexagonal_mmm::Name() const
{
	static const string name = "6/mmm";
	return name;
}

// hkl[3] = -hkl[0]-hkl[1](i.e. hkli)
Int4 LaueHexagonal_mmm::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(hkl[2]==0){
		if( hkl[0]==0 || hkl[1]==0 || hkl[3]==0 ) return 6;
		if( hkl[0]==hkl[1] || hkl[1]==hkl[3] || hkl[0]==hkl[3] ) return 6;
		else return 12;
	}
	else if( hkl[0]==0 ){
		if( hkl[1]==0 ) return 2;
		else return 12;
	}
	else if	( hkl[1]==0 || hkl[3]==0
				|| hkl[0]==hkl[1] || hkl[1]==hkl[3]  || hkl[0]==hkl[3] ) return 12;
	else return 24;
}

LaueHexagonal_mmm::~LaueHexagonal_mmm()
{
}





// Class of rhombohedral(-3) lattice.
LaueRhombohedral_rho::LaueRhombohedral_rho()
{
}

string LaueRhombohedral_rho::CrystalSystemName() const
{
	static const string name = "Trigonal(Rhombohedral-axis)";
	return name;
}

string LaueRhombohedral_rho::Name() const
{
	static const string name = "-3";
	return name;
}

Int4 LaueRhombohedral_rho::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(hkl[0]==hkl[1] && hkl[1]==hkl[2]) return 2;
	else return 6; 
}

bool  LaueRhombohedral_rho::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;
	
	if( length_axis[2] <= 0.0 )
	{
		if( length_axis[0] > 0.0 ) length_axis[2] = length_axis[0];
		else if	( length_axis[1] > 0.0 ) length_axis[2] = length_axis[1];
		else length_axis[2] = 1.0;
		flag = false;
	}
	if( fabs(length_axis[0] - length_axis[2]) > thred * length_axis[2] ) flag = false;
	if( fabs(length_axis[1] - length_axis[2]) > thred * length_axis[2] ) flag = false;
	length_axis[0] = length_axis[2];
	length_axis[1] = length_axis[2];

	if( angle_axis[2] <= 0.0 || angle_axis[2] >= 180.0 )
	{
		if( angle_axis[0] <= 0.0 || angle_axis[0] >= 180.0 ) angle_axis[2] = angle_axis[0];
		else if	( angle_axis[1] <= 0.0 || angle_axis[1] >= 180.0 ) angle_axis[2] = angle_axis[1];
		else angle_axis[2] = 90.0;
		flag = false;
	}
	if( fabs(angle_axis[0] - angle_axis[2]) > thred * angle_axis[2] ) flag = false;
	if( fabs(angle_axis[1] - angle_axis[2]) > thred * angle_axis[2] ) flag = false;
	angle_axis[0] = angle_axis[2];
	angle_axis[1] = angle_axis[2];
	
	return flag;
}

LaueRhombohedral_rho::~LaueRhombohedral_rho()
{
}



const eABCaxis LaueRhombohedral_hex::m_paxis = C_Axis;

// hex is the parameter that gives the coordinate axis.
LaueRhombohedral_hex::LaueRhombohedral_hex()
{
}

string LaueRhombohedral_hex::CrystalSystemName() const
{
	static const string name = "Trigonal(Hexagonal-axis)";
	return name;
}

string LaueRhombohedral_hex::Name() const
{
	static const string name = "-3";
	return name;
}

Int4 LaueRhombohedral_hex::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if( hkl[0]==0 && hkl[1]==0 ) return 2;
	else return 6; 
}

bool  LaueRhombohedral_hex::checkLatticeConstant(VecDat3<Double>& length_axis, VecDat3<Double>& angle_axis, Double thred) const
{
	bool flag = true;

	if( length_axis[m_paxis] <= 0.0 )
	{
		length_axis[m_paxis] = 1.0;
		flag = false;
	}
	if( length_axis[ npaxis(1, m_paxis) ] <= 0.0 )
	{
		if( length_axis[ npaxis(0, m_paxis) ] > 0.0 ) length_axis[ npaxis(1, m_paxis) ] = length_axis[ npaxis(0, m_paxis) ];
		else length_axis[ npaxis(1, m_paxis) ] = 1.0;
		flag = false;
	}
	if( fabs( length_axis[ npaxis(0, m_paxis) ] - length_axis[ npaxis(1, m_paxis) ] ) > thred * length_axis[ npaxis(1, m_paxis) ] ) flag = false;
	length_axis[ npaxis(0, m_paxis) ] = length_axis[ npaxis(1, m_paxis) ];
	
	if( fabs( angle_axis[ m_paxis ] - 120.0 ) > thred * 120.0 ) flag = false;
	if( fabs( angle_axis[ npaxis(0, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	if( fabs( angle_axis[ npaxis(1, m_paxis) ] - 90.0 ) > thred * 90.0 ) flag = false;
	
	angle_axis[ m_paxis ] = 120.0;
	angle_axis[ npaxis(0, m_paxis) ] = 90.0;
	angle_axis[ npaxis(1, m_paxis) ] = 90.0;
	
	return flag;
}

LaueRhombohedral_hex::~LaueRhombohedral_hex()
{
}



// Class of rhombohederal(-3m) lattice(rhombohederal axis).
LaueRhombohedral_rho3m::LaueRhombohedral_rho3m()
{
}

string LaueRhombohedral_rho3m::Name() const
{
	static const string name = "-3m";
	return name;
}

Int4 LaueRhombohedral_rho3m::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if( hkl[0]==hkl[1]){
		if( hkl[1]==hkl[2] ) return 2;
		else return 6;
	}
	else if( hkl[1]==hkl[2] || hkl[0]==hkl[2] ) return 6;
	else if( hkl[0]+hkl[1]+hkl[2]==0){
		for(Int4 i=0; i<ldim; i++)
			if(hkl[i]==0) return 6;
		return 12;
	}
	else return 12; 
}

LaueRhombohedral_rho3m::~LaueRhombohedral_rho3m()
{
}



// Class of rhombohederal(-31m, -3m1) lattice(hexagonal axis).
LaueRhombohedral_hex3m::LaueRhombohedral_hex3m(const Int4& num) : m_flag(num)
{
	assert( 0 <= num || num < 2 );
}

string LaueRhombohedral_hex3m::Name() const
{
	static const string name[] = {"-3.m","-3m."};
	return name[m_flag];
}

// hkl[3] = -hkl[0]-hkl[1](i.e. hkli)
Int4 LaueRhombohedral_hex3m::LaueMultiplicity(const MillerIndex3& hkl) const
{
	if(m_flag){
		if( hkl[0]==0 ){
			if( hkl[1]==0 ) return 2;
			else return 6;
		}
		else if( hkl[1]==0 || hkl[3]==0 ) return 6;
		else if( hkl[2]==0){
			if( hkl[0]==hkl[1] || hkl[0]==hkl[3] || hkl[1]==hkl[3] ) return 6;
			else return 12;
		}
		else return 12; 
	}
	else{
		if( hkl[0]==hkl[1]){
			if( hkl[1]==hkl[3] ) return 2;
			else return 6;
		}
		else if( hkl[1]==hkl[3] || hkl[0]==hkl[3] ) return 6;
		else if( hkl[2]==0 ){
			if( hkl[0]==0 || hkl[1]==0 || hkl[3]==0 ) return 6;
			return 12;
		}
		else return 12; 
	}
}

LaueRhombohedral_hex3m::~LaueRhombohedral_hex3m()
{
}
