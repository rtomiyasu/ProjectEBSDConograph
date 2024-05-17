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
#ifndef CHECK_EQUIV_HH_
#define CHECK_EQUIV_HH_

#include "../utility_data_structure/SymMat.hh"

inline bool equiv_resol(const Double& lhs, const Double& rhs,
		const Double& resol)
{
	Double diff = lhs - rhs;
	Double max_value = max( fabs(lhs), fabs(rhs) );
	if( fabs( diff ) <= max_value * resol ) return true;
	return false;
}

inline bool equiv_zero(const SymMat<Double>& S, const Int4& i, const Int4& j, const Double& resol)
{
	return equiv_resol( S(i,i)+S(j,j)+S(i,j)*2.0, S(i,i)+S(j,j)-S(i,j)*2.0, resol);
}

inline bool equiv_zero(const SymMat<Double>& S, const Int4& i, const Int4& j, const Int4& l, const Double& resol)
{
	// Ki, Kj, Kl, Km, Kn -> Ki, -Kj, Kl, Kj+Km, Kj+Kn.
	// (Sii + Sjj + Sij*2) + (Sjj + Sll + Sjl*2) + (Sii + Sjj + Sll + Sij*2 + Sjl*2 + Sil*2)
	//	= (Sii + Sjj - Sij*2) + (Sjj + Sll - Sjl*2) + (Sii + Sjj + Sll - Sij*2 - Sjl*2 + Sil*2).
	return equiv_resol(
			(S(i,i)+S(j,j)+S(l,l)+S(i,l)+(S(i,j)+S(j,l))*2.0)*2.0+S(j,j),
			(S(i,i)+S(j,j)+S(l,l)+S(i,l)-(S(i,j)+S(j,l))*2.0)*2.0+S(j,j), resol);
}

inline bool equiv_zero(const SymMat<Double>& S, const Int4& i, const Int4& j, const Int4& k, const Int4& l, const Double& resol)
{
	// Ki, Kj, Kk, Kl, Kn -> Ki+Kk, -Kj, -Kk, Kj+Kl, Kj+Kk+Kn
	// (Sii + Sjj + Sij*2) + (Skk + Sll + Skl*2) + (Sii + Sjj + Skk + Sij*2 + Sjk*2 + Sik*2) + (Sjj + Skk + Sll + Sjk*2 + Skl*2 + Sjl*2)
	//	= (Sii + Sjj - Sij*2) + (Skk + Sll - Skl*2) + (Sii + Sjj + Skk - Sij*2 - Sjk*2 + Sik*2) + (Sjj + Skk + Sll - Sjk*2 - Skl*2 + Sjl*2).
	return equiv_resol(
			(S(i,i)+S(j,j)+S(k,k)+S(l,l)+S(i,k)+S(j,l)+(S(i,j)+S(j,k)+S(k,l))*2.0)*2.0+S(j,j)+S(k,k),
			(S(i,i)+S(j,j)+S(k,k)+S(l,l)+S(i,k)+S(j,l)-(S(i,j)+S(j,k)+S(k,l))*2.0)*2.0+S(j,j)+S(k,k), resol);
}


inline bool check_equiv_m(const SymMat<Double>& lhs, //const Double& lhs_detS, const Double& lhs_detS_var,
		const SymMat<Double>& rhs, //const Double& rhs_detS, const Double& rhs_detS_var,
		const Double& resol)
{
	assert( lhs.size() == 3 );
	assert( rhs.size() == 3 );

	if( !equiv_resol( lhs(0,0), rhs(0,0), resol) )
	{
		return false;
	}
	if( !equiv_resol( lhs(1,1), rhs(1,1), resol) )
	{
		return false;
	}
	if( !equiv_resol( lhs(2,2), rhs(2,2), resol) )
	{
		return false;
	}
	if( !equiv_resol( lhs(0,0)+lhs(1,1)+lhs(0,1)*2.0, rhs(0,0)+rhs(1,1)+rhs(0,1)*2.0, resol) )
	{
		return false;
	}
	if( !equiv_resol( lhs(0,0)+lhs(2,2)+lhs(0,2)*2.0, rhs(0,0)+rhs(2,2)+rhs(0,2)*2.0, resol) )
	{
		return false;
	}
	if( !equiv_resol( lhs(1,1)+lhs(2,2)+lhs(1,2)*2.0, rhs(1,1)+rhs(2,2)+rhs(1,2)*2.0, resol) )
	{
		return false;
	}
	return true;
}


inline bool check_equiv_s(const SymMat<Double>& lhs, const SymMat<Double>& rhs,
		const Double& resol)
{
	assert( rhs.size() == lhs.size() );

	const Int4 ISIZE = lhs.size();
	for(Int4 i=0; i<ISIZE; i++)
	{
		if( !equiv_resol( lhs(i,i), rhs(i,i), resol) )
		{
			return false;
		}
		for(Int4 j=1; j<i; j++)
		{
			if( !equiv_resol( lhs(i,i)+lhs(j,j)+lhs(i,j)*2.0, rhs(i,i)+rhs(j,j)+rhs(i,j)*2.0, resol) )
			{
				return false;
			}
		}
	}

	return true;
}


#endif
