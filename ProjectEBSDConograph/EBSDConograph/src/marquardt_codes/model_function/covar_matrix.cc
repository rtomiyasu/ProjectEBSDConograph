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
#include "covar_matrix.hh"

// Returns the inner product.(lhs, rhs are vectors.)
inline Double product_row(const Vec_DP_save& lhs, const SymMat<Double>& rhs, const Int4& row_index)
{
	Double t = 0.0;
	for(UInt4 l=0; l<lhs.size(); l++) t += rhs(row_index, lhs[l].index) * lhs[l].element;
	return t;
}


// Calculate the covariant matrix of the parameters from start'th to (end-1)'th
// using the covariant matrix of the parameters independently fit.
// ( Cmat * covar * transpose(Cmat) from start'th to (end-1)'th gives this matrix. )
void calPartOfCovariantMatrix(const Mat_DP_constr& Cmat,
const SymMat<Double>& covar, const Int4& start, const Int4& end, SymMat<Double>& covar_all)
{
	const Int4 ma2 = max(0, end - start);
    covar_all = SymMat<Double>(ma2, 0.0);

	// Calculate the covariant matrix which takes into account parameters 
    // that are being held fixed or dependent. (For fixed parameters, return zero covaiance.)
	Int4 start2=0;
	for (Int4 k=0; k<start; k++)
		if(Cmat[k].ID==_IDVary) start2++;

	const Int4 ISIZE = covar.size();
	
	Vec_DP tray(ISIZE);
	for (Int4 k=start, k0=0, k2=start2; k<end; k++, k0++)
	{
		if(Cmat.at(k).ID==_IDVary){
			for (Int4 j=start, j0=0, j2=start2; j<=k; j++, j0++){
				if(Cmat.at(j).ID==_IDVary) covar_all(k0,j0) = covar(k2,j2++);
				else if(Cmat.at(j).ID==_IDDepend) covar_all(k0,j0) = product_row( Cmat.at(j).constr, covar, k2 );
			}
			k2++;
		}
		else if(Cmat.at(k).ID==_IDDepend)
			for (Int4 j=start, j0=0, j2=start2; j<=k; j++, j0++){
				if(Cmat.at(j).ID==_IDVary){
					covar_all(k0,j0) = product_row( Cmat.at(k).constr, covar, j2++ );
				}
				else if(Cmat.at(j).ID==_IDDepend)
				{
					for(Int4 d=0; d<ISIZE; d++) tray[d] = product_row( Cmat.at(j).constr, covar, d );
					covar_all(k0,j0) = product_constr( Cmat.at(k).constr, &tray[0] );
				}
			}
	}
}



// Calculate the covariant matrix of the parameters of index in index_tray,
// using the covariant matrix of the parameters independently fit.
// ( Cmat * covar * transpose(Cmat) from start'th to (end-1)'th gives this matrix. )
void calPartOfCovariantMatrix(const Mat_DP_constr& Cmat,
const SymMat<Double>& covar, const Vec_INT& index_tray, SymMat<Double>& covar_all)
{
	
	const Int4 ma2 = index_tray.size();
    covar_all = SymMat<Double>(ma2, 0.0);
    if( ma2 <= 0 ) return;
    
	// Calculate the covariant matrix which takes into account parameters 
    // that are being held fixed or dependent. (For fixed parameters, return zero covaiance.)
    Vec_INT start2(ma2);
	for (Int4 k0=0, start=0, count=0; k0<ma2; k0++)
	{
		for(Int4 k=start; k<index_tray[k0]; k++)
		{
			if(Cmat[k].ID==_IDVary) count++;
		}
		start2[k0] = count;
		start = index_tray[k0];
	}

	const Int4 ISIZE = covar.size();
	
	Vec_DP tray(ISIZE);
	for (Int4 k0=0; k0<ma2; k0++)
	{
		const Int4& k = index_tray[k0];
		const Int4& k2 = start2[k0];
		if(Cmat.at(k).ID==_IDVary){
			for (Int4 j0=0; j0<=k0; j0++)
			{
				const Int4& j = index_tray[j0];
				const Int4& j2 = start2[j0];
				if(Cmat.at(j).ID==_IDVary)
				{
					covar_all(k0,j0) = covar(k2,j2);
				}
				else if(Cmat.at(j).ID==_IDDepend) covar_all(k0,j0) = product_row( Cmat.at(j).constr, covar, k2 );
			}
		}
		else if(Cmat.at(k).ID==_IDDepend)
			for (Int4 j0=0; j0<=k0; j0++)
			{
				const Int4& j = index_tray[j0];
				const Int4& j2 = start2[j0];
				if(Cmat.at(j).ID==_IDVary)
				{
					covar_all(k0,j0) = product_row( Cmat.at(k).constr, covar, j2 );
				}
				else if(Cmat.at(j).ID==_IDDepend)
				{
					for(Int4 d=0; d<ISIZE; d++) tray[d] = product_row( Cmat.at(j).constr, covar, d );
					covar_all(k0,j0) = product_constr( Cmat.at(k).constr, &tray[0] );
				}
			}
	}
}



// Copy the covariant matrix of the parameters from start'th to (end-1)'th.
void copyPartOfCovariantMatrix(const SymMat<Double>& covar, const Int4& start, const Int4& end, SymMat<Double>& ans)
{
	const Int4 ma2 = max(0, end - start);
  	ans = SymMat<Double>(ma2);

	for (Int4 k=start, k0=0; k<end; k++, k0++)
		for (Int4 j=start, j0=0; j<=k; j++, j0++) ans(k0,j0) = covar(k,j);
}


void copyPartOfCovariantMatrix(const SymMat<Double>& covar, const Vec_INT& index_tray, SymMat<Double>& ans)
{
	const Int4 ma2 = index_tray.size();
  	ans = SymMat<Double>(ma2);

	for (Int4 k0=0; k0<ma2; k0++)
		for (Int4 j0=0; j0<=k0; j0++) ans(k0,j0) = covar(index_tray[k0],index_tray[j0]);
}
