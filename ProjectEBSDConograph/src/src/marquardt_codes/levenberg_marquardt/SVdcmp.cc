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
#include <cmath>
#include <algorithm>
#include "SVdcmp.hh"
#include "../zerror_type/error_out.hh"
#include "../../conograph_codes/utility_data_structure/index_set.hh"

using namespace std;

/*
	Decompose a symmetric matrix alpha into transpose(V) * T * V using matrices as followings.
	T is a diagonal matrix. V is a matrix such that	
	transpose(V) * V is a diagonal matrix whose j'th diagonal element equals to the input alpha[j][j].  
    On Output, alpha is replaced by V^-1.
    adiag_in[j] equal the j'th diagonal element of T.
*/   
void SVdcmp_precondition(NRMat<Double>& alpha, NRVec<Double>& adiag_in, const Int4& ibegin, const Int4& iend)
{
    const Int4 asize = adiag_in.size();

	assert( alpha.nrows() == asize && alpha.ncols() == asize );
	assert( 0 <= ibegin && iend <= asize );

	// Set adiag_out as diagonal elements of a diagonal matrix D.
	NRVec<Double> adiag_out(asize);
    for (Int4 j=ibegin; j<iend; j++)
    {
		if( alpha[j][j] <= 0.0 ) adiag_out[j] = 0.0;
		else adiag_out[j] = 1.0/sqrt(alpha[j][j]);
    }
 
	// Decompose alpha into D^-1 * M * D^-1.
	// M is a matrix whose diagonal elements all equal 1.
    // alpha is replace by M.
	for (Int4 j=ibegin; j<iend; j++)
	{
		if(adiag_out[j] <= 0.0)
		{
			for (Int4 k=ibegin; k<j; k++) alpha[j][k] = 0.0;
		    alpha[j][j] = 0.0;
		}
		else
		{
 			for (Int4 k=ibegin; k<j; k++) alpha[j][k] *= adiag_out[j]*adiag_out[k];
		    alpha[j][j] = 1.0;
		}
    }
    for (Int4 j=ibegin; j<iend; j++)
	    for (Int4 k=j+1; k<iend; k++) alpha[j][k] = alpha[k][j];

    // Singular Value Decomposition ( alpha = U * adiag_in * Transpose(U) ).
    // Replace alpha by U. 
   	SVdcmp(alpha, adiag_in, ibegin, iend);

   	// alpha is replace by D*U.
   	for (Int4 j=ibegin; j<iend; j++)
   		for (Int4 k=ibegin; k<iend; k++) alpha[j][k] *= adiag_out[j];

}


// Calculate the inverse matrix of transpose(V^-1) * adiag * V^-1. (adiag is a diagonal matrix.)
void getInverseMatrix(const NRMat<Double>& V, const NRVec<Double>& adiag, 
		const Int4& ibegin, const Int4& iend,
		const Double& thred, SymMat<Double>& InvMat)
{
	assert( V.nrows() == adiag.size() && V.ncols() == adiag.size() );
	assert( 0 <= ibegin && iend <= adiag.size() );
	assert( InvMat.size() + ibegin == iend );

	InvMat = 0.0;

	Double t;
	for(int l=ibegin; l<iend; l++)
	{
		if( adiag[l] <= thred ) continue;

		t = 1.0/adiag[l];
		for (int j=ibegin,j2=0; j<iend; j++,j2++)
			for (int k=ibegin,k2=0; k<j+1; k++,k2++) InvMat(j2,k2)+=V[j][l]*V[k][l]*t;
	}
}


/*
    alpha is decomposed into U * T * Transpose(U) such that
    	U:Orthogonal matrix
    	d: Diagonal elements of the diagonal matrix T.
    On Output, alpha is replaced by U,
*/   
void SVdcmp(NRMat<Double>& alpha, NRVec<Double>& d, const Int4& ibegin, const Int4& iend)
{
	if( ibegin >= iend ) return;
	const Int4 asize = d.size();

	assert( alpha.nrows() == asize && alpha.ncols() == asize );
	assert( 0 <= ibegin && iend <= asize );

	
	NRVec<Double> e(asize);

    // Householeder method.
    tred2(alpha, d, e, ibegin, iend);

    // QL method.
    tqli(alpha, d, e, ibegin, iend);
}

/*
    Dicompose the symmetric matrix alpha into U * T * Transpose(U) such that
    U:Orthogonal matrix,
    d: Diagonal elements of the tridiagonal matrix T,
    e: Offdiagonal elements of the tridiagonal matrix T with e[0] = 0.
    The return value of alpha is U.
*/   
void tred2(NRMat<Double>& alpha, NRVec<Double>& d, NRVec<Double>& e, const Int4& ibegin, const Int4& iend)
{
	assert( e.size() == d.size() );
	assert( alpha.nrows() == d.size() );
	assert( alpha.ncols() == d.size() );
	assert( 0 <= ibegin && ibegin < iend && iend <= d.size() );

	int l;
    Double scale,hh,h,g,f;

//	int asize = d.size();

	for (int i=iend-1;i>ibegin;i--)
	{
		l=i-1;
		h=0.0;
		scale=0.0;
		if (l > ibegin) {
			for (int k=ibegin;k<l+1;k++) scale += fabs(alpha[i][k]);
			
			if (scale == 0.0) e[i]=alpha[i][l];
			else {
				for (int k=ibegin;k<l+1;k++) {
					alpha[i][k] /= scale;
					h += alpha[i][k]*alpha[i][k];
				}
				f=alpha[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				alpha[i][l]=f-g;
				f=0.0;
				for (int j=ibegin;j<l+1;j++) {
					alpha[j][i]=alpha[i][j]/h;
					g=0.0;
					for (int k=ibegin;k<j+1;k++)
						g += alpha[j][k]*alpha[i][k];
					for (int k=j+1;k<l+1;k++)
						g += alpha[k][j]*alpha[i][k];
					e[j]=g/h;
					f += e[j]*alpha[i][j];
				}
				hh=f/(h+h);

				for (int j=ibegin;j<l+1;j++) {
					f=alpha[i][j];
					g = e[j]-hh*f;
					e[j] = g;
					for (int k=ibegin;k<j+1;k++)
						alpha[j][k] -= (f*e[k]+g*alpha[i][k]);
				}
			}
		}
		else e[i]=alpha[i][l];
		
		d[i]=h;
	}
	d[ibegin]=0.0;
	e[ibegin]=0.0;
	for (int i=ibegin;i<iend;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (int j=ibegin;j<l;j++) {
				g=0.0;
				for (int k=ibegin;k<l;k++)
					g += alpha[i][k]*alpha[k][j];
				for (int k=ibegin;k<l;k++)
					alpha[k][j] -= g*alpha[k][i];
			}
		}
		d[i]=alpha[i][i];
		alpha[i][i]=1.0;
		for (int j=ibegin;j<l;j++){
			alpha[j][i]=0.0;
			alpha[i][j]=0.0;
		}
	}
}

/*
    This function dicompose the tridiagonal symmetric matrix T into U * D * Transpose(U) such that
    U:Orthogonal matrix,
    D:Diagonal matrix D.

    On input, d contains the diagonal elements of T.
    e contains the offdiagonal elements of T with e[0] = 0.
    On output, d contains the diagonal elements of D.
    Ort is multiplied by U on the right side.
*/   
void tqli(NRMat<Double>& z, NRVec<Double>& d, NRVec<Double>& e, const Int4& ibegin, const Int4& iend)
{
	assert( e.size() == d.size() );
	assert( z.nrows() == d.size() );
	assert( z.ncols() == d.size() );
	assert( 0 <= ibegin && iend <= d.size() );

	int m,l,iter,i,k;
	Double s,r,p,g,f,dd,c,b;

	const Int4 max_it = max(30, iend-ibegin);
	
	for (i=ibegin+1;i<iend;i++) e[i-1]=e[i];
	e[iend-1]=0.0;
	for (l=ibegin;l<iend;l++)
	{
		iter=0;
		do {
			for (m=l;m<iend-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if(m != l){
				if ( iter++ == max_it)
				{
					throw ZErrorMessage(ZErrorManyIteration, __FILE__, __LINE__, __FUNCTION__);
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=ibegin;k<iend;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}


void putMatrixDiagonalAscendingOrder(const NRMat<Double>& alpha, Vec_INT& row_order)
{
    const int asize = alpha.nrows();
    assert( alpha.ncols() == asize );
    
    vector< index_set<Double> > tray(asize);
    for(Int4 k=0; k<asize; k++){
    	tray[k].index = k;
    	tray[k].element = fabs(alpha[k][k]);
    }

	// Sort elements into ascending order.
	stable_sort(tray.begin(), tray.end());
	
    row_order.resize(asize);
    for(Int4 k=0; k<asize; k++) row_order[k] = tray[k].index;
}


void permSquareMatrixElement(const Vec_INT& row_order, NRMat<Double>& alpha)
{
    const int isize = row_order.size();
	
	Vec_BOOL vflag(isize+1, true);	// The last index is against array overflow.
	
	Int4 row_start = 0, i, i2;
	while(row_start < isize){
		i = row_start;
		i2 = row_order[i];
		vflag[i] = false;
		if(i != i2)
			while( vflag[i2] ){
				for(Int4 j=0; j<isize; j++) swap(alpha[i][j], alpha[i2][j]);
				for(Int4 j=0; j<isize; j++) swap(alpha[j][i], alpha[j][i2]);
				vflag[i2] = false;
				i = row_order[i];
				i2 = row_order[i];
			}

		row_start++;
		while( !vflag[row_start] ) row_start++;
	}
}


void permSquareMatrixRow(const Vec_INT& row_order, NRMat<Double>& alpha)
{
    const int isize = row_order.size();

	Vec_BOOL vflag(isize+1, true);	// The last index is against array overflow.
	
	Int4 row_start = 0, i, i2;
	while(row_start < isize)
	{
		i = row_start;
		i2 = row_order[i];
		vflag[i] = false;
		if(i != i2)
			while( vflag[i2] )
			{
				for(Int4 j=0; j<isize; j++) swap(alpha[i2][j], alpha[row_start][j]);
				vflag[i2] = false;
				i = row_order[i];
				i2 = row_order[i];
			}

		row_start++;
		while( !vflag[row_start] ) row_start++;
	}
}
