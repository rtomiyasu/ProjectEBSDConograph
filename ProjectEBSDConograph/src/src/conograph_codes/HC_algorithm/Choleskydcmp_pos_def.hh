/*
 Choleskydcmp.hh

 Author:  JAMA_CHOLESKY_H (modified by R. Tomiyasu)

 Copyright 2010 KEK. All rights reserved.
 */

#ifndef JAMA_CHOLESKY_H
#define JAMA_CHOLESKY_H

#include "../../conograph_codes/utility_data_structure/nrutil_nr.hh"

/**
	Constructs a lower triangular matrix L, such that L*L'= A.
	On input, if L_ is not symmetric positive-definite (SPD), only a
	partial factorization is performed.  If is_spd()
	evalutate true (1) then the factorizaiton was successful.
*/
template <class Real>
bool Choleskydcmp(NRMat<Real>& L_)
{
	assert(L_.nrows() == L_.ncols());

      // Main loop.
     for (int j = 0; j < L_.ncols(); j++)
	 {
        Real d(0.0);
        for (int k = 0; k < j; k++)
		{
            Real s(0.0);
            for (int i = 0; i < k; i++)
			{
               s += L_[k][i]*L_[j][i];
            }
            s = (L_[j][k] - s)/L_[k][k];
            d = d + s*s;

            L_[j][k] = s;
            
         }

         d = L_[j][j] - d;
         if( d > 0.0 )
         {
             L_[j][j] = sqrt(d);
         }
         else return false;

         for (int k = j+1; k < L_.ncols(); k++)
	     {
         	L_[j][k] = 0.0;
         }
	}
     return true;
}

/**

	Solve a linear system A*x = b, using the previously computed
	cholesky factorization of A: L*L'.

   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     x so that L*L'*x = b.  If b is nonconformat, or if A
   				was not symmetric posidtive definite, a null (0x0)
   						array is returned.
*/
template <class Real>
void solve(const NRMat<Real>& L_, const NRVec<Real> &b, NRVec<Real>& x, const size_t& ibegin, const size_t& iend)
{
	assert(b.size() == L_.nrows());
	assert(L_.ncols() == L_.nrows());
	assert(0 <= ibegin && iend <= L_.nrows());

	x = b;

    // Solve L*y = b;
    for (int k = ibegin; k < iend; k++)
	{
        if( L_[k][k] <= 0.0 )
        {
        	x[k] = 0.0;
        }
		else
		{
	        for (int i = ibegin; i < k; i++)
	        {
	        	x[k] -= x[i]*L_[k][i];
	        }
			x[k] /= L_[k][k];
		}
    }

    // Solve L'*X = Y;
    for (int k = iend-1; k >= ibegin; k--)
	{
        if( L_[k][k] <= 0.0 )
        {
        	x[k] = 0.0;
        }
        else
        {
            for (int i = k+1; i < iend; i++)
            {
            	x[k] -= x[i]*L_[i][k];
            }
        	x[k] /= L_[k][k];
        }
    }
}

/**

	Solve a linear system A*X = I, using the previously computed
	cholesky factorization of A: L*L'.

   @return     X so that L*L'*X = B.  If B is nonconformat, or if A
   				was not symmetric posidtive definite, a null (0x0)
   						array is returned.
*/
template <class Real>
void getInverseMatrix(const NRMat<Real>& L_,
		NRMat<Real>& InvL_)
{
	int n = L_.nrows();
	assert(L_.ncols() == n);
	assert(InvL_.nrows() == n);
	assert(InvL_.ncols() == n);

	InvL_ = 0.0;
	for (int j=0; j<n; j++) InvL_[j][j] = 1.0;
		    
	 // The lower triangle part of X equals L^{-1}.
  	 for (int k=0; k<n; k++)
	 {
  		if( L_[k][k] <= 0.0 )
  		{
  			InvL_[k][k] = 0.0;
  			continue;
  		}
	    for (int j=k; j<n; j++)
		{
	  		if( L_[j][j] <= 0.0 )
	  		{
	  			continue;
	  		}
			for (int i = k; i < j; i++) InvL_[j][k] -= L_[j][i]*InvL_[i][k];
			InvL_[j][k] /= L_[j][j];
		}
     }
}

#endif

// JAMA_CHOLESKY_H
