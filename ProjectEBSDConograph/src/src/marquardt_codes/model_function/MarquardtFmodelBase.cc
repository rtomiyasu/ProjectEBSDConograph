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
#include "MarquardtFmodelBase.hh"

void calculate_dparam(const NRVec<Double>& adiag, const NRVec<Double>& beta, NRVec<Double>& dbeta,
		const Double& alamda, const Double& eps1, const Int4& ibegin, const Int4& iend)
{
	assert( 0 <= ibegin && iend <= adiag.size() );
	
   	// dbeta[j] = beta[j] / (adiag[j] + alamda).
	Double s;
   	for (int j=ibegin; j<iend; j++)
	{
   		if( adiag[j] <= eps1 )
   		{
   	   		dbeta[j] = 0.0;
   			continue;
   		}
//   		s = adiag[j];
//   		if( j < num_indep_nonlinear_param_all )
   		s = adiag[j] + alamda;
//   		if( s <= eps1 ) continue;
   		
   		dbeta[j] = beta[j] / s;
	}
}


void check_diagonal(const NRVec<Double>& adiag, const Double& eps1)
{
	const Int4 num_indep_param = adiag.size();
	
   	// dbeta[j] = beta[j] / (adiag[j] + alamda).
	Double mind = 1.0, maxd = 1.0;
	Vec_DP tray;
   	for (int j=0; j<num_indep_param; j++)
	{
   		if( adiag[j] <= eps1 )
   		{
   			Vec_DP::iterator it = upper_bound(tray.begin(), tray.end(), adiag[j]);
   	   		tray.insert(it, adiag[j]);
   		}
   		else if( adiag[j] < mind ) mind = adiag[j]; 
   		else if( adiag[j] > maxd ) maxd = adiag[j]; 
	}
//#ifdef DEBUG
//   	cout << "The maximum eigenvalue d_max : " << maxd << endl;
//   	cout << "The minimum eigenvalue d_min > Neps = " << eps1 << " : " << mind << endl;
//   	if( !tray.empty() )
//   	{
//   		cout << "The eigenvalues d < Neps : ";
//   		for(Vec_DP::const_iterator it=tray.begin(); it!=tray.end(); it++) cout << *it << " ";
//   		cout << endl;
//   	}
//   	cout << "The condition number d_max/d_min : " << maxd / mind << endl;
//   	cout << endl;
//#endif 
}
