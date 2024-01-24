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
#ifndef SYMMATNPLUS1N_HH_
#define SYMMATNPLUS1N_HH_

#include <map>
#include "../utility_data_structure/SymMat.hh"

typedef pair< SymMat<double>, NRMat<int> > SymMatNplus1N_Double;
// Assume the first element is a symmetric matrix of size 3 and the second element is a 4-by-3 matrix.
typedef pair< SymMat<double>, NRMat<int> > SymMat43_Double;

inline bool operator<(const SymMat<double>& lhs, const SymMat<double>& rhs)
{
	assert( lhs.size() == rhs.size() );
	const int ISIZE = lhs.size();

	for(int i=0; i<ISIZE; i++)
		for(int j=i; j<ISIZE; j++)
		{
			if( lhs(i,j) < rhs(i,j) ) return true;
			if( lhs(i,j) > rhs(i,j) ) return false;
		}
	return false;
};

#endif /*SYMMATNPLUS1N_HH_*/
