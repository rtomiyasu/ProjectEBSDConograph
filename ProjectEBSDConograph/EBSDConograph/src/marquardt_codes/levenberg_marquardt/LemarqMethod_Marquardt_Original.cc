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
#include "LemarqMethod.hh"

// From the argument adiag, alpha beta and the initial parameters p1, ..., pm in paratry0,
// this method calculates dp1, ...., dpm.
// If the chi-squares of p1+dp1, ..., pm+dpm is smaller than the argument chisq0_all,
// new parameters are set in paratry0.
// Returns 0 : The argument paratry0 are replaced by new parameters.
// Returns 1 : The new parameters are removed and the old ones in argument paratry0 are reserved.
Int4 LemarqMethod::execute_Marquardt_Original(const Double& eps1,
const Double& chisq0_all, const Double& chisq_all,
Double& alamda) const
{
	if( chisq_all <= chisq0_all )
    {
		alamda *= 0.1;
		if( alamda < eps1 ) alamda = eps1;
   		return 0;
    }

	alamda = alamda*10.0;
   	return 1;
}
