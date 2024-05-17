/*
 * The MIT License

   EBSD-Conograph (Conograph software for EBSD ab-initio indexing)

Copyright (c) <2020> <Ryoko Oishi-Tomiyasu, Kyushu University>

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
#ifndef ZMATH_HH_
#define ZMATH_HH_

#include"../ZAnalysisTypes.hh"
#include"../HC_algorithm/GCD_1_3_1.hh"

template<class T>
inline T iround_half_up(const double& value)
{
	if (value >= 0.0)
   		return T(floor(+value + 0.5));
	else
		return - T(floor(-value + 0.5));
}

template<class T>
inline T ifloor(const double& value)
{
	return T(floor(value));
}

template<class T>
inline T iceil(const double& value)
{
	return T(ceil(value));
}

template<class T>
inline bool dbl2fraction(const Double& dbl, pair<T, T>& frac)
{
	const T k = iround_half_up<T>(dbl*48);
	if( fabs( k - dbl*48 ) >= 1.0e-8 ) return false;
	
	if( k == 0 )
	{
		frac.first = 0;
		frac.second = 1;
	}
	else
	{
		const T common_divisor = GCD_Euclid(k, 48);
		frac.first = k/common_divisor;
		frac.second = 48/common_divisor;
	}
	return true;
}

#endif /*GCD_HH_*/
