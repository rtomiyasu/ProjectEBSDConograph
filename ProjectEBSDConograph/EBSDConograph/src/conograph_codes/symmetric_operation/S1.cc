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
#include <cassert>
#include "../utility_func/zstring.hh"
#include "../utility_func/zmath.hh"
#include "S1.hh"

const Double S1::EPS = 1.0e-10;
const CMPX_DP S1::PI2i(0, 2.0 * M_PI);

S1::S1()
{
	mem = 1.0;
}

S1::S1(const S1& rhs)
{
	mem = rhs.mem;
}

S1::S1(const Double& dbl)
{
	mem = exp(PI2i*dbl);
}

S1::S1(const Int4& num, const Int4& den)
{
	static const CMPX_DP i(0.0, 1.0);
	static const CMPX_DP omega(-0.5, sqrt(0.75));

	if( num == 0 ) mem = 1.0;
	else if( 4 % den == 0 ){
		Int4 r = num * (4 / den) % 4;
		if(r < 1) mem = 1.0;
		else if(r < 2) mem = i;
		else if(r < 3) mem = -1.0;
		else mem = -i;
	}
	else if( 6 % den == 0 ){
		Int4 r = num * (6 / den) % 6;
		if(r < 1) mem = 1.0;
		else if(r < 2) mem = omega + 1.0;
		else if(r < 3) mem = omega;
		else if(r < 4) mem = -1.0;
		else if(r < 5) mem = conj(omega);
		else mem = conj(omega) + 1.0;
	}
	else mem=exp( PI2i*Double(num)/Double(den) );
}

S1::S1(const CMPX_DP&arg) {
	assert(arg != 0.0);
	mem = arg / abs(arg);
}



S1::~S1()
{
}

S1& S1::operator=(const S1& rhs)
{
	if(this != &rhs) mem = rhs.mem;
	return *this;
}



Double S1::toDouble() const
{
	static const Double PI2 = 2.0*M_PI;
	Double dbl = atan2( imag(mem), real(mem) ) / PI2;
	if(dbl < 0.0) dbl += 1.0;
	if(abs(dbl)<1e-10 || abs(1-dbl)<1e-10) { dbl = 0.0; } // "-0.0"‚â"1"‚ð"0"‚É‚·‚é

	return dbl;
}

string S1::toString() const
{
	const Double dbl = this->toDouble();
	pair<Int4, Int4> frac;
	if( !dbl2fraction(dbl, frac) ) return num2str<Double>(dbl);
	if( frac.first == 0 || frac.second == 1 ) return "0";
	return num2str<Int4>(frac.first)+"/"+num2str<Int4>(frac.second);
}
