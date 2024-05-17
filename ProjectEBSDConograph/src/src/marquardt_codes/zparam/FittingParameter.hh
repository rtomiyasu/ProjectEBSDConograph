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

#ifndef FittingParameter_H_
#define FittingParameter_H_

#include "../../conograph_codes/ZAnalysisTypes.hh"

class FittingParameter 
{
public:
    double value; //Parameter's value
    double error; //Error value. (standard deviation)
	
	FittingParameter(){ value = 0.0; error = 0.0; };
	FittingParameter(const FittingParameter& rhs){ value = rhs.value; error = rhs.error; };
	FittingParameter(const Double& theValue){ value = theValue; error = 0.0; };
    FittingParameter(const Double& theValue, const Double& theError){ value = theValue; error = theError; };
	
    virtual ~FittingParameter(){};

    FittingParameter& operator=(const FittingParameter& rhs){ if(this != &rhs){ value = rhs.value; error = rhs.error; } return *this; };
    FittingParameter& operator=(const Double& t){ value = t; error = 0.0; return *this; };
	void Descriptor();
	inline FittingParameter& operator*=(const Double& rhs){ value *= rhs; error *= rhs; return *this; };
	inline bool operator<(const Double& rhs) const { return value < rhs; };
};

inline FittingParameter operator*(const FittingParameter& lhs, const Double& rhs)
{ 
	FittingParameter ans = lhs;
	ans *= rhs;
	return ans; 
}

#endif /*FittingParameter_H_*/
