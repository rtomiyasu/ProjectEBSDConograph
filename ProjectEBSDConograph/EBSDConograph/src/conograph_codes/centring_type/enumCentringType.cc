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

#ifdef DEBUG
	#include <iostream>
#endif
#include "enumCentringType.hh"
using namespace std;

string CenteringType::ToString(const CenteringType::Enum e)
{
	static const size_t NUM_C = 7;
	static const string name[NUM_C] = { "(P)", "(A)", "(B)",
											"(C)", "(F)",
											"(I)", "(hexagonal-axis)" };

	if( 0 <= (size_t)e && (size_t)e < NUM_C ) return name[(size_t)e];
	return "(unknown)";
}
