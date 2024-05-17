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

#ifndef _GCD_1_3_1_HH_
#define _GCD_1_3_1_HH_

template<class T>
T GCD_Euclid_For_NonNegative_Integers(T a,  T b)
{
	T r;
	while (true)
	{
		if (b == 0) return a;

		r = a % b;
		a = b;
		b = r;
	}
	return r;
}

template<class T>
T GCD_Euclid(T a,  T b)
{
	if( a < 0 ) a *= -1;
	if( b < 0 ) b *= -1;
	return GCD_Euclid_For_NonNegative_Integers(a, b);
}


//int main ()
//{
//	int a = 36;
//	int b = 24;
//
//	cout<<"y = "<<GCD_Euclid(a,b)<<endl;
//	return 0;
//}

#endif

