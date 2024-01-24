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

#ifndef _POWER_1_2_3_HH_
#define _POWER_1_2_3_HH_

static int int_bit(int n)
{
	int e = 0;
	n >>= 1;
	while( n )
	{
		e++;
		n >>= 1;
	}
	return e;
}

// Algorithm 1.2.3 page#9
template<class T>
static T power_LRB_bit(const T& g, const int& n)
{
	if (n==0) return T(1);

	int N;
	double z;
	if (n<0)
	{
		z=T(1)/g;
		N=-n;
	}
	else{
		z=g;
		N=n;
	}

	//2^f <= |n| < 2^(f+1)
	int f=int_bit(N);
	T y=z;

	while (f!=0)
	{
		f-=1;
		y*=y;
		if ( (N >> f) & 1 ) //bit number f of N
		{
			y*=z;
		}
	}
	return y;
}

//int main ()
//{
//	int n=-7;
//	double g=2;
//
//	cout<<"y = "<<power_LRB_bit(g,n)<<endl;
//	return 0;
//}

#endif

