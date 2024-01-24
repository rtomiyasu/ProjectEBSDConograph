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
#ifndef ENUMAXIS_HH_
#define ENUMAXIS_HH_

#include <map>
#include "../ZAnalysisTypes.hh"

enum eABCaxis{ UndefinedABC_Axis=-1, A_Axis=0, B_Axis=1, C_Axis=2 };
enum eRHaxis{ UndefinedRho_Axis=-1, Rho_Axis=0, Hex_Axis=1 };

inline map<eABCaxis, string> putABCaxisLabel()
{
	map<eABCaxis, string> ans;
	ans[A_Axis] = "A";
	ans[B_Axis] = "B";
	ans[C_Axis] = "C";
	return ans;
}

inline map<eRHaxis, string> putRHaxisLabel()
{
	map<eRHaxis, string> ans;
	ans[Rho_Axis] = "Rhombohedral";
	ans[Hex_Axis] = "Hexagonal";
	return ans;
}

template<class T>
inline string find_data(const map<T, string>& tray, const T& arg)
{
	typename map<T, string>::const_iterator it=tray.find(arg);
	if( it == tray.end() ) return "Undefined";
	else return it->second;
}

inline string str_axis(const eABCaxis& k)
{
	static map<eABCaxis, string> ABCaxis_Str = putABCaxisLabel();
	return find_data(ABCaxis_Str, k);
}

inline string str_axis(const eRHaxis& k)
{
	static map<eRHaxis, string> RHaxis_Str = putRHaxisLabel();
	return find_data(RHaxis_Str, k);
}

// non-principal axis.
inline Int4 npaxis(const Int4& b, const eABCaxis& i)
{
	static const Int4 npaxis[][3] = { {1,2,0}, {2,0,1} };
	return npaxis[b][(size_t)i];
}

#endif /*ENUMAXIS_HH_*/
