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
#ifndef _ERROR_OUT_HH_
#define _ERROR_OUT_HH_

// error_out.hh

#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../../conograph_codes/utility_func/zstring.hh"
#include "ZErrorType.hh"
#include "error_mes.hh"


inline ZErrorMessage nerror(const ZErrorType& arg, const string& filename, const Int4& line_number, const string& funcname)
{
	return ZErrorMessage(arg, filename, line_number, funcname);
};

inline ZErrorMessage nerror(const ZErrorType& arg, const string& mess, const string& filename, const Int4& line_number, const string& funcname)
{
	return ZErrorMessage(arg, mess, filename, line_number, funcname);
};

template<class T>
inline ZErrorMessage nerror_arg(const T& c, const string& filename, const Int4& line_number, const string& funcname)
{
	return ZErrorMessage(ZErrorArgument,
							"Wrong argument : "+num2str<T>(c), filename, line_number, funcname);
};

template<class T, class S>
inline ZErrorMessage nerror_out_range(const T& c1, const string& filename, const Int4& line_number, const string& funcname)
{
	return ZErrorMessage(ZErrorOutRange,
				num2str<T>(c1)+" is out of range", filename, line_number, funcname);
};


template<class T, class S>
inline ZErrorMessage nerror_out_range(const T& c1, const S& c2, const string& filename, const Int4& line_number, const string& funcname)
{
	return ZErrorMessage(ZErrorOutRange,
				num2str<T>(c1)+" is out of range : "+num2str<S>(c2), filename, line_number, funcname);
};


inline ZErrorMessage nerror(const bad_alloc&, const string& filename, const Int4& line_number, const string& funcname)
{
    return ZErrorMessage(ZErrorFailedMemoryAllocate, filename, line_number, funcname);
};

#endif
