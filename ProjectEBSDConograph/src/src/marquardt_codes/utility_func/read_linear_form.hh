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
#ifndef _READ_LINEAR_FORM_HH_
#define _READ_LINEAR_FORM_HH_

#include <map>
#include "zmath.hh"
#include "zstring.hh"
#include "../zerror_type/error_mes.hh"

using namespace std;

template<class T>
static ZErrorMessage divide_term(const string& rhs_str, pair<string, T>& coef_term)
{
	if( rhs_str.empty() ) ZErrorMessage("EMPTY LINEAR FORM", __FILE__, __LINE__, __FUNCTION__);

	vector<string> term_tray, term_tray2;
	if( rhs_str.at(0) == '-' )
	{
		coef_term.second = -1;
		split(rhs_str.substr(1), term_tray, '*');
	}
	else
	{
		coef_term.second = 1;
		split(rhs_str, term_tray, '*');
	}

	Double coef2;
	bool flag = false;
	for(size_t l=0; l<term_tray.size(); l++)
	{
		split(term_tray[l], term_tray2, '/');
		if( term_tray2.empty() ) return ZErrorMessage("UNPROPER LINEAR FORM: "+ rhs_str, __FILE__, __LINE__, __FUNCTION__);

		if( isdigit( term_tray2[0].at(0) ) )
		{
			if( !str2num(term_tray2[0], coef2) ) return ZErrorMessage("UNPROPER LINEAR FORM: "+ rhs_str, __FILE__, __LINE__, __FUNCTION__);
			coef_term.second *= coef2;
		}
		else
		{
			if( flag ) coef_term.first += "*"+term_tray2[0];
			else coef_term.first = term_tray2[0];
			flag = true;
		}

		for(size_t l2=1; l2<term_tray2.size(); l2++)
		{
			if( term_tray2[l2].empty() ) return ZErrorMessage("UNPROPER LINEAR FORM: "+ rhs_str, __FILE__, __LINE__, __FUNCTION__);
			if( isdigit( term_tray2[l2].at(0) ) )
			{
				if( !str2num(term_tray2[l2], coef2) ) return ZErrorMessage("UNPROPER LINEAR FORM: "+ rhs_str, __FILE__, __LINE__, __FUNCTION__);
				coef_term.second /= coef2;
			}
			else
			{
				return ZErrorMessage("DIVISION BY A VARIABLE IS NOT ALLOWED IN LINEAR FORMS: "+ rhs_str, __FILE__, __LINE__, __FUNCTION__);
			}
		}
	}
	if( !flag ) coef_term.first = ""; // return as a constant term.
    return ZErrorMessage();
}


// On input, the first argument should be a linear form like "a*x*y + b*y + z/2 + 2/3".
// (Note: this function is not completely made so that it can handle all the polynomials.)
// In the above case, the second argument consists of the entries <"x*y", a>, <"y", b>, <"z", 1/2>, <"", 2/3> on output.
template<class T>
ZErrorMessage read_linear_form(const string& arg, map<string, T>& coef_term)
{
	vector<string> rhs_str;
	istringstream iss(remove_blank(arg));
	if( !split_term(iss, rhs_str) ) throw nerror_arg( arg, __FILE__, __LINE__, __FUNCTION__);

	pair<string, T> term;
	ZErrorMessage zerr;
	coef_term.clear();
	for(vector<string>::const_iterator it=rhs_str.begin(); it!=rhs_str.end(); it++)
	{
		zerr = divide_term(*it, term);
		if( ZErrorNoError != zerr.putErrorType() ) return zerr;
		typename map<string, T>::iterator it_map = coef_term.find(term.first);
		if( it_map == coef_term.end() ) coef_term.insert(term);
		else
		{
			it_map->second +=term.second;
		}
	}

	return ZErrorMessage();
}



#endif /* _READ_LINEAR_FORM_HH_ */
