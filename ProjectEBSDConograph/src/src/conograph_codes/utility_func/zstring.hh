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
#ifndef _ZSTRING_HH_
#define _ZSTRING_HH_

#include <sstream>
#include "../ZAnalysisTypes.hh"

// Change the argument type from T to string.
template<class T>
inline bool str2num(const string& str, T& ans)
{
  istringstream iss(str);
  iss >> ans;
  if( iss.fail() ) return false;
  return true;
}


template<class T>
inline string num2str(const T& num)
{
  stringstream strstream;
  strstream << num;
  return strstream.str();
}


template<class T>
inline string num2str(const T& num, const Int4& precision)
{
  stringstream strstream;
  strstream.precision(precision);
  strstream << num;
  return strstream.str();
}

template<class T>
inline string vec2str(const vector<T>& num)
{
	stringstream strstream;
	if( num.empty() ) return "";
	if( num.size() < 2 ) return num2str(num[0]);

	strstream << "(" << *num.begin();
	for(typename vector<T>::const_iterator it=num.begin()+1; it<num.end(); it++)
	{
		strstream << "," << *it;
	}
	strstream << ")";
	return strstream.str();
}

// Check if the string is blank.
bool is_blank(const string&);

// Split the argument string by the delimiter.
void split(const string&, vector<string>&, const char&);

bool split_term(istream&, vector<string>&);

string getFileExtension (const string& inputFilePathString);

void removeFileExtension(const string&, string&);

void removeFilePath(const string&, string&);

inline string remove_blank(const string& str)
{
	string ans, str2;
	istringstream iss(str);
	
	iss >> str2;
	while( !iss.fail() )
	{
		ans += str2;
		iss >> str2;
	}
	return ans;
}

inline bool getnewline(istream& ifs, string& s)
{
	s.clear();
	char c;
	while( ifs.get(c) )
	{
		if( c == '\n' ) return true;
		if( c == '\r' )
		{
			if( ifs.eof() ) return true;
			c = ifs.peek();
			if( c == '\n' ) ifs.get(c);
			return true;
		}
		s += c;
	}
	return false;
}


#endif /*STRING_HH_*/
