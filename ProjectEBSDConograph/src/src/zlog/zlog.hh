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
#ifndef ZLOG_HH_
#define ZLOG_HH_

#include <fstream>

#include "../zlog/rlog.hh"

class CCoutListner : public CRLog::Listner
{
public:
	CCoutListner(){};
	virtual ~CCoutListner(){};

	void log(const zLevel& cLevel, const string& cpStr);
	void log(const zLevel& cLevel, const int& threadID, const string& cpStr) { log(cLevel, cpStr); };
   	zListnerID putId(void) const{ return zListnerID(0); };
};

class FileoutListner : public CRLog::Listner
{
private:
	const zListnerID m_ID;
	ofstream* ofs;

public:
	FileoutListner(const string& fname, const zListnerID& ID)
		: m_ID(ID), ofs(new ofstream(fname.c_str())){};
	virtual ~FileoutListner() { if( ofs != NULL ){ ofs->close(); delete ofs; } };

	void log(const zLevel& cLevel, const string& cpStr)
    {
        *ofs << cpStr;
    };
	void log(const zLevel& cLevel, const int& threadID, const string& cpStr) { log(cLevel, cpStr); };

   	zListnerID putId(void) const{ return m_ID; };
};

#endif /* ZLOG_HH_ */
