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
#ifndef _ERROR_MESS_HH_
#define _ERROR_MESS_HH_

// error_out.hh
#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "ZErrorType.hh"

class ZErrorMessage
{
private:
	ZErrorType Type;
	Int4 SRC_LINE_NUMBER;
	string SRC_FILE_NAME;
	string SRC_FUNC_NAME;
	string ERR_MASSAGE;

	template<class T>
	inline void setErrorMessage(const T& mess);
	template<class T>
	inline void setErrorMessage(const ZErrorType&,
								const T& mess,
								const string& filename, const Int4& line_number, const string& funcname);

public:
	ZErrorMessage(){ Type = ZErrorNoError; SRC_LINE_NUMBER = 0; };

	template<class T>
	ZErrorMessage(const ZErrorType& arg,
					const T& mess,
					const string& filename,
					const Int4& line_number,
					const string& funcname) { setErrorMessage(arg, mess, filename, line_number, funcname); };

	ZErrorMessage(const ZErrorType&, const string& filename, const Int4& line_number, const string& funcname);
	template <class T>
	ZErrorMessage(const T& mess, const string& filename, const Int4& line_number, const string& funcname);

	virtual ~ZErrorMessage(){};

	inline void setErrorType(const ZErrorType& arg){ Type = arg; };
	inline void setErrorMessage(const string& arg){ ERR_MASSAGE = arg; };
	inline void setSrcLineNumber(const Int4& line_number);
	inline void setSrcFileName(const string& filename);
	inline void setSrcFuncName(const string& funcname);

	inline const ZErrorType& putErrorType() const { return Type; };
	inline const Int4& putSrcLineNumber() const;
	inline const string& putSrcFileName() const;
	inline const string& putSrcFuncName() const;
	inline const string& putErrorMessage() const;

	string printErrorLog() const;
};

template<class T>
ZErrorMessage::ZErrorMessage(const T& mess, const string& filename, const Int4& line_number, const string& funcname)
{
	setErrorMessage(ZErrorUndefined, mess, filename, line_number, funcname);
}

inline void ZErrorMessage::setSrcLineNumber(const Int4& line_number)
{
	SRC_LINE_NUMBER = line_number;
}

inline void ZErrorMessage::setSrcFileName(const string& filename)
{
	SRC_FILE_NAME = filename;
}

inline void ZErrorMessage::setSrcFuncName(const string& funcname)
{
	SRC_FUNC_NAME = funcname;
}

template<class T>
inline void ZErrorMessage::setErrorMessage(const T& mess)
{
	ERR_MASSAGE = string(mess);
}

template<class T>
inline void ZErrorMessage::setErrorMessage(const ZErrorType& arg, const T& mess,
		const string& filename, const Int4& line_number, const string& funcname)
{
	setErrorType(arg);
	setErrorMessage(mess);
	setSrcLineNumber(line_number);
	setSrcFileName(filename);
	setSrcFuncName(funcname);
}

inline const Int4& ZErrorMessage::putSrcLineNumber() const
{
	return SRC_LINE_NUMBER;
}

inline const string& ZErrorMessage::putSrcFileName() const
{
	return SRC_FILE_NAME;
}

inline const string& ZErrorMessage::putSrcFuncName() const
{
	return SRC_FUNC_NAME;
}

inline const string& ZErrorMessage::putErrorMessage() const
{
	return ERR_MASSAGE;
}

class ZErrorMessageReadingFile : public ZErrorMessage
{
private:
	string READING_FILE_NAME;

public:
	ZErrorMessageReadingFile() : ZErrorMessage() {};
	ZErrorMessageReadingFile(const string& filename, const ZErrorMessage& arg) : ZErrorMessage(arg) { READING_FILE_NAME = filename; this->setErrorMessage("Program failed to read the file : " + filename + ".\n" + putErrorMessage()); };

	inline void setReadingFileName(const string& fname);
	inline const string& putReadingFileName() const;
};

inline void ZErrorMessageReadingFile::setReadingFileName(const string& fname)
{
	setErrorType(ZErrorNoError);
	READING_FILE_NAME = fname;
}


const string& ZErrorMessageReadingFile::putReadingFileName() const
{
	return READING_FILE_NAME;
}

#endif
