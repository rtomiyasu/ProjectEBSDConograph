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
#include "error_mes.hh"

ZErrorMessage::ZErrorMessage(const ZErrorType& arg, const string& filename, const Int4& line_number, const string& funcname)
{
	string err_message;
	switch (arg) {
		case ZErrorNoError:{ err_message = "NO ERROR"; }break;
		case ZErrorFailedMemoryAllocate:{ err_message = "CANNOT ALLOCATE RAM"; }break;
		case ZErrorArgmentSize:{ err_message = "THE SIZE OF THE ARGUMENTS ARE NOT PROPER"; }break;
		case ZErrorArgument:{ err_message = "THE ARGUMENTS ARE NOT PROPER"; }break;
		case ZErrorArrayOverFlow:{ err_message = "ARRAY OVERFLOW HAS HAPPENED"; }break;
		case ZErrorEscapeFunction:{ err_message = "ESCAPE FROM THE NON OVERLOADED FUNCTION"; }break;
		case ZErrorFileNotFound:{ err_message = "FILE IS NOT FOUND"; }break;
		case ZErrorFunctionIsNotImplemented:{ err_message = "SOLLY THIS FUNCTION IS NOT IMPLEMENTED"; }break;
		case ZErrorInterrupted:{ err_message = "THE PROGRAM WAS INTERRUPTED"; }break;
		case ZErrorInvalidKey:{ err_message = "INVALID KEY IN FILE"; }break;
		case ZErrorManyIteration:{ err_message = "TOO MANY ITERATIONS"; }break;
		case ZErrorNoValue:{ err_message = "NO VALUE"; }break;
		case ZErrorNullPointer:{ err_message = "NULL POINTER"; }break;
		case ZErrorOutRange:{ err_message = "SOME PARAMETER IS OUT OF RANGE"; }break;
		case ZErrorFileFormatBroken:{ err_message = "FILE FORMAT IS BROKEN"; }break;
		case ZErrorZeroDivision:{ err_message = "ZERO DIVISION OCCURED HERE"; }break;
		case ZErrorUndefined:
		default:{ err_message = "SOME ERROR HAS OCCURED. FOR DETAILS, PLEASE REFER THE LOG MESSAGES"; } break;
	}

	setErrorMessage(arg, err_message, filename, line_number, funcname);
}


string ZErrorMessage::printErrorLog() const
{
   	stringstream os;
	os << putSrcFileName() + "(line " << putSrcLineNumber() << ") : " + putSrcFuncName() + " fails.\n"
		+ putErrorMessage() + ".\n";
	return os.str();
}
