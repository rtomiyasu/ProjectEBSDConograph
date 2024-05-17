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
#ifndef _ERROR_TYPE_H_
#define _ERROR_TYPE_H_

typedef enum {
	ZErrorNoError = 0, //  No error.
	
	//Standard Error
	ZErrorFailedMemoryAllocate,
	ZErrorArgmentSize,
	ZErrorArgument,
	ZErrorArrayOverFlow,
	ZErrorDelimiterNotFound,
	ZErrorDouplicatedLabel,
	ZErrorEscapeFunction,
	ZErrorFileNotFound,
	ZErrorFunctionIsNotImplemented,
	ZErrorInterrupted,
	ZErrorInvalidKey,
	ZErrorLabelNotFound,
	ZErrorManyIteration,
	ZErrorNullPointer,
	ZErrorNoValue,
	ZErrorOutRange,
	ZErrorZeroDivision,

	//General Error
	ZErrorFileFormatBroken,
	ZErrorNotIGORFile,
	ZErrorErrorsAreNotContained,
	ZErrorPeakPositionsAreNotContained,

	//Undefined
	ZErrorUndefined,

} ZErrorType;

#endif
