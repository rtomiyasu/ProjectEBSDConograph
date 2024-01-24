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
#include <limits>
#include "LemarqMethod.hh"

const Double LemarqMethod::EPS=numeric_limits<double>::epsilon();
const Double LemarqMethod::MAX_LAMBDA=1.0e+45;

LemarqMethod::LemarqMethod()
{
	m_limiter = 0.01;
    m_max_itnum = 1000;
	m_initial_lambda = 0.001;

	m_output_view_flag = false;
}

LemarqMethod::~LemarqMethod()
{
}

void LemarqMethod::setParam(
		const bool& output_view_flag,
		const Double& limiter, const Int4& itnum, const Double& initial_lambda)
{
	m_output_view_flag = output_view_flag;

	m_limiter = limiter;
	m_max_itnum = itnum;
	m_initial_lambda = initial_lambda;
}
