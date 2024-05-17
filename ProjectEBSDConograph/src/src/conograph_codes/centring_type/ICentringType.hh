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
#ifndef ICENTRINGTYPE_HH_
#define ICENTRINGTYPE_HH_

#include"../ZAnalysisTypes.hh"
#include"../symmetric_operation/S1.hh"
#include"../utility_data_structure/VecDat3.hh"
#include"enumCentringType.hh"

typedef VecDat3<Int4> MillerIndex3;

// Class for centring types.
class ICentringType
{
public:
	virtual const string& Name() const = 0;
	virtual const VecDat3<S1>* putTranslationVector(Int4&) const = 0;
	virtual Int4 putSTFactorBlavaisLatticeCoef(const MillerIndex3&) const = 0;
	virtual ~ICentringType(){};

	template<class T>
	Int4 putEquivalentPosition(const VecDat3<T>&, vector< VecDat3<T> >&) const;
};

template<class T>
Int4 ICentringType::putEquivalentPosition(const VecDat3<T>& xyz, vector< VecDat3<T> >& xyz_equiv) const
{
	Int4 num; 
	const VecDat3<S1>* trans_vec = putTranslationVector(num);

	xyz_equiv.resize(num);
	for(Int4 k=0; k<num; k++) xyz_equiv[k] = xyz + trans_vec[k];
	return num; 
}

#endif /*ICentringType_H_*/
