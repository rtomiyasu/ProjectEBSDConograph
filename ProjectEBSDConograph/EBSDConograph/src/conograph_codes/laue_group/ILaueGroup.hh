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
#ifndef ILAUEGROUP_HH_
#define ILAUEGROUP_HH_

#include "../ZAnalysisTypes.hh"
#include "../utility_data_structure/index_set.hh"
#include "../utility_data_structure/VecDat3.hh"
#include "../point_group/enumPointGroup.hh"
#include "../symmetric_operation/MillerIndex3.hh"

// Class of Laue Group.
class ILaueGroup
{
public:
	virtual string CrystalSystemName() const = 0;
	virtual string Name() const = 0;
	virtual Int4 LaueMultiplicity(const MillerIndex3&) const = 0;
	virtual bool checkLatticeConstant(VecDat3<Double>&, VecDat3<Double>&, Double thred = 1.0e-3) const = 0;
	virtual ~ILaueGroup(){};
};


#endif /*ILAUEGROUP_H_*/
