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
#ifndef HKL_Q_HH_
#define HKL_Q_HH_

#include "../ZAnalysisTypes.hh"
#include "../utility_data_structure/nrutil_nr.hh"
#include "../utility_data_structure/VecDat3.hh"

class HKL_Q
{
private:
	VecDat3<Int4> hkl;
	Double q;

public:
	HKL_Q(){ q = 0.0; };
	HKL_Q(const VecDat3<Int4>& arg1, const Double& arg2): hkl(arg1), q(arg2){};
	HKL_Q(const NRVec<Int4>& arg1, const Double& arg2): hkl(arg1[0], arg1[1], arg1[2]), q(arg2){ assert(arg1.size()==3); };

	inline const VecDat3<Int4>& HKL() const { return hkl; };
	inline void setHKL(const VecDat3<Int4>& arg) { hkl[0] = arg[0]; hkl[1] = arg[1]; hkl[2] = arg[2]; };
	inline const Double& Q() const { return q; };

	static bool cmp(const HKL_Q& lhs, const HKL_Q& rhs) { return lhs.q < rhs.q; };

	// For GUI
	const VecDat3<Int4> &getref_hkl() const {return hkl;}
		  VecDat3<Int4> &getref_hkl()       {return hkl;}
	const double        &getref_q()   const {return q;} 
	      double        &getref_q()         {return q;} 

};

inline bool operator<(const HKL_Q& lhs, const HKL_Q& rhs)
{
	return lhs.Q() < rhs.Q();
};

#endif /* HKL_Q_HH_ */
