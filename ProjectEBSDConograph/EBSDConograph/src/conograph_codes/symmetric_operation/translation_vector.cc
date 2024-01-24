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
#include"translation_vector.hh"

static const S1& S1_00()
{
	static const S1 S1_00(0,0);
	return S1_00;
}

static const S1& S1_12()
{
	static const S1 S1_12(1,2);
	return S1_12;
}

static const S1& S1_13()
{
	static const S1 S1_13(1,3);
	return S1_13;
}

static const S1& S1_23()
{
	static const S1 S1_23(2,3);
	return S1_23;
}

static const S1& S1_14()
{
	static const S1 S1_14(1,4);
	return S1_14;
}

static const S1& S1_34()
{
	static const S1 S1_34(3,4);
	return S1_34;
}

const VecDat3<S1>& exp2PIi000000()
{
	static const VecDat3<S1> exp2PIi000000(S1_00(),S1_00(),S1_00());
	return exp2PIi000000;
}

const VecDat3<S1>& exp2PIi121212()
{
	static const VecDat3<S1> exp2PIi121212(S1_12(),S1_12(),S1_12());
	return exp2PIi121212;
}

const VecDat3<S1>& exp2PIi141414()
{
	static const VecDat3<S1> exp2PIi141414(S1_14(),S1_14(),S1_14());
	return exp2PIi141414;
}

const VecDat3<S1>& exp2PIi343434()
{
	static const VecDat3<S1> exp2PIi343434(S1_34(),S1_34(),S1_34());
	return exp2PIi343434;
}

const VecDat3<S1>& exp2PIi120000()
{
	static const VecDat3<S1> exp2PIi120000(S1_12(),S1_00(),S1_00());
	return exp2PIi120000;
}

const VecDat3<S1>& exp2PIi001200()
{
	static const VecDat3<S1> exp2PIi001200(S1_00(),S1_12(),S1_00());
	return exp2PIi001200;
}

const VecDat3<S1>& exp2PIi000012()
{
	static const VecDat3<S1> exp2PIi000012(S1_00(),S1_00(),S1_12());
	return exp2PIi000012;
}

const VecDat3<S1>& exp2PIi121200()
{
	static const VecDat3<S1> exp2PIi121200(S1_12(),S1_12(),S1_00());
	return exp2PIi121200;
}

const VecDat3<S1>& exp2PIi120012()
{
	static const VecDat3<S1> exp2PIi120012(S1_12(),S1_00(),S1_12());
	return exp2PIi120012;
}

const VecDat3<S1>& exp2PIi001212()
{
	static const VecDat3<S1> exp2PIi001212(S1_00(),S1_12(),S1_12());
	return exp2PIi001212;
}

const VecDat3<S1>& exp2PIi140000()
{
	static const VecDat3<S1> exp2PIi140000(S1_14(),S1_00(),S1_00());
	return exp2PIi140000;
}

const VecDat3<S1>& exp2PIi001400()
{
	static const VecDat3<S1> exp2PIi001400(S1_00(),S1_14(),S1_00());
	return exp2PIi001400;
}

const VecDat3<S1>& exp2PIi000014()
{
	static const VecDat3<S1> exp2PIi000014(S1_00(),S1_00(),S1_14());
	return exp2PIi000014;
}

const VecDat3<S1>& exp2PIi141400()
{
	static const VecDat3<S1> exp2PIi141400(S1_14(),S1_14(),S1_00());
	return exp2PIi141400;
}

const VecDat3<S1>& exp2PIi140014()
{
	static const VecDat3<S1> exp2PIi140014(S1_14(),S1_00(),S1_14());
	return exp2PIi140014;
}

const VecDat3<S1>& exp2PIi001414()
{
	static const VecDat3<S1> exp2PIi001414(S1_00(),S1_14(),S1_14());
	return exp2PIi001414;
}

const VecDat3<S1>& exp2PIi340000()
{
	static const VecDat3<S1> exp2PIi340000(S1_34(),S1_00(),S1_00());
	return exp2PIi340000;
}

const VecDat3<S1>& exp2PIi003400()
{
	static const VecDat3<S1> exp2PIi003400(S1_00(),S1_34(),S1_00());
	return exp2PIi003400;
}

const VecDat3<S1>& exp2PIi000034()
{
	static const VecDat3<S1> exp2PIi000034(S1_00(),S1_00(),S1_34());
	return exp2PIi000034;
}

const VecDat3<S1>& exp2PIi343400()
{
	static const VecDat3<S1> exp2PIi343400(S1_34(),S1_34(),S1_00());
	return exp2PIi343400;
}

const VecDat3<S1>& exp2PIi340034()
{
	static const VecDat3<S1> exp2PIi340034(S1_34(),S1_00(),S1_34());
	return exp2PIi340034;
}

const VecDat3<S1>& exp2PIi003434()
{
	static const VecDat3<S1> exp2PIi003434(S1_00(),S1_34(),S1_34());
	return exp2PIi003434;
}

const VecDat3<S1>& exp2PIi341414()
{
	static const VecDat3<S1> exp2PIi341414(S1_34(),S1_14(),S1_14());
	return exp2PIi341414;
}

const VecDat3<S1>& exp2PIi143414()
{
	static const VecDat3<S1> exp2PIi143414(S1_14(),S1_34(),S1_14());
	return exp2PIi143414;
}

const VecDat3<S1>& exp2PIi141434()
{
	static const VecDat3<S1> exp2PIi141434(S1_14(),S1_14(),S1_34());
	return exp2PIi141434;
}

const VecDat3<S1>& exp2PIi343414()
{
	static const VecDat3<S1> exp2PIi343414(S1_34(),S1_34(),S1_14());
	return exp2PIi343414;
}

const VecDat3<S1>& exp2PIi341434()
{
	static const VecDat3<S1> exp2PIi341434(S1_34(),S1_14(),S1_34());
	return exp2PIi341434;
}

const VecDat3<S1>& exp2PIi143434()
{
	static const VecDat3<S1> exp2PIi143434(S1_14(),S1_34(),S1_34());
	return exp2PIi143434;
}

const VecDat3<S1>& exp2PIi003414()
{
	static const VecDat3<S1> exp2PIi003414(S1_00(),S1_34(),S1_14());
	return exp2PIi003414;
}

const VecDat3<S1>& exp2PIi140034()
{
	static const VecDat3<S1> exp2PIi140034(S1_14(),S1_00(),S1_34());
	return exp2PIi140034;
}

const VecDat3<S1>& exp2PIi341400()
{
	static const VecDat3<S1> exp2PIi341400(S1_34(),S1_14(),S1_00());
	return exp2PIi341400;
}

const VecDat3<S1>& exp2PIi001434()
{
	static const VecDat3<S1> exp2PIi001434(S1_00(),S1_14(),S1_34());
	return exp2PIi001434;
}

const VecDat3<S1>& exp2PIi340014()
{
	static const VecDat3<S1> exp2PIi340014(S1_34(),S1_00(),S1_14());
	return exp2PIi340014;
}

const VecDat3<S1>& exp2PIi143400()
{
	static const VecDat3<S1> exp2PIi143400(S1_14(),S1_34(),S1_00());
	return exp2PIi143400;
}

const VecDat3<S1>& exp2PIi123414()
{
	static const VecDat3<S1> exp2PIi123414(S1_12(),S1_34(),S1_14());
	return exp2PIi123414;
}

const VecDat3<S1>& exp2PIi141234()
{
	static const VecDat3<S1> exp2PIi141234(S1_14(),S1_12(),S1_34());
	return exp2PIi141234;
}

const VecDat3<S1>& exp2PIi341412()
{
	static const VecDat3<S1> exp2PIi341412(S1_34(),S1_14(),S1_12());
	return exp2PIi341412;
}

const VecDat3<S1>& exp2PIi121434()
{
	static const VecDat3<S1> exp2PIi121434(S1_12(),S1_14(),S1_34());
	return exp2PIi121434;
}

const VecDat3<S1>& exp2PIi341214()
{
	static const VecDat3<S1> exp2PIi341214(S1_34(),S1_12(),S1_14());
	return exp2PIi341214;
}

const VecDat3<S1>& exp2PIi143412()
{
	static const VecDat3<S1> exp2PIi143412(S1_14(),S1_34(),S1_12());
	return exp2PIi143412;
}

const VecDat3<S1>& exp2PIi001234()
{
	static const VecDat3<S1> exp2PIi001234(S1_00(),S1_12(),S1_34());
	return exp2PIi001234;
}

const VecDat3<S1>& exp2PIi120034()
{
	static const VecDat3<S1> exp2PIi120034(S1_12(),S1_00(),S1_34());
	return exp2PIi120034;
}

const VecDat3<S1>& exp2PIi121234()
{
	static const VecDat3<S1> exp2PIi121234(S1_12(),S1_12(),S1_34());
	return exp2PIi121234;
}

const VecDat3<S1>& exp2PIi001214()
{
	static const VecDat3<S1> exp2PIi001214(S1_00(),S1_12(),S1_14());
	return exp2PIi001214;
}

const VecDat3<S1>& exp2PIi120014()
{
	static const VecDat3<S1> exp2PIi120014(S1_12(),S1_00(),S1_14());
	return exp2PIi120014;
}

const VecDat3<S1>& exp2PIi121214()
{
	static const VecDat3<S1> exp2PIi121214(S1_12(),S1_12(),S1_14());
	return exp2PIi121214;
}

const VecDat3<S1>& exp2PIi000013()
{
	static const VecDat3<S1> exp2PIi000013(S1_00(),S1_00(),S1_13());
	return exp2PIi000013;
}

const VecDat3<S1>& exp2PIi000023()
{
	static const VecDat3<S1> exp2PIi000023(S1_00(),S1_00(),S1_23());
	return exp2PIi000023;
}

const VecDat3<S1>& exp2PIi231313()
{
	static const VecDat3<S1> exp2PIi231313(S1_23(),S1_13(),S1_13());
	return exp2PIi231313;
}

const VecDat3<S1>& exp2PIi132323()
{
	static const VecDat3<S1> exp2PIi132323(S1_13(),S1_23(),S1_23());
	return exp2PIi132323;
}
