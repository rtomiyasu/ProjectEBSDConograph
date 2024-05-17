/*
 * The MIT License

   EBSD-Conograph (Conograph software for EBSD ab-initio indexing)

Copyright (c) <2020> <Ryoko Oishi-Tomiyasu, Kyushu University>

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

#ifndef UTILITY_RATIONALNUMBER_HH_
#define UTILITY_RATIONALNUMBER_HH_

#include <string>
#include "../HC_algorithm/GCD_1_3_1.hh"

using namespace std;

template<class Integer>
class RationalNumber
{
	template<class Integer2>
	friend inline Integer2 floor(const RationalNumber<Integer2>& rhs);

	template<class Integer2>
	friend inline Integer2 lrint(const RationalNumber<Integer2>& rhs);

private:
	// A number represented by the fraction m_numer / m_denom.
	Integer m_numer;
	Integer m_denom;

	static Integer floor_int(const Integer& a, const Integer& b)
	{
		assert (b != 0);
		const Integer c = abs(a) / abs(b);
		if ( (a >= 0 && b >0) || (a < 0 && b < 0)  ) return c;
		else if (a % b == 0) return -c;
		else return -c - 1;
	}

	static Integer lrint_int(const Integer& a, const Integer& b)
	{
		return floor_int(a*2+b, b*2);
	}

public:
	RationalNumber(){};
	RationalNumber(const Integer& arg1);
	RationalNumber(const Integer& arg1, const Integer& arg2);
	~RationalNumber(){};

	inline RationalNumber<Integer>& operator+=(const RationalNumber<Integer>& rhs);
	inline RationalNumber<Integer>& operator-=(const RationalNumber<Integer>& rhs);
	inline RationalNumber<Integer>& operator*=(const Integer& rhs);
	inline RationalNumber<Integer>& operator/=(const Integer& rhs);
	inline RationalNumber<Integer>& operator*=(const RationalNumber<Integer>&);
	inline RationalNumber<Integer>& operator/=(const RationalNumber<Integer>&);

	inline const Integer& Numerator() const { return m_numer; };
	inline const Integer& Denominator() const { return m_denom; };
	inline const RationalNumber<Integer>& Reduce() { const Integer gcd = GCD_Euclid(m_numer, m_denom)*(m_denom < 0?-1:1); m_numer/=gcd; m_denom/=gcd; return *this; };
	inline RationalNumber<Integer> moduloZ() const;
	double toDouble() const { return double(m_numer) / double(m_denom); };
	string toString() const;

	inline bool operator==(const Integer& rhs) const { return m_numer == m_denom*rhs; };
	inline bool operator!=(const Integer& rhs) const { return m_numer != m_denom*rhs; };
	inline bool operator<(const Integer& rhs) const { if( m_denom < 0 ) return m_numer > m_denom*rhs; else return m_numer < m_denom*rhs; };
	inline bool operator<=(const Integer& rhs) const { if( m_denom < 0 ) return m_numer >= m_denom*rhs; else return m_numer <= m_denom*rhs; };

	inline bool operator==(const RationalNumber<Integer>& rhs) const { return *this - rhs == 0; };
	inline bool operator!=(const RationalNumber<Integer>& rhs) const { return *this - rhs != 0; };
	inline bool operator<=(const RationalNumber<Integer>& rhs) const { return *this - rhs <= 0; };
	inline bool operator<(const RationalNumber<Integer>& rhs) const { return *this - rhs < 0; };

	inline Integer floor() const { return floor_int(m_numer, m_denom); };
	inline Integer lrint() const { return lrint_int(m_numer, m_denom); };
};


template<class Integer>
RationalNumber<Integer>::RationalNumber(const Integer& arg1) : m_numer(arg1), m_denom(1)
{
}

template<class Integer>
RationalNumber<Integer>::RationalNumber(const Integer& arg1, const Integer& arg2) : m_numer(arg1), m_denom(arg2)
{
}

template<class Integer>
string RationalNumber<Integer>::toString() const
{
	const Integer igcd = GCD_Euclid_For_NonNegative_Integers(m_numer, m_denom);
	const Integer a = m_numer / igcd;
	const Integer b = m_denom / igcd;

	if(a % b != 0)
	{
		 if( b < 0 )
		 {
			 return num2str(-a) + "/" + num2str(abs(b));
		 }
		 else
		 {
			 return num2str(a) + "/" + num2str(b);
		 }
	}
	return num2str(a/b);
}


template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator+=(const RationalNumber<Integer>& rhs)
{
	// Compute this->m_numer / this->m_denom + rhs.m_numer / rhs.denom.
	const Integer igcd = GCD_Euclid_For_NonNegative_Integers(this->m_denom, rhs.m_denom); // Calculate the greatest common devisor.
	const Integer c = this->m_denom / igcd;
	this->m_numer = this->m_numer * (rhs.m_denom / igcd) + rhs.m_numer * c;
	this->m_denom = c * rhs.m_denom;
	return *this;
}

template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator-=(const RationalNumber<Integer>& rhs)
{
	// Compute this->m_numer / this->m_denom - rhs.m_numer / rhs.denom.
	const Integer igcd = GCD_Euclid_For_NonNegative_Integers(this->m_denom, rhs.m_denom); // Calculate the greatest common devisor.
	const Integer c = this->m_denom / igcd;
	this->m_numer = this->m_numer * (rhs.m_denom / igcd) - rhs.m_numer * c;
	this->m_denom = c * rhs.m_denom;
	return *this;
}

template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator*=(const Integer& rhs)
{
	const Integer igcd = GCD_Euclid_For_NonNegative_Integers(this->m_denom, rhs); // Calculate the greatest common devisor.
	this->m_numer=this->m_numer * (rhs / igcd);
	this->m_denom/=igcd;
	return *this;
}

template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator/=(const Integer& rhs)
{
	const Integer igcd = GCD_Euclid_For_NonNegative_Integers(this->m_numer, rhs); // Calculate the greatest common devisor.
	this->m_numer/=igcd;
	this->m_denom=this->m_denom * (rhs / igcd);
	return *this;
}

template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator*=(const RationalNumber<Integer>& rhs)
{
	*this *= rhs.m_numer;
	*this /= rhs.m_denom;
	return *this;
}

template<class Integer>
inline RationalNumber<Integer>& RationalNumber<Integer>::operator/=(const RationalNumber<Integer>& rhs)
{
	*this *= rhs.m_denom;
	*this /= rhs.m_numer;
	return *this;
}

template<class Integer>
inline RationalNumber<Integer> RationalNumber<Integer>::moduloZ() const
{
	Integer numer_new = m_numer*(m_denom<0?-1:1);
	const Integer denom_new = abs(m_denom);
	numer_new = numer_new % denom_new;
	if( numer_new >= 0 ) return RationalNumber(numer_new, denom_new);
	else return RationalNumber(denom_new + numer_new, denom_new);
}

template<class Integer>
inline RationalNumber<Integer> operator+(const RationalNumber<Integer>& lhs, const RationalNumber<Integer>& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans += rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator-(const RationalNumber<Integer>& lhs, const RationalNumber<Integer>& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans -= rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator*(const RationalNumber<Integer>& lhs, const Integer& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans *= rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator/(const RationalNumber<Integer>& lhs, const Integer& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans /= rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator*(const RationalNumber<Integer>& lhs, const RationalNumber<Integer>& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans *= rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator/(const RationalNumber<Integer>& lhs, const RationalNumber<Integer>& rhs)
{
	RationalNumber<Integer> ans(lhs);
	ans /= rhs;
	return ans;
}

template<class Integer>
inline RationalNumber<Integer> operator-(const RationalNumber<Integer>& rhs)
{
	return RationalNumber<Integer>(-rhs.Numerator(), rhs.Denominator());
}

template<class Integer>
inline RationalNumber<Integer> abs(const RationalNumber<Integer>& rhs)
{
	return RationalNumber<Integer>(abs(rhs.Numerator()), abs(rhs.Denominator()));
}

template<class Integer>
inline string num2str(const RationalNumber<Integer>& rhs)
{
	return rhs.chToString();
}

template<class Integer>
inline Integer floor(const RationalNumber<Integer>& rhs)
{
	return RationalNumber<Integer>::floor_int(rhs.Numerator(), rhs.Denominator());
}

template<class Integer>
inline Integer lrint(const RationalNumber<Integer>& rhs)
{
	return RationalNumber<Integer>::lrint_int(rhs.Numerator(), rhs.Denominator());
}

#endif /* UTILITY_RATIONALNUMBER_HH_ */
