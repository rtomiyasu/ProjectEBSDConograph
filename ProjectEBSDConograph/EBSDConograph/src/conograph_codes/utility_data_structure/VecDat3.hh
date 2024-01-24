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
#ifndef VECDAT3_H_
#define VECDAT3_H_

#include <complex>

using namespace std;

// Class of a vector of size 3 with some operators.
template <class T>
class VecDat3
{
protected:
	enum{ ISIZE = 3 };
	T m_vec[ISIZE];
	
public:
	VecDat3();
	VecDat3(const T&);
	VecDat3(const T&, const T&, const T&);
	VecDat3(const VecDat3<T>&); // copy constructor
	virtual ~VecDat3();
	VecDat3<T>& operator=(const T&);
	VecDat3<T>& operator=(const VecDat3<T>&);
	inline VecDat3<T>& operator+=(const VecDat3<T>&);
	inline VecDat3<T>& operator-=(const VecDat3<T>&);
	inline VecDat3<T>& operator*=(const T&);
	inline VecDat3<T> operator*(const T&) const;
	inline T& operator[](const int&);
	inline const T& operator[](const int&) const;
};

template <class T>
VecDat3<T>::VecDat3()
{
//	for(int k=0; k<ISIZE; k++) m_vec[k]=0;
}

template <class T>
VecDat3<T>::VecDat3(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=a;
}

template <class T>
VecDat3<T>::VecDat3(const T& a, const T& b, const T& c)
{
	m_vec[0]=a;
	m_vec[1]=b;
	m_vec[2]=c;
}

template <class T>
VecDat3<T>::VecDat3(const VecDat3<T> &rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=rhs[k];
}

template <class T>
VecDat3<T>::~VecDat3()
{
}

template <class T>
VecDat3<T>& VecDat3<T>::operator=(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]=a;
	return *this;
}

template <class T>
VecDat3<T>& VecDat3<T>::operator=(const VecDat3<T>& rhs)
{
	if(this != &rhs)
		for(int k=0; k<ISIZE; k++) m_vec[k]=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator+=(const VecDat3<T>& rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]+=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator-=(const VecDat3<T>& rhs)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]-=rhs[k];
	return *this;
}

template <class T>
inline VecDat3<T>& VecDat3<T>::operator*=(const T& a)
{
	for(int k=0; k<ISIZE; k++) m_vec[k]*=a;
	return *this;
}

template <class T>
inline VecDat3<T> VecDat3<T>::operator*(const T& a) const
{
	VecDat3<T> ans(*this);
	for(int k=0; k<ISIZE; k++) ans[k]*=a;
	return ans;
}

template <class T>
inline T& VecDat3<T>::operator[](const int& k)
{
	return m_vec[k];
}

template <class T>
inline const T& VecDat3<T>::operator[](const int& k) const
{
	return m_vec[k];
}

template <class T>
inline VecDat3<T> operator+(const VecDat3<T>& lhs, const VecDat3<T>& rhs)
{
	VecDat3<T> ans(lhs);
	for(int k=0; k<3; k++) ans[k]+=rhs[k];
	return ans;
}

template <class T>
inline VecDat3<T> operator-(const VecDat3<T>& rhs)
{
	VecDat3<T> ans;
	for(int k=0; k<3; k++) ans[k]=-rhs[k];
	return ans;
}

template <class T>
inline VecDat3<T> operator-(const VecDat3<T>& lhs, const VecDat3<T>& rhs)
{
	VecDat3<T> ans(lhs);
	for(int k=0; k<3; k++) ans[k]-=rhs[k];
	return ans;
}

template <class T>
inline VecDat3<T> real(const VecDat3< complex<T> >& rhs)
{
	VecDat3<T> ans;
	ans[0] = real(rhs[0]);
	ans[1] = real(rhs[1]);
	ans[2] = real(rhs[2]);
	return ans;
}

template <class T>
inline VecDat3< complex<T> > conj(const VecDat3< complex<T> >& rhs)
{
	VecDat3< complex<T> > ans;
	ans[0] = conj(rhs[0]);
	ans[1] = conj(rhs[1]);
	ans[2] = conj(rhs[2]);
	return ans;
}

#endif /*VECDAT3_H_*/
