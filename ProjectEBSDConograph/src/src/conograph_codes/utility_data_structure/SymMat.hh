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
#ifndef SYMMAT_H_
#define SYMMAT_H_

#include <assert.h>
#include "../utility_data_structure/nrutil_nr.hh"

using namespace std;

// Class of a symmmetric matrix (a_ij) determined by m_mat as below.
//          m_mat[0] m_mat[1] m_mat[2]
// (a_ij) = m_mat[1] m_mat[3] m_mat[4]
//          m_mat[2] m_mat[4] m_mat[5]
template <class T>
class SymMat
{
private:
	int ISIZE;
	int NUM_ELEMENT;
	T* m_mat;
	
	inline T& operator[](const int&);
	inline const T& operator[](const int&) const;

public:
	SymMat(){ ISIZE=0; NUM_ELEMENT=0; m_mat=0; };
	SymMat(const int& isize);
	SymMat(const int& isize, const T&);
	SymMat(const SymMat<T>&); // copy constructor
	~SymMat();

	SymMat<T>& operator=(const T&);
	SymMat<T>& operator=(const SymMat<T>&);
	SymMat<T>& operator+=(const SymMat<T>&);
	SymMat<T>& operator-=(const SymMat<T>&);
	SymMat<T>& operator*=(const T&);
	template<class U>
	SymMat<T>& operator/=(const U&);
	
	inline T& operator()(const int&, const int&);
	inline const T& operator()(const int&, const int&) const;

	inline const int& size() const { return ISIZE; };

	// for GUI
	const T *getref_m_mat() const {return m_mat;}
	      T *getref_m_mat()       {return m_mat;}
};

template <class T>
SymMat<T>::SymMat(const int& isize) : ISIZE(isize), NUM_ELEMENT(isize*(isize+1)/2), m_mat(NULL)
{
	assert( isize >= 0 );
	if( isize > 0 ) m_mat = new T[NUM_ELEMENT];
}

template <class T>
SymMat<T>::SymMat(const int& isize, const T& a) : ISIZE(isize), NUM_ELEMENT(isize*(isize+1)/2), m_mat(NULL)
{
	assert( isize >= 0 );
	if( isize > 0 ) m_mat = new T[NUM_ELEMENT];
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]=a;
}

template <class T>
SymMat<T>::SymMat(const SymMat<T>&rhs) : ISIZE(rhs.ISIZE), NUM_ELEMENT(ISIZE*(ISIZE+1)/2), m_mat(NULL)
{
	if( ISIZE > 0 ) m_mat = new T[NUM_ELEMENT];
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]=rhs.m_mat[k];
}

template <class T>
SymMat<T>::~SymMat()
{
	delete[] (m_mat);
}

template <class T>
SymMat<T>& SymMat<T>::operator=(const T& a)
{
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]=a;
	return *this;
}

template <class T>
SymMat<T>& SymMat<T>::operator=(const SymMat<T> &rhs)
{
	if (this != &rhs)
	{
		if (ISIZE != rhs.ISIZE) {
			delete [] m_mat;
			m_mat = 0;
			ISIZE=rhs.ISIZE;
			NUM_ELEMENT=rhs.NUM_ELEMENT;
			if( NUM_ELEMENT > 0 ) m_mat = new T[NUM_ELEMENT];
		}
		for (int i=0; i<NUM_ELEMENT; i++)
			m_mat[i]=rhs[i];
	}
	return *this;
}

template <class T>
SymMat<T>& SymMat<T>::operator+=(const SymMat<T>& rhs)
{
	assert( ISIZE == rhs.ISIZE );
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]+=rhs.m_mat[k];
	return *this;
}

template <class T>
SymMat<T>& SymMat<T>::operator-=(const SymMat<T>& rhs)
{
	assert( ISIZE == rhs.ISIZE );
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]-=rhs.m_mat[k];
	return *this;
}

template <class T>
SymMat<T>& SymMat<T>::operator*=(const T& a)
{
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]*=a;
	return *this;
}

template <class T> template<class U>
SymMat<T>& SymMat<T>::operator/=(const U& a)
{
	for(int k=0; k<NUM_ELEMENT; k++) m_mat[k]/=a;
	return *this;
}

template <class T>
inline SymMat<T> operator+(const SymMat<T>& lhs, const SymMat<T>& rhs)
{
	SymMat<T> ans(lhs);
	ans += rhs;
	return ans;
}

template <class T>
inline SymMat<T> operator-(const SymMat<T>& lhs, const SymMat<T>& rhs)
{
	SymMat<T> ans(lhs);
	ans -= rhs;
	return ans;
}

template <class T>
inline SymMat<T> operator*(const SymMat<T>& lhs, const T& a)
{
	SymMat<T> ans(lhs);
	ans *= a;
	return ans;
}

template <class T, class U>
inline SymMat<T> operator/(const SymMat<T>& lhs, const U& a)
{
	SymMat<T> ans(lhs);
	ans /= a;
	return ans;
}

template <class T>
inline T& SymMat<T>::operator[](const int& j)
{
	assert( 0 <= j && j < NUM_ELEMENT );
	return m_mat[j];
}

template <class T>
inline const T& SymMat<T>::operator[](const int& j) const
{
	assert( 0 <= j && j < NUM_ELEMENT );
	return m_mat[j];
}


inline int put_index(const int& ISIZE, const int& j, const int& k)
{
	if( j < k ) return j*(ISIZE*2-j-1)/2 + k;
	else return k*(ISIZE*2-k-1)/2 + j;
}


template <class T>
inline T& SymMat<T>::operator()(const int& j, const int& k)
{
//	assert( 0 <= j && j < ISIZE );
//	assert( 0 <= k && k < ISIZE );
	return m_mat[put_index(ISIZE, j, k)];
}

template <class T>
inline const T& SymMat<T>::operator()(const int& j, const int& k) const
{
//	assert( 0 <= j && j < ISIZE );
//	assert( 0 <= k && k < ISIZE );
	return m_mat[put_index(ISIZE, j, k)];
}


inline double Determinant3(const SymMat<double>& rhs)
{
	assert( rhs.size() == 3 );

	return rhs(0,0)*rhs(1,1)*rhs(2,2)
			- rhs(1,2)*rhs(1,2)*rhs(0,0) - rhs(0,2)*rhs(0,2)*rhs(1,1) - rhs(0,1)*rhs(0,1)*rhs(2,2)
			+ rhs(0,1)*rhs(0,2)*rhs(1,2)*2.0;
}


inline SymMat<double> Inverse3(const SymMat<double>& rhs)
{
	assert( rhs.size() == 3 );

	const double det01 = rhs(0,0)*rhs(1,1)-rhs(0,1)*rhs(1,0);
	const double det02 = rhs(0,0)*rhs(2,2)-rhs(0,2)*rhs(2,0);
	const double det12 = rhs(1,1)*rhs(2,2)-rhs(1,2)*rhs(2,1);
	const double det01_02 = rhs(0,0)*rhs(1,2)-rhs(0,2)*rhs(1,0);
	const double det01_12 = rhs(0,1)*rhs(1,2)-rhs(0,2)*rhs(1,1);
	const double det02_12 = rhs(0,1)*rhs(2,2)-rhs(0,2)*rhs(2,1);
	
	const double det = 1.0/(rhs(0,0)*det12 - rhs(0,1)*det02_12 + rhs(0,2)*det01_12);
	
	SymMat<double> ans(3);
	ans(0,0) = det12;
	ans(1,1) = det02;
	ans(2,2) = det01;
	ans(0,1) = -det02_12;
	ans(0,2) = det01_12;
	ans(1,2) = -det01_02;
	
	ans *= det;
	return ans;
}

template <class U, class T>
inline SymMat<T> transform_sym_matrix(const NRMat<U>& lhs, const SymMat<T>& rhs)
{
	const int icol = lhs.ncols();
	const int isize = lhs.nrows();

	assert( rhs.size() == icol );

	SymMat<T> ans(isize, 0);
	for(int i=0; i<isize; i++)
	{
		for(int i2=i; i2<isize; i2++)
		{
			for(int j=0; j<icol; j++)
			{
				T& ans_r = ans(i,i2);
				ans_r += rhs(j,j)*( lhs[i][j]*lhs[i2][j] );
				for(int j2=j+1; j2<icol; j2++)
				{
					ans_r += rhs(j,j2)*( lhs[i][j]*lhs[i2][j2] + lhs[i][j2]*lhs[i2][j] );
				}
			}
		}
	}
	return ans;
}


#endif /*SymMat_H_*/
