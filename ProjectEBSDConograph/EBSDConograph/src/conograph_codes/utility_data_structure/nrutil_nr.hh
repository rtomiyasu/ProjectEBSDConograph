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
#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace std;

template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(int n, const T &a);	//initialize to constant value
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();

	NRVec& operator+=(const NRVec& rhs);
	NRVec& operator-=(const NRVec& rhs);
//	NRVec& operator*=(const NRVec&);
	NRVec& operator*=(const T&);
	NRVec& operator/=(const T&);

//	inline T& at(const int i);	//subscripting: pointer to row i
//	inline const T& at(const int i) const;
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(0)
{
	if( n > 0 )	v = new T[n];
}

template <class T>
NRVec<T>::NRVec(int n, const T& a) : nn(n), v(0)
{
	if( n > 0 )	v = new T[n];
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(0)
{
	if( rhs.nn > 0 ) v = new T[nn];
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			delete [] (v);
			v = 0;
			nn=rhs.nn;
			if( nn > 0 ) v = new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
//	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
	assert( 0 <= i && i < nn );
	return v[i];
}

//template <class T>
//inline T& NRVec<T>::at(const int i)	//subscripting: pointer to row i
//{
//	if( i < 0 || i >= nn ) throw nerror(ZErrorArrayOverFlow, __FILE__, __LINE__, __FUNCTION__);
//	return v[i];
//}
//
//template <class T>
//inline const T& NRVec<T>::at(const int i) const
//{
//	if( i < 0 || i >= nn ) throw nerror(ZErrorArrayOverFlow, __FILE__, __LINE__, __FUNCTION__);
//	return v[i];
//}

template <class T>
inline int NRVec<T>::size() const
{
	return nn;
}

template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		delete[] (v);
}

template <class T>
NRVec<T>& NRVec<T>::operator+=(const NRVec<T>& rhs)
{
	assert( nn ==rhs.size() );
	for(int k=0; k<nn; k++) v[k]+=rhs[k];
	return *this;
}

template <class T>
NRVec<T> operator+(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	NRVec<T> ans(lhs);
	ans+=rhs;
	return ans;
}

template <class T>
NRVec<T>& NRVec<T>::operator-=(const NRVec<T>& rhs)
{
	assert( nn ==rhs.size() );
	for(int k=0; k<nn; k++) v[k]-=rhs[k];
	return *this;
}

template <class T>
NRVec<T> operator-(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	NRVec<T> ans(lhs);
	ans-=rhs;
	return ans;
}

//template <class T>
//NRVec<T>& NRVec<T>::operator*=(const NRVec<T>& rhs)
//{
//	assert( nn ==rhs.size() );
//	for(int k=0; k<nn; k++) v[k]*=rhs[k];
//	return *this;
//}


template <class T>
NRVec<T>& NRVec<T>::operator*=(const T& rhs)
{
	for(int k=0; k<nn; k++) v[k]*=rhs;
	return *this;
}

template <class T>
NRVec<T> operator*(const NRVec<T>& lhs, const T& rhs)
{
	NRVec<T> ans(lhs);
	ans*=rhs;
	return ans;
}

template <class T>
NRVec<T>& NRVec<T>::operator/=(const T& rhs)
{
	for(int k=0; k<nn; k++) v[k]/=rhs;
	return *this;
}

template <class T>
T inner_product(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	const int isize = lhs.size();
	assert( isize==rhs.size() );
	
	T ans = 0;
	for(int k=0; k<isize; k++) ans += lhs[k]*rhs[k];
	return ans;
}

template <class T>
class NRMat {
private:
	int nn;
	int mm;
	T **v;
public:
	NRMat();
	NRMat(int n, int m);			// Zero-based array
	NRMat(int n, int m, const T &a);	//Initialize to constant
	NRMat(const NRMat &rhs);		// Copy constructor
	NRMat & operator=(const NRMat &rhs);	//assignment
	NRMat & operator=(const T &a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	~NRMat();

	NRMat& operator+=(const NRMat&);
	NRMat& operator-=(const NRMat&);
	NRMat& operator*=(const T&);
};

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}

template <class T>
NRMat<T>::NRMat(int n, int m) : nn(n), mm(m), v(0)
{
	if(n > 0)
	{
		v = new T*[n];
		if(m > 0)
		{
			v[0] = new T[m*n];
			for (int i=1; i< n; i++)
				v[i] = v[i-1] + m;
		}
		else for (int i=0; i< n; i++) v[i] = 0;
	}
}

template <class T>
NRMat<T>::NRMat(int n, int m, const T &a) : nn(n), mm(m), v(0)
{
	int i,j;
	if(n > 0)
	{
		v = new T*[n];
		if(m > 0)
		{
			v[0] = new T[m*n];
			for (i=1; i< n; i++)
				v[i] = v[i-1] + m;
			for (i=0; i< n; i++)
				for (j=0; j<m; j++)
					v[i][j] = a;
		}
		else for (int i=0; i< n; i++) v[i] = 0;
	}
}

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(0)
{
	int i,j;
	if(nn > 0)
	{
		v = new T*[nn];
		if(mm > 0)
		{
			v[0] = new T[mm*nn];
			for (i=1; i< nn; i++)
				v[i] = v[i-1] + mm;
			for (i=0; i< nn; i++)
				for (j=0; j<mm; j++)
					v[i][j] = rhs[i][j];
		}
		else for (int i=0; i<nn; i++) v[i] = 0;
	}
}


template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if(nn > 0) delete[] (v[0]);
			delete[] (v);
			v = 0;

			nn=rhs.nn;
			mm=rhs.mm;
			if( nn > 0 )
			{
				v = new T*[nn];
				if( mm > 0 )
				{
					v[0] = new T[mm*nn];
					for (i=1; i< nn; i++)
						v[i] = v[i-1] + mm;
				}
				else for (i=0; i< nn; i++) v[i] = 0;
			}
		}
		for (i=0; i< nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}



template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* NRMat<T>::operator[](const int i)	//subscripting: pointer to row i
{
//	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
	assert( 0 <= i && i < nn );
	return v[i];
}

template <class T>
inline int NRMat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
	return mm;
}

template <class T>
NRMat<T>::~NRMat()
{
	if(nn > 0) delete[] (v[0]);
	delete[] (v);
}

template <class T>
NRMat<T>& NRMat<T>::operator+=(const NRMat<T>& rhs)
{
	assert( nn==rhs.nrows() );
	assert( mm!=rhs.ncols() );

	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] += rhs[k][j];
	return *this;
}

template <class T>
NRMat<T>& NRMat<T>::operator-=(const NRMat<T>& rhs)
{
	assert( nn==rhs.nrows() );
	assert( mm!=rhs.ncols() );

	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] -= rhs[k][j];
	return *this;
}

template <class T>
NRMat<T>& NRMat<T>::operator*=(const T& rhs)
{
	for(int k=0; k<nn; k++)
		for(int j=0; j<mm; j++) v[k][j] *= rhs;
	return *this;
}

template <class T>
NRMat<T> operator*(const NRMat<T>& lhs, const T& rhs)
{
	NRMat<T> ans = lhs;
	for(int k=0; k<ans.nrows(); k++)
		for(int j=0; j<ans.ncols(); j++) ans[k][j] *= rhs;
	return ans;
}

template <class T>
NRMat<T> mprod(const NRMat<T>& lhs, const NRMat<T>& rhs)
{
	const int irow = lhs.nrows();
	const int icol = rhs.ncols();
	const int n = lhs.ncols();
	assert( lhs.ncols()==rhs.nrows() );

	
	NRMat<T> ans( irow, icol, 0 );
	for(int k=0; k<irow; k++)
	{
		for(int j=0; j<icol; j++)
		{
			for(int i=0; i<n; i++)
			{
				ans[k][j] += lhs[k][i]*rhs[i][j];
			}
		}
	}

	return ans;
}

template<class T>
inline void transpose_square_matrix(NRMat<T>& rhs)
{
	const int isize = rhs.nrows();
	assert( isize == rhs.ncols() );

	for(int i=0; i<isize; i++)
		for(int j=0; j<i; j++) swap(rhs[i][j], rhs[j][i]);
}

template<class T>
inline T Determinant3(const NRMat<T>& rhs)
{
	assert( rhs.nrows() == 3 && rhs.ncols() == 3 );

	const T det12 = rhs[1][1]*rhs[2][2]-rhs[1][2]*rhs[2][1];
	const T det12_02 = rhs[1][0]*rhs[2][2]-rhs[1][2]*rhs[2][0];
	const T det12_01 = rhs[1][0]*rhs[2][1]-rhs[1][1]*rhs[2][0];

	return rhs[0][0]*det12 - rhs[0][1]*det12_02 + rhs[0][2]*det12_01;
}

inline NRMat<int> Inverse3(const NRMat<int>& rhs)
{
	assert( rhs.nrows() == 3 && rhs.ncols() == 3 );

	const int det01 = rhs[0][0]*rhs[1][1]-rhs[0][1]*rhs[1][0];
	const int det02 = rhs[0][0]*rhs[2][2]-rhs[0][2]*rhs[2][0];
	const int det12 = rhs[1][1]*rhs[2][2]-rhs[1][2]*rhs[2][1];
	const int det01_02 = rhs[0][0]*rhs[1][2]-rhs[0][2]*rhs[1][0];
	const int det02_01 = rhs[0][0]*rhs[2][1]-rhs[0][1]*rhs[2][0];
	const int det01_12 = rhs[0][1]*rhs[1][2]-rhs[0][2]*rhs[1][1];
	const int det02_12 = rhs[0][1]*rhs[2][2]-rhs[0][2]*rhs[2][1];
	const int det12_02 = rhs[1][0]*rhs[2][2]-rhs[1][2]*rhs[2][0];
	const int det12_01 = rhs[1][0]*rhs[2][1]-rhs[1][1]*rhs[2][0];

	const int det = rhs[0][0]*det12 - rhs[0][1]*det12_02 + rhs[0][2]*det12_01;
	assert( abs(det) == 1 );

	NRMat<int> ans(3,3);
	ans[0][0] = det12;
	ans[0][1] = -det02_12;
	ans[0][2] = det01_12;
	ans[1][0] = -det12_02;
	ans[1][1] = det02;
	ans[1][2] = -det01_02;
	ans[2][0] = det12_01;
	ans[2][1] = -det02_01;
	ans[2][2] = det01;

	if( det == 1 )
	{
		return ans;
	}
	else
	{
		return ans * (-1);
	}
}

template <class T>
inline NRMat<T> transpose(const NRMat<T>& rhs)
{
	NRMat<T> ans(rhs.ncols(), rhs.nrows());
	for(int i=0; i<rhs.nrows(); i++)
	{
		for(int j=0; j<rhs.ncols(); j++)
		{
			ans[j][i] = rhs[i][j];
		}
	}
	return ans;
}


template <class T>
NRMat<T> transpose_product(const NRVec<T>& lhs, const NRVec<T>& rhs)
{
	const int irow = lhs.size();
	const int icol = lhs.size();
	
	NRMat<T> ans( irow, icol );
	for(int k=0; k<irow; k++)
		for(int j=0; j<icol; j++) ans[k][j] = lhs[k]*rhs[j];

	return ans;
}

inline bool operator==(const NRMat<int> &m1, const NRMat<int> &m2) {
	if(m1.ncols() != m2.ncols()) { return false; }
	if(m1.nrows() != m2.nrows()) { return false; }
	for(int i=0; i<m1.nrows(); i++) {
		for(int j=0; j<m1.ncols(); j++) {
			if(m1[i][j] != m2[i][j]) { return false; }
		}
	}
	return true;
}

template <class T>
inline NRMat<T> identity_matrix(const int& isize)
{
	NRMat<T> TransMat(isize,isize,0);
	for(int k=0; k<isize; k++) TransMat[k][k] = 1;

	return TransMat;
}

#endif /* _NR_UTIL_H_ */
