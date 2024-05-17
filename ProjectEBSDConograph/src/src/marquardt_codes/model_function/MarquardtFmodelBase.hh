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
#ifndef MARQUARDTFMODELBASE_HH_
#define MARQUARDTFMODELBASE_HH_

#include <algorithm>
#ifndef __GNUC__
	#include <float.h>
#endif
#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../levenberg_marquardt/SVdcmp.hh"
#include "../zerror_type/error_out.hh"
#include "IMarquardtFmodel.hh"

template<class T>
class MarquardtFmodelBase
{
protected:
//	bool m_param_view_flag;
//    bool m_weight_flag;          // 0:weight proportional(not equal) to standard error.
	
	IMarquardtFmodel<T>* m_func;
	Mat_DP_constr m_Cmat;
	
	vector< vector<T> > m_xdata;
	vector< Vec_DP > m_ydata;
	vector< Vec_DP > m_inv_edata;
	
	Int4 m_degf_all;
	Int4 m_num_indep_param_all;
	Int4 m_num_data_all;
	Vec_INT m_num_indep_param;
	Vec_INT m_degf;
//	Vec_INT m_num_data;

	bool calJacobian(NRMat<Double>& Jacob, NRVec<Double>& R,
								Double& chisq_all, Vec_DP& chisq) const;
//	bool calJacobianWithWeight(const IMarquardtFmodel& func, NRMat<Double>& Jacob, NRVec<Double>& R);

	static void calLMMatrixJ(const NRMat<Double>& Jacob, const NRVec<Double>& R,
								SymMat<Double>& alpha, NRVec<Double>& beta);

public:
	MarquardtFmodelBase(IMarquardtFmodel<T>& func);
	virtual ~MarquardtFmodelBase(){};

	inline Int4 putHistogramNum() const { return m_xdata.size(); };
	virtual bool putSeparableLSFlag(){ return false; };
	
	inline Int4 putParamNum() const { return m_func->putParamNum(); };
	inline Int4 putNumberOfIndependentParamAll() const { return m_num_indep_param_all; };
	inline Int4 putNumberOfIndependentParam(const Int4& index) const { assert( 0 <= index&& index < putHistogramNum() ); return m_num_indep_param[index]; };
	inline Int4 putDegreeOfFreedomAll() const { return m_degf_all; };
	inline Int4 putDegreeOfFreedom(const Int4&) const;
	
	virtual Int4 putNumberOfNonlinearIndependentParamAll() const { return m_num_indep_param_all; };

	inline const Mat_DP_constr& putConstraint() const { return m_Cmat; };
	
	virtual Int4 setParam(NRVec<Double>& paratry,
							ostringstream& outPutString){ return m_func->setParamErrOut(&paratry[0], outPutString); };
	virtual Int4 setParamAll(NRVec<Double>& paratry0,
							ostringstream& outPutString){ return setParam(paratry0, outPutString); };
	virtual Int4 setParam(const NRVec<Double>& dpara, NRVec<Double>& paratry,
							ostringstream& outPutString);
	
	virtual bool optimizeLinearParam(const NRVec<Double>& paratry0, const Double& eps1,
									SymMat<Double>& alpha,
									NRVec<Double>& beta, 
									NRMat<Double>& Lmat, 
									NRVec<Double>& adiag,
									Double& chisq_all,
									Vec_DP& chisq){ throw ZErrorMessage(ZErrorEscapeFunction, __FILE__, __LINE__, __FUNCTION__); return 0; };

	virtual bool calLMMatrix(SymMat<Double>& alpha, NRVec<Double>& beta, 
								Double& chisq_all, Vec_DP& chisq, ostringstream& outPutString) const;
	
	virtual void setCovariantMatrix(const NRMat<Double>& Lmat, 
									const NRVec<Double>& adiag,
									const Double& eps1);
	
	virtual void setWeightMatrixOnIntensities(const SymMat<Double>& alpha);

	virtual ZErrorMessage addHistogram(const vector<T>& xdata,
										const Vec_DP& ydata,
										const Vec_DP& inv_wty,
										const Vec_BOOL& nfit);
};


template<class T>
inline Int4 MarquardtFmodelBase<T>::putDegreeOfFreedom(const Int4& index) const
{
	assert( 0 <= index && index < (Int4)m_xdata.size() );
	assert( m_degf[index] >= 0 );
	return m_degf[index];
}


template<class T>
MarquardtFmodelBase<T>::MarquardtFmodelBase(IMarquardtFmodel<T>& func)
{
//    m_weight_flag = true;          // 0:weight proportional(not equal) to standard error.

    m_func = &func;

	// Constraints(a Jacobian matrix) J.
	m_Cmat.resize(m_func->putParamNum());
	m_func->putConstraint(&m_Cmat[0]);
	
	const Int4 profile_num = m_func->putProfileNum();
	m_degf.resize(profile_num);
	m_num_indep_param.resize(profile_num, 0);

	m_num_indep_param_all = m_func->putNumberOfIndependentParamAll();
	m_degf_all = - m_num_indep_param_all;
	m_num_data_all = 0;
	for(Int4 l=0; l<profile_num; l++)
	{
		m_func->setTemporaryProfile(l);
		m_num_indep_param[l] = m_func->putNumberOfIndependentParam();
		m_degf[l] = - m_num_indep_param[l];
	}
}


template<class T>
ZErrorMessage MarquardtFmodelBase<T>::addHistogram(const vector<T>& xdata, const Vec_DP& ydata,
		const Vec_DP& edata, const Vec_BOOL& nfit)
{
	const Int4 index = m_xdata.size();
	assert( index < m_func->putProfileNum() );

	m_xdata.resize( index + 1 );
	m_ydata.resize( index + 1 );
	m_inv_edata.resize( index + 1 );
	
   	const Int4 isize = xdata.size();
	for(int k=0; k<isize; k++)
	{
		if( !nfit[k] ) continue;	// Excluded range;
		if( edata[k] <= 0.0 ) return ZErrorMessage("SOME ERROR VALUE <= 0 IN HISTOGRAM", __FILE__, __LINE__, __FUNCTION__);

		m_xdata[index].push_back( xdata[k] );
		m_ydata[index].push_back( ydata[k] );
		m_inv_edata[index].push_back( 1.0 / edata[k] );
   	}
	
	// Degree of Freedom.
	const Int4 isize2 = m_xdata[index].size();
//	m_num_data[index] = isize2;
	m_num_data_all += isize2;
	m_degf[index] += isize2;
	m_degf_all += isize2;
	if( m_degf[index] <= 0 ) return ZErrorMessage("THE NUMBER OF DATA IS NOT MORE THAN THE PARAMETERS", __FILE__, __LINE__, __FUNCTION__);
	if( index + 1 >= m_func->putProfileNum() )
	{
		if( m_degf_all <= 0 ) return ZErrorMessage("THE NUMBER OF DATA IS NOT MORE THAN THE PARAMETERS", __FILE__, __LINE__, __FUNCTION__);
	}
	
	return ZErrorMessage();
}


template<class T>
Int4 MarquardtFmodelBase<T>::setParam(const NRVec<Double>& dpara, 
		NRVec<Double>& paratry, ostringstream& outPutString)
{ 
	const Int4 ParaNum = m_func->putParamNum();
    	
	Int4 d2=0;
    for (Int4 d=0; d<ParaNum; d++)
    {
       	if(m_Cmat[d].ID==_IDVary) paratry[d] += dpara[d2++];
	   	else if(m_Cmat[d].ID==_IDDepend) paratry[d] += product_constr(m_Cmat[d].constr, &dpara[0]);
    }
	
	return setParam(paratry, outPutString);
}


template<class T>
bool MarquardtFmodelBase<T>::calJacobian(
		NRMat<Double>& Jacob, NRVec<Double>& R, Double& chisq_all, Vec_DP& chisq) const
{
	assert( Jacob.nrows() == m_num_data_all );
	assert( R.size() == m_num_data_all );
	
	const Int4& fitparam_num = m_func->putNumberOfIndependentParamAll();
	const Int4& profile_num = m_func->putProfileNum();
	assert( fitparam_num == Jacob.ncols() );

   	chisq_all=0.0;
   	
   	chisq.clear();
   	chisq.resize(profile_num, 0.0);
   	
   	Double ymod;
   	
	for (Int4 l=0,i2=0; l<profile_num; l++)
	{
		m_func->setTemporaryProfile(l);

		const vector<T>& xdata = m_xdata[l];
		const Vec_DP& ydata = m_ydata[l];
		const Vec_DP& inv_edata = m_inv_edata[l];
		const Int4 data_num = xdata.size();

		for (Int4 i=0; i<data_num; i++, i2++)
		{
			if( !m_func->putFittingCoefErrOut(xdata[i], ymod, Jacob[i2]) ) return false;

			R[i2] = (ydata[i] - ymod)*inv_edata[i];
	   		for (int j=0;j<fitparam_num;j++) Jacob[i2][j] *= inv_edata[i];
	   		chisq[l] += R[i2]*R[i2];
	   	}

		// Calculate chisq_all;
		chisq_all += chisq[l];
		chisq[l] /= m_degf[l];
	}
	chisq_all /= m_degf_all;

#ifndef __GNUC__
	if( !_finite(chisq_all) ) return false;
#else
	if( !finite(chisq_all) ) return false;
#endif
	return true;
}


//bool MarquardtFmodelBase::calJacobianWithWeight(IMarquardtFmodel& func,
//		NRMat<Double>& Jacob, NRVec<Double>& R, Double& chisq_all, Vec_DP& chisq)
//{
//	assert( Jacob.nrows() == m_num_data_all );
//	assert( Jacob.ncols() == m_num_indep_param_all );
//	assert( R.size() == m_num_data_all );
//	
//	const Int4& fitparam_num = m_num_indep_param_all;
//	const Int4& profile_num = m_func->putProfileNum();
//
//   	Double t;
//   	Double ymod, dy, wt;
//   	
//   	chisq_all=0.0;
//
//   	chisq.clear();
//   	chisq.resize(profile_num, 0.0);
//	
//	vector<double> dyda(fitparam_num);
//	for (Int4 l=0, i2=0; l<profile_num; l++)
//	{
//		func.setTemporaryProfile(l);
//
//		const Int4& data_num = m_num_data[l];
//		const Vec_DP& xdata = m_xdata[l];
//		const Vec_DP& ydata = m_ydata[l];
//		const Vec_DP& inv_edata = m_inv_edata[l];
//
//		const Int4& start = i2;
//		for (Int4 i=0; i<data_num; i++, i2++)
//		{
//			if( !func.putFittingCoefErrOut(xdata[i], ymod, Jacob[i2]) ) return false;
//
//	   		R[i] = (ydata[i] - ymod)*inv_edata[i];
//	   		for (int j=0;j<fitparam_num;j++) Jacob[i2][j] *= inv_edata[i];
//	   		chisq[l] += R[i]*R[i];
//	   	}
//
//   		t = sqrt( m_num_data[l] / chisq[l] );
//   		chisq_all += chisq[l];
//   		chisq[l] /= m_degf[l]; 
//
//   		i2=start;
//   		for(int i=0; i<data_num; i++, i2++)
//		{
//		   	R[i] *= t; 
//	   		for (int j=0;j<fitparam_num;j++) Jacob[i2][j] *= t;
//		}
//	}
//	chisq_all /= m_degf_all;
//
//	return true;
//}


template<class T>
void MarquardtFmodelBase<T>::calLMMatrixJ(const NRMat<Double>& Jacob, const NRVec<Double>& R,
		SymMat<Double>& alpha, NRVec<Double>& beta)
{
	assert( R.size() == Jacob.nrows() );
	
	const Int4& num_data_all = Jacob.nrows();
	const Int4& fitparam_num = Jacob.ncols();

	assert( alpha.size() >= fitparam_num );
	assert( beta.size() >= fitparam_num );
	
   	alpha = 0.0;
   	beta = 0.0;
	
	for (Int4 i=0; i<num_data_all;i++)
	{
		const Double& dy = R[i];
		const Double* dyda = Jacob[i];
   		for (int j=0;j<fitparam_num;j++)
   		{
       		for (int k=0; k<j+1;k++) alpha(j,k) += dyda[j]*dyda[k];
       		beta[j] += dy*dyda[j];
  		}
   	}
}


//  For the calculation when m_weight_flag is true.
//  allparam_num : the number of all the parameters.
//  fitparam_num : the number of all the fitting parameters.
//  Evaluate the linearized fitting matrix and calculate chi-square.
//  The solution of the linear equation alpha*x = beta, gives approximated values of parameters.
template<class T>
bool MarquardtFmodelBase<T>::calLMMatrix(SymMat<Double>& alpha, NRVec<Double>& beta,
		Double& chisq_all, Vec_DP& chisq, ostringstream& outPutString) const
{
	const Int4& fitparam_num = beta.size();
	assert( fitparam_num == alpha.size() );

	const Int4& profile_num = m_func->putProfileNum();
   	alpha = 0.0;
   	beta = 0.0;
   	
   	chisq_all=0.0;
	chisq.clear();
	chisq.resize(profile_num, 0.0);
	
   	Double ymod, dy, wt;
   	
	vector<double> dyda(fitparam_num);
	for (Int4 l=0; l<profile_num; l++)
	{
		m_func->setTemporaryProfile(l);

		const vector<T>& xdata = m_xdata[l];
		const Vec_DP& ydata = m_ydata[l];
		const Vec_DP& inv_edata = m_inv_edata[l];
		const Int4& data_num = xdata.size();

		for (Int4 i=0; i<data_num;i++)
		{
			if( !m_func->putFittingCoefErrOut(xdata[i], ymod, &dyda[0], outPutString) ) return false;

			const Double wtdata = inv_edata[i]*inv_edata[i];
			dy = ydata[i] - ymod;
	   		for (int j=0;j<fitparam_num;j++)
	   		{
	       		wt = dyda[j]*wtdata;
	       		for (int k=0; k<j+1;k++) alpha(j,k) += wt*dyda[k];
	       		beta[j] += dy*wt;
	  		}
	   		chisq[l] += dy*dy*wtdata;
	   	}

		// Calculate chisq_all;
		chisq_all += chisq[l];
		chisq[l] /= m_degf[l];
	}
	chisq_all /= m_degf_all;
#ifndef __GNUC__
	if( !_finite(chisq_all) ) return false;
#else
	if( !finite(chisq_all) ) return false;
#endif
	return true;
}




template<class T>
void MarquardtFmodelBase<T>::setCovariantMatrix(const NRMat<Double>& Lmat, const NRVec<Double>& adiag,
		const Double& eps1)
{
	SymMat<Double> covar(m_num_indep_param_all, m_num_indep_param_all);

	getInverseMatrix(Lmat, adiag, 0, m_num_indep_param_all, eps1, covar);
	
	m_func->setCovariantMatrix(covar);
}


template<class T>
void MarquardtFmodelBase<T>::setWeightMatrixOnIntensities(const SymMat<Double>& alpha)
{
	m_func->setWeightMatrixOnIntensities(alpha);
}



inline void copy(const SymMat<Double>& lhs, NRMat<Double>& rhs,
		const Int4& ibegin, const Int4& iend)
{
	assert( rhs.nrows() == lhs.size() && rhs.ncols() == lhs.size() );
	
	for(Int4 i=ibegin; i<iend; i++)
	{
		for(Int4 j=ibegin; j<iend; j++)
		{
			rhs[i][j] = lhs(i,j);
		}
	}
}


void calculate_dparam(const NRVec<Double>& adiag, const NRVec<Double>& beta, NRVec<Double>& dbeta,
		const Double& alamda, const Double& eps1,
		const Int4& ibegin, const Int4& iend);


//inline Double norm2(const NRVec<Double>& adiag, const NRVec<Double>& vec)
//{
//	const Int4 isize = vec.size();
//	assert( isize > 0 && isize == adiag.size() );
//	
//	Double ans = 0.0;
//	for(Int4 i=0; i<isize; i++)
//	{
//		ans += vec[i]*vec[i]*adiag[i];
//	}
//	return ans;
//}




inline void product_right(NRVec<Double>& beta, const NRMat<Double>& Lmat,
		const Int4& ibegin, const Int4& iend)
{
	assert( beta.size() == Lmat.nrows() );
	assert( beta.size() == Lmat.ncols() );
	assert( 0 <= ibegin && iend <= beta.size() );

	NRVec<Double> ans = beta;


	for(Int4 k=ibegin; k<iend; k++) 
	{
		beta[k] = 0.0;
		for(Int4 j=ibegin; j<iend; j++)
			beta[k] += ans[j]*Lmat[j][k];
	}
}

inline void product_left(const NRMat<Double>& Lmat, const NRVec<Double>& beta,
		NRVec<Double>& ans, const Int4& ibegin, const Int4& iend)
{
	assert( beta.size() == Lmat.nrows() && beta.size() == Lmat.ncols() );
	assert( beta.size() == ans.size() );
	assert( 0 <= ibegin && iend <= beta.size() );

	for(Int4 k=ibegin; k<iend; k++)
	{
		ans[k] = 0.0;
		for(Int4 j=ibegin; j<iend; j++) ans[k] += Lmat[k][j]*beta[j];
	}
}

void check_diagonal(const NRVec<Double>& adiag, const Double& eps1);

#endif /*MARQUARDTFMODELBASE_HH_*/
