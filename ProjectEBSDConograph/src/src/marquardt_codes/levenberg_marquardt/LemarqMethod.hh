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
#ifndef _LemarqMethod_h_
#define _LemarqMethod_h_
// LemarqMethod.hh

#ifndef __GNUC__
	#include <float.h>
#endif
#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../../conograph_codes/utility_data_structure/nrutil_nr.hh"
#include "../../zlog/zlog.hh"
#include "../zparam/FittingParameter.hh"
#include "../model_function/MarquardtFmodelBase.hh"
#include "stopx.hh"

using namespace std;

class LemarqMethod
{
private:
	static const Double EPS;
	static const Double MAX_LAMBDA;

	bool m_output_view_flag;

	// The method to calculate inverse matrix 
	// (Cholesky dcmp : 0, Singular Value dcmp : 1)
//	bool m_cal_inverse_by_SVD_flag;
	
	// Calculate the weight matrix of intenities?
	// (No : 0, Yes : 1)
//	bool m_cal_intensity_weight_flag;
	
	Double m_initial_lambda;			// Initial value of lambda.
//	Int4 m_method_to_calculate_lambda;	// 0:original method by Marquardt, 1:Method by Flecher, 2:Trust region method.
//	Double m_ita;

    Double m_limiter;       // Convergence judgment value.
    Int4 m_max_itnum;       // Iteration number.

	Int4 execute_Marquardt_Original(const Double& thred,
									const Double& chisq0_all, const Double& chisq_all,
									Double& alamda
								) const;

public:
	LemarqMethod();
	~LemarqMethod();
	
	// Set the member variables.
	void setParam(const bool& output_view_flag,
					const Double& limiter, const Int4& itnum,
					const Double& initial_lambda);
	
    // reason for terminating: 1 - stopped by small gradient J^T e
	//                         2 - stopped by small Dp
	//                         3 - stopped by itmax
	//                         4 - no further error reduction is possible. Restart with increased mu
	//                         5 - stopped by small chi-square.
	template<class T>
	pair<bool, ZErrorMessage> execute(const vector<Double>& i_paratry,
					MarquardtFmodelBase<T>& fmod, Int4& itnum, Double& chisq0_init, Double& chisq_all,
					Vec_DP& chisq, Int4& reason_for_terminating) const;
};

inline Double norm2(const NRVec<Double>& adiag, const NRVec<Double>& vec,
		const Int4& ibegin, const Int4& iend)
{
	if( ibegin >= iend ) return 0.0;
	assert( vec.size() == adiag.size() );
	assert( 0 <= ibegin && iend <= vec.size() );
	
	Double ans = 0.0;
	for(Int4 i=ibegin; i<iend; i++)
	{
		ans += vec[i]*vec[i]*adiag[i];
	}
	return ans;
}


inline Double norm2(const NRMat<Double>& L, const NRVec<Double>& dpara,
	const Int4& ibegin, const Int4& iend)
{
	assert( L.nrows() == L.ncols() && L.nrows() == dpara.size() );
	assert( 0 <= ibegin && iend <= L.nrows() );
	
	Double t, ans = 0.0;
	for(Int4 i=ibegin; i<iend; i++)
	{
		// t[i] = transpose(L[*][i])*dpara
		t = 0.0;
		for(Int4 j=i; j<iend; j++)
		{
			t += L[j][i]*dpara[j];
		}
		
		ans += t*t;
	}
	return ans;
}



// This method executes Levenberg-Marquardt method.
// nx, ny are the x-values and y-values of histograms to fit.
// nfit gives the excluded region.
// If nfit[*][i] is true, ny[*][i] is included in the region to fit.
// If nfit[*][i] is false, ny[*][i] is excluded from the region. 
// inv_wty[*][i] is the error value of ny[*][i + irange_start[[*]]]
// irange_start and irange_end gives the range to fit.
// The elements in iparatry are the initial parameters p0_1, ..., p0_m.
// fmod is the model function.
// chisq : Chi-square value as a result of chi-square fitting.
// degf : Degree of freedom of chi-square.
// itnum : Iteration number
// inv_wty is used as weights in least square method. 
template<class T>
pair<bool, ZErrorMessage> LemarqMethod::execute(const vector<Double>& i_paratry,
MarquardtFmodelBase<T>& fmod, Int4& itnum, Double& chisq0_init, Double& chisq0_all, Vec_DP& chisq0,
Int4& reason_terminate_LM) const
{
	reason_terminate_LM = 0;
	try{
    	// Degree of Freedom.
		const Int4 num_nonlinear_indep_param_all = fmod.putNumberOfNonlinearIndependentParamAll();
		const Int4 num_indep_param_all = fmod.putNumberOfIndependentParamAll();
		if( num_indep_param_all <= 0 )
		{
			return pair<bool, ZErrorMessage>(false, ZErrorMessage("NUMBER OF DATA IS TOO SMALL", __FILE__, __LINE__, __FUNCTION__));
		}

		const Double eps1 = num_indep_param_all * EPS;
		
		SymMat<Double> alpha(num_indep_param_all);
    	NRVec<Double> beta(num_indep_param_all);
    	NRVec<Double> adiag(num_indep_param_all), dbeta(num_indep_param_all);
       	NRMat<Double> Lmat(num_indep_param_all, num_indep_param_all);

		const Int4 ma2 = fmod.putParamNum();
		assert( ma2 == (Int4)i_paratry.size() );

		NRVec<Double> paratry0( ma2 );	// initial parameters p0_1, ..., p0_m.
		for(Int4 l=0; l<ma2; l++) paratry0[l] = i_paratry[l];

		ostringstream outPutString;
		if( fmod.setParam(paratry0, outPutString) <= 0 )
		{
			return pair<bool, ZErrorMessage>(false, ZErrorMessage("THE PROGRAM COULD NOT SET WRONG INITIAL PARAMETER", __FILE__, __LINE__, __FUNCTION__));
		}

if( m_output_view_flag ) ZLOG_INFO( outPutString.str() );
outPutString.str("");
outPutString.clear();

		// Calculate the matrix for least squares.
		if( fmod.putSeparableLSFlag() )
		{
			if( !fmod.optimizeLinearParam(paratry0, eps1, alpha, beta, Lmat, adiag,
											chisq0_all, chisq0) )
			{
				return pair<bool, ZErrorMessage>(false, ZErrorMessage("WRONG INITIAL PARAMETER", __FILE__, __LINE__, __FUNCTION__));
			}
		}
		else
		{
			if( !fmod.calLMMatrix(alpha, beta, chisq0_all, chisq0, outPutString) )
			{
				return pair<bool, ZErrorMessage>(false, ZErrorMessage("WRONG INITIAL PARAMETER", __FILE__, __LINE__, __FUNCTION__));
			}
		}
		chisq0_init = chisq0_all;

if( m_output_view_flag ) ZLOG_INFO( outPutString.str() );
outPutString.str("");
outPutString.clear();

#ifndef __GNUC__
		if( !_finite(chisq0_all) )
#else
		if( !finite(chisq0_all) )
#endif
		{
			return pair<bool, ZErrorMessage>(false, ZErrorMessage("Chi-square diverged", __FILE__, __LINE__, __FUNCTION__));
		}
		if( chisq0_all < eps1 ) reason_terminate_LM = 5;
		

		copy(alpha, Lmat, 0, num_nonlinear_indep_param_all);
		SVdcmp_precondition(Lmat, adiag, 0, num_nonlinear_indep_param_all);

		product_right(beta, Lmat, 0, num_nonlinear_indep_param_all);	// beta = transpose(V^(-1)) * beta
		if( fmod.putSeparableLSFlag() ) product_right(beta, Lmat, num_nonlinear_indep_param_all, num_indep_param_all);	// beta = transpose(V^(-1)) * beta
//   		check_diagonal(adiag, eps1);

		itnum = 0;

		NRVec<Double> paratry;
		Double chisq_all;
		Vec_DP chisq;

		SymMat<Double> Mtemp(num_indep_param_all);
		NRVec<Double> Vtemp(num_indep_param_all);
		NRMat<Double> Ltemp(num_indep_param_all, num_indep_param_all);
		NRVec<Double> Dtemp(num_indep_param_all);
		NRVec<Double> dpara(num_indep_param_all);

		bool set_flag = true; // Flag that indicates whether the parameters in paratray0 are set in fmod.
		Double alamda = m_initial_lambda;

		while ( !reason_terminate_LM )
		{
//			if( putInterruptionSignal() )
//			{
//				reason_terminate_LM = 1;
//				break;
//			}

			if( itnum >= m_max_itnum )
			{
				reason_terminate_LM = 3;
				break;
			}

			itnum++;
//if( m_output_view_flag )
//{
//	cout << "\n" << itnum << " : " << chisq0_all << endl;
//}

			// Calculate delta2 and differences between old parameters and new ones.
			calculate_dparam(adiag, beta, dbeta, alamda, eps1, 0, num_indep_param_all);

			// Get dpara = alpha * dbeta as the solution of the linear equation of the least-squares.
			product_left(Lmat, dbeta, dpara, 0, num_nonlinear_indep_param_all);
			if( fmod.putSeparableLSFlag() ) product_left(Lmat, dbeta, dpara, num_nonlinear_indep_param_all, num_indep_param_all);

			set_flag = false;

			paratry = paratry0;
			if( fmod.setParam(dpara, paratry, outPutString) )
			{
				bool flag = false;
				if( fmod.putSeparableLSFlag() )
				{
					flag = fmod.optimizeLinearParam(paratry, eps1, Mtemp, Vtemp, Ltemp, Dtemp,
														chisq_all, chisq);
				}
				else
				{
					flag = fmod.calLMMatrix(Mtemp,Vtemp,chisq_all,chisq, outPutString);
				}
if( m_output_view_flag ) ZLOG_INFO( outPutString.str() );
outPutString.str("");
outPutString.clear();

#ifndef __GNUC__
				if( !_finite(chisq_all) )
#else
				if( !finite(chisq_all) )
#endif
				{
					if( m_output_view_flag ) ZLOG_INFO( "Chi-square diverged.\n" );
				}
				else if( flag )
				{
				    // reason for terminating: 0 - The argument paratry0 are replaced by new parameters.
					//                         1 - The new parameters are removed and the old ones in argument paratry0 are reserved.
					//						   4 - stopped by too large lambda.
					//                         5 - stopped by small chi-square.
					if( chisq_all < eps1 )
					{
						reason_terminate_LM = 5;
						break;
					}
					else // if( m_method_to_calculate_lambda == 0 )	// Original Marquardt method.
					{
						reason_terminate_LM = execute_Marquardt_Original(eps1, chisq0_all, chisq_all, alamda);
					}

					if( alamda > MAX_LAMBDA )
					{
						reason_terminate_LM = 4;
						break;
					}
					else if( reason_terminate_LM == 1 )
					{
						reason_terminate_LM = 0;
						continue;
					}
					else // if( reason_terminate_LM == 0 )
	        		{
	        	   		if( norm2( adiag, dbeta, 0, num_indep_param_all) <= m_limiter * num_indep_param_all )
	        	   		{
	        	   			reason_terminate_LM = 2;
	        	   			break;
	        	   		}

	        	   		alpha = Mtemp;
	        			beta = Vtemp;
	        			Lmat = Ltemp;
	        			adiag = Dtemp;

	        	   		paratry0=paratry;
	        			set_flag = true;

	        	   		chisq0_all=chisq_all;
	        	   		chisq0 = chisq;

        				copy(alpha, Lmat, 0, num_nonlinear_indep_param_all);
        				SVdcmp_precondition(Lmat, adiag, 0, num_nonlinear_indep_param_all);

        				product_right(beta, Lmat, 0, num_nonlinear_indep_param_all);
        				if( fmod.putSeparableLSFlag() ) product_right(beta, Lmat, num_nonlinear_indep_param_all, num_indep_param_all);	// beta = transpose(V^(-1)) * beta
//       				check_diagonal(adiag, eps1);

	        	   		continue;
	        		}
				}
			}
			alamda = alamda*10.0;
			// The Marquart parameter lambda has deverged.  Quit the calculation.
			if( alamda > MAX_LAMBDA ) reason_terminate_LM = 4;
		}

//if( m_output_view_flag )
//{
//	assert( 1 <= reason_terminate_LM && reason_terminate_LM <= 6 );
//	cout << "Reason for terminating Levenberg-Marquardt method : ";
//	if( reason_terminate_LM <= 1 ) cout << "LM was interrupted by user's operation.\n";
//	else if( reason_terminate_LM <= 2 ) cout << "stopped by small Dp.\n";
//	else if( reason_terminate_LM <= 3 ) cout << "stopped by itmax.\n";
//	else if( reason_terminate_LM <= 4 ) cout << "stopped by too large lambda.\n";
//	else if( reason_terminate_LM <= 5 ) cout << "stopped by small chi-square.\n";
//	else cout << "stopped by failure of Cholesky decomposition.(This caused by some bad conditions in a non-positive definite matrix.)\n\n";
//}

		if( !set_flag || fmod.putSeparableLSFlag() ) fmod.setParamAll(paratry0, outPutString);
		fmod.setCovariantMatrix(Lmat, adiag, eps1);
	}
	catch(bad_alloc& ball)
	{
		return pair<bool, ZErrorMessage>(false, nerror(ball, __FILE__, __LINE__, __FUNCTION__));
	}
	return pair<bool, ZErrorMessage>(true, ZErrorMessage());
}

#endif
