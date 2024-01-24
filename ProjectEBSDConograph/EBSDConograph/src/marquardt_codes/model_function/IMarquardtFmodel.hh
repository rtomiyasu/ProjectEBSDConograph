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
#ifndef _IMarquardtFmodel_hh_
#define _IMarquardtFmodel_hh_
// IMarquardtFmodel.hh
#include <vector>
#include <typeinfo>
#include "../../conograph_codes/ZAnalysisTypes.hh"
#include "../../conograph_codes/utility_data_structure/SymMat.hh"
#include "../zerror_type/error_out.hh"
#include "covar_matrix.hh"

using namespace std;

// Interface class for LemarqMethod class.

template<class T>
class IMarquardtFmodel
{
protected:
//	static bool m_param_view_flag;

	// Set the values of parameters to fit by the argument.
	// If the parameters in the argument are out of range,
	//     returns 0, when this method quits to set them.
	//     returns 2, when this method moves the parameters into the range and sets them.
	// Otherwise returns 1.
	virtual Int4 setParam(Double*, ostringstream& outPutString) = 0;

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, also return the derivatives at X by dYdX.
	virtual bool putFittingCoef(const T& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const = 0;

public:
    virtual ~IMarquardtFmodel(){};
    
	// Return the number of all the parameters.(This includes the parameters fixed.)
	virtual Int4 putParamNum() const = 0;
	
	// Returns all the number of parameters independently fit.
    virtual Int4 putProfileNum() const{ return 1; };

    // Sets the temporary profile function.
	virtual void setTemporaryProfile(const Int4& iprof){ return; };

	// Returns all the number of parameters independently fit.
    virtual Int4 putNumberOfIndependentParamAll() const{ return putNumberOfIndependentParam(); };

    // Return the number of parameters independently fit.
    virtual Int4 putNumberOfIndependentParam() const;

    // Returns the flags whether the parameter is linear or not.
    virtual void putLinearParamFlag(bool* ans) const{ for(Int4 i=0; i<putParamNum(); i++) ans[i] = false; };

	// Returns constraints on the parameters.
	virtual void putConstraint(constr_DP* const) const = 0;

	// Sets constraints on the parameters.
	virtual void setConstraint(const constr_DP*) = 0;

	// Sets the covariant matrix. (The arugument is the the covariant matrix for the independently fit parameters.)
	virtual void setCovariantMatrix(const SymMat<Double>&); 

	// Sets the covariant matrix among linear parameters.
	virtual void setWeightMatrixOnIntensities(const SymMat<Double>& weight_mat) { throw ZErrorMessage(ZErrorEscapeFunction, __FILE__, __LINE__, __FUNCTION__); };

	// Sets the covariant matrix. (The arugument is the the covariant matrix for all the parameters.)
	virtual void setCovariantMatrixAll(const SymMat<Double>&) = 0;

	// Return the value(Ymod) of the model function at X.
	inline void putYmodel(const T& X, Double& Ymod) const;

	// Set the values of parameters by the argument.
	// If the parameters in the argument is out of range,
	//     returns 0, when this method quits to set them.
	//     returns 2, when this method moves the parameters into the range and sets them.
	// Otherwise returns 1.
	virtual Int4 setParamView(Double*, ostringstream& outPutString);

	inline Int4 setParamErrOut(Double*, ostringstream& outPutString, const bool flag = false);

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, the entries of dYdP are derivatives at X by the fitting parameters.
	inline bool putFittingCoefErrOut(const T& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;

// Approximate the x-value where y-value equals the argument Y near ini_x by Newton-Raphson method. 
//	bool putXValueForYValue(const Double& ini_x, const Double& goal_Y, 
//								Double& ans, const Double thred = 1.0e-8) const;
};

template<class T>
inline void IMarquardtFmodel<T>::putYmodel(const T& x, Double& Ymod) const
{
    putFittingCoef(x, Ymod, NULL);
}


// Count the number of the parameters independently fit.
inline Int4 countNumberOfIndependentParam(
const Mat_DP_constr::const_iterator& it_begin,const Mat_DP_constr::const_iterator& it_end)
{
	Mat_DP_constr::const_iterator it = it_begin;
	Int4 sum = 0;
	while( it != it_end )
	{
		if(it->ID==_IDVary) sum++;
		it++;
	}
	return sum;
}


template<class T>
inline Int4 IMarquardtFmodel<T>::setParamErrOut(Double* param, ostringstream& outPutString, const bool flag)
{
	if(flag) return setParamView(param, outPutString);
	else
	{
		Int4 ans = setParam(param, outPutString);
		if( ans == 0 ) outPutString << string( typeid(*this).name() ) + "::setParam failed.\n";
		else if( ans == 2 ) outPutString << string( typeid(*this).name() ) + "::setParam reports that some parameter is near the border.\n";
		
		return ans;
	}
}

template<class T>
inline bool IMarquardtFmodel<T>::putFittingCoefErrOut(const T& X, Double& Ymod, Double* dYdX, ostringstream& outPutString) const
{
	if( putFittingCoef(X, Ymod, dYdX, outPutString) ) return true;
	outPutString << typeid(*this).name() << "::putFittingCoef failed." << endl;
	
	return false;
}

template<class T>
Int4 IMarquardtFmodel<T>::putNumberOfIndependentParam() const
{
	Mat_DP_constr cmat( putParamNum() );
	putConstraint(&cmat[0]);
	return countNumberOfIndependentParam(cmat.begin(),cmat.end());
}


template<class T>
Int4 IMarquardtFmodel<T>::setParamView(Double* Param, ostringstream& outPutString)
{
	const Int4 param_num = putParamNum();
	outPutString << typeid(*this).name() << " : ";
	for(Int4 k=0; k<param_num; k++) outPutString << Param[k] << " ";
	outPutString << endl;
	
	return setParamErrOut(Param, outPutString);
}

// covar : the covariant matrix for the parameters independently fit.
// Set covariant matrix by the member mothod setCovariantMatrixAll.
template<class T>
void IMarquardtFmodel<T>::setCovariantMatrix(const SymMat<Double>& covar)
{
	const Int4 ma2 = putParamNum();
	Mat_DP_constr Cmat( ma2 );
	this->putConstraint( &Cmat[0] );
	
	SymMat<Double> covar_all(ma2);
	calPartOfCovariantMatrix(Cmat, covar, 0, ma2, covar_all);
	setCovariantMatrixAll(covar_all);
/*
cout << "Cmat:\n";
for(size_t i=0; i<Cmat.size(); i++)
{
	cout << Cmat[i].ID << " ";
	for(size_t k=0; k<Cmat[i].constr.size(); k++)
	{
		cout << "(" << Cmat[i].constr[k].index << ", " << Cmat[i].constr[k].element << ") ";
	}
	cout << endl;
}
*/
}


inline void putOrderOfIndependentParam(const Mat_DP_constr& cmat, Vec_INT& number_indep_param)
{
    number_indep_param.resize( cmat.size() );
    for(UInt4 i=0, i2=0; i<cmat.size(); i++)
		if(cmat[i].ID==_IDVary) number_indep_param[i]=i2++;
		else number_indep_param[i]=-1;
}


inline void setIndex(Mat_DP_constr& Cmat)
{
	Vec_INT indep_index;
	putOrderOfIndependentParam(Cmat, indep_index);
	
	const Int4 ParamNum = Cmat.size();
	Vec_DP_save::iterator it2;
	Vec_DP_save tray;

	// Resets indices in Cmat using indep_index.
	for(Int4 i=0; i<ParamNum; i++)
	{
		if( Cmat[i].ID < _IDDepend ) continue;

		it2 = Cmat[i].constr.begin();
		while( it2 != Cmat[i].constr.end() )
		{
			if( indep_index[ it2->index ] >= 0 )
			{
				it2->index = indep_index[ it2->index ];
			}
			else if( Cmat[ it2->index ].ID == _IDFixed )
			{
				it2 = Cmat[i].constr.erase(it2);
				continue;
			}
			else // if( Cmat[ it2->index ].flag == 2 )
			{
				assert(false);
//				throw ZErrorMessage( "The "+num2str(i+1)+"\'th parameters is constrained by other constrained parameters",
//								__FILE__, __LINE__, __FUNCTION__);
			}
			it2++;
		}
	}
}


#endif
