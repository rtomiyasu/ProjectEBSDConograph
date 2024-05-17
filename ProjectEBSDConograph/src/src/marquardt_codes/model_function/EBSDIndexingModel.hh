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
#ifndef _EBSDIndexingModel_hh_
#define _EBSDIndexingModel_hh_

// EBSDIndexingModel.hh

#include "../zparam/FittingParameter.hh"
#include "../../conograph_codes/bravais_type/BravaisType.hh"
#include "../../conograph_codes/utility_data_structure/VecDat3.hh"
#include "IMarquardtFmodel.hh"


class EBSDIndexingModel : public IMarquardtFmodel< pair< VecDat3<Int4>, char > >
{
private:
	enum{ ParaNum = 21, MainParaNum = 12 };
//	static bool m_param_view_flag;

	BravaisType m_bravais_type;

	// Fitting parameters:
	// (Main parameters)
	// m_param[0]: Delta_x,
	// m_param[1]: Delta_y,
	// m_param[2]: Delta_z,
	// m_param[3]: a21,
	// m_param[4]: a22,
	// m_param[5]: a31,
	// m_param[6]: a32,
	// m_param[7]: a33,
	// m_param[8]: theta (rad),
	// m_param[9]: sigma (rad),
	// m_param[10]: phi (rad).
	// m_param[11]: Scale factor.
	// ( cos(theta) sin(theta)  0)   (1          0          0 )   ( cos(phi) sin(phi)  0)
	// (-sin(theta) cos(theta)  0) * (0  cos(sigma) sin(sigma)) * (-sin(phi) cos(phi)  0)
	// (        0           0   1)   (0 -sin(sigma) cos(sigma))   (        0       0   1)
	// The direction of (h k l) = (0 0 1) is set to (sin(sigma)*sin(phi), -sin(sigma)*cos(phi), cos(sigma)).
	vector<FittingParameter> m_param;

	// (The other fitting parameters)
	NRMat<FittingParameter> m_reciprocal_lattice_basis; // The entries are fitting parameters dependent on the above.
	SymMat<double> m_S_covar; // Covariance matrix between the entries refined_basis*transpose(refined_basis).

	double m_wavelength;

	Int4 m_num_indep;
	Vec_INT m_indep_index;	// Index of independently fit parameters.
	Mat_DP_constr m_cmat;

	bool putFittingCoef_X_Z(const VecDat3<Int4>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;
	bool putFittingCoef_Y_Z(const VecDat3<Int4>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;
	bool putFittingCoef_BandWidth(const VecDat3<Int4>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;

	static void putCellOrientation(const NRMat<double>& reciprocal_lattice_basis, NRMat<double>& basis_lower_triangle, double& theta, double& sigma, double& phi);

protected:
    void setBravaisType(const BravaisType& arg){ m_bravais_type = arg; };

    Int4 setParam(Double*, ostringstream& outPutString);

	// Return the value(Ymod) of the model function at X.
	// If dYdX is not NULL, the entries of dYdP are the derivatives at X by the fitting parameters.
	// hkl[i].second == false -> phi^{obs}
	// hkl[i].second == true -> phi^{obs}
	bool putFittingCoef(const pair<VecDat3<Int4>, char>& X, Double& Ymod, Double* dYdP, ostringstream& outPutString) const;

public:
	EBSDIndexingModel();
	virtual ~EBSDIndexingModel();
	
	inline void setWaveLength(const double& arg){ m_wavelength = arg; };

	static bool putInitialParameters(const BravaisType& bravais_type, const FittingParameter* ProjectionCenterShift, const FittingParameter& ScaleFactor,
									 const NRMat<double>& reciprocal_lattice_basis, vector<double>& init_param);

	Int4 putParamNum() const;
	
	void putResult(FittingParameter* ProjectionCenterShift, FittingParameter* EulerAngles, FittingParameter& ScaleFactor, NRMat<FittingParameter>& ans, SymMat<Double>& S_covar) const;
	
	void setCovariantMatrixAll(const SymMat<Double>&); 

    // Return the number of parameters independently fit.
    Int4 putNumberOfIndependentParam() const;

	// Returns constraints on the parameters.
	void putConstraint(constr_DP* const) const;

	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const constr_DP*);

	// Set the constraints which indicates which parameters are fixed or independently fit.
	void setConstraint(const BravaisType& arg, const bool fitProjectionCenterShiftShift[3], const bool& fitScaleFactor);

	pair<bool, ZErrorMessage> setFittedParam(const vector< pair<VecDat3<Int4>, char> >& hkl,
						const vector<Double>& lclparam,
						const vector<Double>& lclvar,
						const vector<bool>& nxfit,
						const bool& output_view_flag,
						const Double&, const Int4&, const vector<Double>&, Double& chisq_init, Double& chisq_all);
};

#endif
