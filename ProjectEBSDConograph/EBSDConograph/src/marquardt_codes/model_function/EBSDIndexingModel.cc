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
#include "EBSDIndexingModel.hh"

#include "../../conograph_codes/bravais_type/CrystalSystem.hh"
#include "../zparam/etype_ID.hh"
#include "../levenberg_marquardt/Choleskydcmp.hh"
#include "../levenberg_marquardt/LemarqMethod.hh"

namespace {
	template<class T>
	inline T square(const T& arg)
	{
		return arg * arg;
	}

	SymMat<double> putMetricTensor(const NRMat<double>& lhs)
	{
		const int irow = lhs.nrows();
		const int icol = lhs.ncols();
		SymMat<double> ans(irow, 0.0);
		for(int k=0; k<irow; k++)
		{
			for(int j=k; j<irow; j++)
			{
				for(int i=0; i<icol; i++) ans(k,j) += lhs[k][i]*lhs[j][i];
			}
		}
		return ans;
	}

	template<class T>
	static NRMat<T> putMatrix(
			const T& a11, const T& a12, const T& a13,
			const T& a21, const T& a22, const T& a23,
			const T& a31, const T& a32, const T& a33)
	{
		NRMat<T> ans(3,3);
		ans[0][0] = a11;
		ans[0][1] = a12;
		ans[0][2] = a13;
		ans[1][0] = a21;
		ans[1][1] = a22;
		ans[1][2] = a23;
		ans[2][0] = a31;
		ans[2][1] = a32;
		ans[2][2] = a33;
		return ans;
	}

	static void setOrthoMatrix(const Double& theta, const Double& sigma, const Double& phi, NRMat<double>& ortho)
	{
		// ( cos(theta) sin(theta)  0)   (1          0          0 )   ( cos(phi) sin(phi)  0)
		// (-sin(theta) cos(theta)  0) * (0  cos(sigma) sin(sigma)) * (-sin(phi) cos(phi)  0)
		// (        0           0   1)   (0 -sin(sigma) cos(sigma))   (        0       0   1)
		ortho = mprod( mprod( putMatrix(cos(theta),sin(theta),0.0,  -sin(theta),cos(theta),0.0,  0.0,0.0,1.0),
							  putMatrix(              1.0,0.0,0.0,   0.0,cos(sigma),sin(sigma),  0.0,-sin(sigma),cos(sigma)) ),
							  putMatrix(    cos(phi),sin(phi),0.0,      -sin(phi),cos(phi),0.0,  0.0,0.0,1.0) );
	}
}

EBSDIndexingModel::EBSDIndexingModel()
	: m_param(MainParaNum, 0.0), m_reciprocal_lattice_basis(3,3,0.0), m_S_covar(6, 0.0)
{
	m_cmat.resize(ParaNum);
	for(Int4 i=0; i<MainParaNum; i++) m_cmat[i].ID = _IDFixed;
	for(Int4 i=MainParaNum; i<ParaNum; i++) m_cmat[i].ID = _IDDepend;
	m_num_indep = 0;
}

EBSDIndexingModel::~EBSDIndexingModel()
{
}

void EBSDIndexingModel::putCellOrientation(const NRMat<double>& reciprocal_lattice_basis,
		NRMat<double>& basis_lower_triangle, double& theta, double& sigma, double& phi)
{
	const SymMat<double> S(putMetricTensor(reciprocal_lattice_basis));
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			basis_lower_triangle[i][j] = S(i,j);
		}
	}

	// Cholesky decomposition.
	Choleskydcmp(basis_lower_triangle, 0, 3);

	NRMat<double> Inv_L(3,3);	// = inverse of S_mat.
	getInverseMatrix(basis_lower_triangle, Inv_L);

	const NRMat<double> ortho( mprod(Inv_L, reciprocal_lattice_basis) );

	// Set theta, sigma, phi (rad) from the argument ortho (an orthogonal matrix).
	// ( cos(theta) sin(theta)  0)   (1          0          0 )   ( cos(phi) sin(phi)  0)
	// (-sin(theta) cos(theta)  0) * (0  cos(sigma) sin(sigma)) * (-sin(phi) cos(phi)  0)
	// (        0           0   1)   (0 -sin(sigma) cos(sigma))   (        0       0   1)
	// The third row = (sin(sigma)*sin(phi), -sin(sigma)*cos(phi), cos(sigma)).
	// The third column = (sin(sigma)*sin(theta), sin(sigma)*cos(theta), cos(sigma))^T.
	if ( ortho[2][0]*ortho[2][0] + ortho[2][1] * ortho[2][1] <= 0.0 ) // sin(sigma) = 0.0
	{
		// Only theta + phi can be determined.
		theta =  0.; // Set theta to 0.
		sigma  = 0.; // Set sigma to 0
		phi = atan2( ortho[0][1], ortho[0][0] );  // = phi
    }
	else
	{
		// sin(sigma) > 0.0 is assumed here.
		theta = atan2( ortho[0][2], ortho[1][2] );
		phi = atan2( ortho[2][0],-ortho[2][1] );
		sigma   = atan2( sqrt(ortho[2][0]*ortho[2][0] + ortho[2][1]*ortho[2][1]), ortho[2][2] );  // = sigma
	}

#ifdef DEBUG
	NRMat<double> ortho2(3,3);
	setOrthoMatrix(theta, sigma, phi, ortho2);
	if( fabs( ortho[0][0] - ortho2[0][0] ) > max(max(fabs(ortho[0][0]), fabs(ortho2[0][0]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[0][1] - ortho2[0][1] ) > max(max(fabs(ortho[0][1]), fabs(ortho2[0][1]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[0][2] - ortho2[0][2] ) > max(max(fabs(ortho[0][2]), fabs(ortho2[0][2]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[1][0] - ortho2[1][0] ) > max(max(fabs(ortho[1][0]), fabs(ortho2[1][0]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[1][1] - ortho2[1][1] ) > max(max(fabs(ortho[1][1]), fabs(ortho2[1][1]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[1][2] - ortho2[1][2] ) > max(max(fabs(ortho[1][2]), fabs(ortho2[1][2]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[2][0] - ortho2[2][0] ) > max(max(fabs(ortho[2][0]), fabs(ortho2[2][0]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[2][1] - ortho2[2][1] ) > max(max(fabs(ortho[2][1]), fabs(ortho2[2][1]))*1.0e-4, 1.0e-16)
		&& fabs( ortho[2][2] - ortho2[2][2] ) > max(max(fabs(ortho[2][2]), fabs(ortho2[2][2]))*1.0e-4, 1.0e-16) )
	{
ZLOG_INFO(num2str(ortho[0][0]) + "," + num2str(ortho[0][1]) + "," + num2str(ortho[0][2])+"\n");
ZLOG_INFO(num2str(ortho[1][0]) + "," + num2str(ortho[1][1]) + "," + num2str(ortho[1][2])+"\n");
ZLOG_INFO(num2str(ortho[2][0]) + "," + num2str(ortho[2][1]) + "," + num2str(ortho[2][2])+"\n\n");
ZLOG_INFO(num2str(ortho2[0][0]) + "," + num2str(ortho2[0][1]) + "," + num2str(ortho2[0][2])+"\n");
ZLOG_INFO(num2str(ortho2[1][0]) + "," + num2str(ortho2[1][1]) + "," + num2str(ortho2[1][2])+"\n");
ZLOG_INFO(num2str(ortho2[2][0]) + "," + num2str(ortho2[2][1]) + "," + num2str(ortho2[2][2])+"\n\n\n");

		assert( fabs( ortho[0][0] - ortho2[0][0] ) <= max(max(fabs(ortho[0][0]), fabs(ortho2[0][0]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[0][1] - ortho2[0][1] ) <= max(max(fabs(ortho[0][1]), fabs(ortho2[0][1]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[0][2] - ortho2[0][2] ) <= max(max(fabs(ortho[0][2]), fabs(ortho2[0][2]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[1][0] - ortho2[1][0] ) <= max(max(fabs(ortho[1][0]), fabs(ortho2[1][0]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[1][1] - ortho2[1][1] ) <= max(max(fabs(ortho[1][1]), fabs(ortho2[1][1]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[1][2] - ortho2[1][2] ) <= max(max(fabs(ortho[1][2]), fabs(ortho2[1][2]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[2][0] - ortho2[2][0] ) <= max(max(fabs(ortho[2][0]), fabs(ortho2[2][0]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[2][1] - ortho2[2][1] ) <= max(max(fabs(ortho[2][1]), fabs(ortho2[2][1]))*1.0e-4, 1.0e-16) );
		assert( fabs( ortho[2][2] - ortho2[2][2] ) <= max(max(fabs(ortho[2][2]), fabs(ortho2[2][2]))*1.0e-4, 1.0e-16) );
	}

#endif
}

bool EBSDIndexingModel::putInitialParameters(const BravaisType& arg,
		const FittingParameter* ProjectionCenterShift, const FittingParameter& ScaleFactor, const NRMat<double>& reciprocal_lattice_basis,
		vector<double>& init_param)
{
	assert( Determinant3(reciprocal_lattice_basis) > 0.0 );

	NRMat<double> basis_lower_triangle(3,3), ortho(3,3);
	init_param.resize(ParaNum);
	putCellOrientation(reciprocal_lattice_basis, basis_lower_triangle, init_param[8], init_param[9], init_param[10]);

	// Pattern center.
//	assert(ProjectionCenterShift.size());
	init_param[0] = ProjectionCenterShift[0].value;
	init_param[1] = ProjectionCenterShift[1].value;
	init_param[2] = ProjectionCenterShift[2].value;
	init_param[11] = ScaleFactor.value * basis_lower_triangle[0][0];

	// a21, a22, a31, a32, a33.
	init_param[3] = basis_lower_triangle[1][0] / basis_lower_triangle[0][0];
	init_param[4] = basis_lower_triangle[1][1] / basis_lower_triangle[0][0];
	init_param[5] = basis_lower_triangle[2][0] / basis_lower_triangle[0][0];
	init_param[6] = basis_lower_triangle[2][1] / basis_lower_triangle[0][0];
	init_param[7] = basis_lower_triangle[2][2] / basis_lower_triangle[0][0];

	if(arg.enumLatticeSymmetry() == C2h_X)
	{
		// a21, a31 = 0.
		init_param[3] = 0.0;
		init_param[5] = 0.0;
	}
	else if(arg.enumLatticeSymmetry() == C2h_Y)
	{
		// a21, a32 = 0.
		init_param[3] = 0.0;
		init_param[6] = 0.0;
	}
	else if(arg.enumLatticeSymmetry() == C2h_Z)
	{
		// a31, a32 = 0.
		init_param[5] = 0.0;
		init_param[6] = 0.0;
	}
	else if(arg.enumLatticeSymmetry() == D2h)
	{
		// a21, a31, a32 = 0.
		init_param[3] = 0.0;
		init_param[5] = 0.0;
		init_param[6] = 0.0;
	}
	else if(arg.enumLatticeSymmetry() == D4h_Z)
	{
		init_param[3] = 0.0;
		init_param[4] = 1.0;
		init_param[5] = 0.0;
		init_param[6] = 0.0;
		init_param[7] = basis_lower_triangle[2][2] * 2.0 / (basis_lower_triangle[0][0] + basis_lower_triangle[1][1]);
	}
	else if(arg.enumLatticeSymmetry() == D31d_rho)
	{
		const double a = (init_param[3] + init_param[5])*0.5;
		if( -0.5 >= a || a >= 1.0 ) return false;
		init_param[3] = a;
		init_param[4] = sqrt(1 - a*a);
		init_param[5] = a;
		init_param[6] = init_param[4]*a/(1+a);
		init_param[7] = sqrt((1-a)*(1+2*a)/(1+a));
	}
	else if(arg.enumLatticeSymmetry() == D3d_1_hex || arg.enumLatticeSymmetry() == D6h)
	{
		init_param[3] = 0.5;
		init_param[4] = sqrt(3.0)*0.5;
		init_param[5] = 0.0;
		init_param[6] = 0.0;
		init_param[7] = basis_lower_triangle[2][2] / basis_lower_triangle[0][0];
	}
	else if(arg.enumLatticeSymmetry() == Oh)
	{
		init_param[3] = 0.0;
		init_param[4] = 1.0;
		init_param[5] = 0.0;
		init_param[6] = 0.0;
		init_param[7] = 1.0;
	}

	return true;
}

Int4 EBSDIndexingModel::putParamNum() const
{ 
	return ParaNum;
};


void EBSDIndexingModel::putResult(FittingParameter* ProjectionCenterShift, FittingParameter* EulerAngles,
		FittingParameter& ScaleFactor, NRMat<FittingParameter>& ans, SymMat<Double>& S_covar) const
{
	ProjectionCenterShift[0] = m_param[0];
	ProjectionCenterShift[1] = m_param[1];
	ProjectionCenterShift[2] = m_param[2];

	EulerAngles[0] = m_param[8];
	EulerAngles[1] = m_param[9];
	EulerAngles[2] = m_param[10];

	ScaleFactor = m_param[11];

	ans = m_reciprocal_lattice_basis;

	S_covar = m_S_covar;
/*
cout << "Resulting parameters: \n";
for (int i=0; i<MainParaNum; i++){
	  cout << m_param[i].value << " ";
}
cout << endl;
*/
}

Int4 EBSDIndexingModel::putNumberOfIndependentParam() const
{
	return m_num_indep;
}


void EBSDIndexingModel::putConstraint(constr_DP* const tray) const
{
	for(Int4 k=0; k<ParaNum; k++)
	{
		tray[k] = m_cmat[k];
	}
}


void EBSDIndexingModel::setConstraint(const constr_DP* tray)
{
	for(Int4 k=0; k<ParaNum; k++)
	{
		m_cmat[k] = tray[k];
	}
	m_num_indep = countNumberOfIndependentParam(m_cmat.begin(), m_cmat.end());
}


void EBSDIndexingModel::setConstraint(const BravaisType& arg, const bool fitProjectionCenterShiftShift[3], const bool& fitScaleFactor)
{
	m_bravais_type = arg;

	m_cmat.clear();
	m_cmat.resize(ParaNum);

	for(Int4 i=0; i<MainParaNum; i++) m_cmat[i].ID = _IDVary;
	for(Int4 i=MainParaNum; i<ParaNum; i++) m_cmat[i].ID = _IDDepend;

	if( !fitProjectionCenterShiftShift[0] ) m_cmat[0].ID = _IDFixed;
	if( !fitProjectionCenterShiftShift[1] ) m_cmat[1].ID = _IDFixed;
	if( !fitProjectionCenterShiftShift[2]
		|| (!fitScaleFactor && arg.enumBravaisType() == BravaisType::Triclinic) ) m_cmat[2].ID = _IDFixed;
	if( !fitScaleFactor ) m_cmat[11].ID = _IDFixed;

	if( arg.enumLatticeSymmetry() == C2h_X ) // Monoclinic
	{
		// a21, a31 = 0.
		m_cmat[3].ID = _IDFixed;
		m_cmat[5].ID = _IDFixed;
	}
	if( arg.enumLatticeSymmetry() == C2h_Y ) // Monoclinic
	{
		// a21, a32 = 0.
		m_cmat[3].ID = _IDFixed;
		m_cmat[6].ID = _IDFixed;
	}
	if( arg.enumLatticeSymmetry() == C2h_Z ) // Monoclinic
	{
		// a31, a32 = 0.
		m_cmat[5].ID = _IDFixed;
		m_cmat[6].ID = _IDFixed;
	}
	else if( arg.enumLatticeSymmetry() == D2h ) // Orthorhombic
	{
		// a21, a31, a32 = 0.
		m_cmat[3].ID = _IDFixed;
		m_cmat[5].ID = _IDFixed;
		m_cmat[6].ID = _IDFixed;
	}
	else if( arg.enumLatticeSymmetry() == D4h_Z ) // Tetragonal
	{
		m_cmat[3].ID = _IDFixed;
		m_cmat[4].ID = _IDFixed;
		m_cmat[5].ID = _IDFixed;
		m_cmat[6].ID = _IDFixed;
	}
	else if( arg.enumLatticeSymmetry() == D31d_rho ) // Rhombohedral
	{
		m_cmat[4].ID = _IDDepend;
		m_cmat[5].ID = _IDDepend;
		m_cmat[6].ID = _IDDepend;
		m_cmat[7].ID = _IDDepend;

		m_cmat[4].constr.push_back(index_set<double>(3, 0.0)); // The differential is set in the member function setParam.
		m_cmat[5].constr.push_back(index_set<double>(3, 1.0)); // a31 = a21.
		m_cmat[6].constr.push_back(index_set<double>(3, 0.0)); // The differential is set in the member function setParam.
		m_cmat[7].constr.push_back(index_set<double>(3, 0.0)); // The differential is set in the member function setParam.
	}
	else if( arg.enumLatticeSymmetry() == D3d_1_hex || arg.enumLatticeSymmetry() == D6h ) // Hexagonal.
	{
		// a21: fixed to 1/2.
		// a22: fixed to 1.
		// a31, a32: fixed to 0,
		for(Int4 i=3; i<7; i++) m_cmat[i].ID = _IDFixed;
	}
	else if( arg.enumLatticeSymmetry() == Oh ) // Cubic
	{
		// a21, a31, a32: fixed to 0,
		// a22, a33: fixed to 1.
		for(Int4 i=3; i<8; i++) m_cmat[i].ID = _IDFixed;
	}

	m_indep_index.clear();
	m_indep_index.resize(ParaNum, -1);
	for(Int4 i=0, i2=0; i<MainParaNum; i++)
	{
		if( m_cmat[i].ID == _IDVary ) m_indep_index[i] = i2++;
	}
	m_num_indep = countNumberOfIndependentParam(m_cmat.begin(), m_cmat.end());

	if( arg.enumLatticeSymmetry() == D31d_rho ) // Rhombohedral
	{
		assert(m_cmat[3].ID == _IDVary);
		m_cmat[4].constr.front().index = m_indep_index[3];
		m_cmat[5].constr.front().index = m_indep_index[3];
		m_cmat[6].constr.front().index = m_indep_index[3];
		m_cmat[7].constr.front().index = m_indep_index[3];
	}
}


void EBSDIndexingModel::setCovariantMatrixAll(const SymMat<Double>& cov)
{
	assert( cov.size() == ParaNum );
	for(Int4 k=0; k<MainParaNum; k++)
	{
		m_param[k].error = sqrt_d( cov(k,k) );
	}

	for(int i=0,k=MainParaNum; i<3; i++)
	{
		for(int j=0; j<3; j++,k++)
		{
			m_reciprocal_lattice_basis[i][j].error = sqrt_d( cov(k, k) );
		}
	}

	SymMat<Double> covar_Amat(5);
	copyPartOfCovariantMatrix(cov, 3, 8, covar_Amat);

	const double& a21 = m_param[3].value;
	const double& a22 = m_param[4].value;
	const double& a31 = m_param[5].value;
	const double& a32 = m_param[6].value;
	const double& a33 = m_param[7].value;

	NRMat<double> M_sij_akl(6, 5, 0.);

	M_sij_akl[1][0] = 2. * a21; //=d (a21*a21 + a22*a22) / d a21
	M_sij_akl[1][1] = 2. * a22; //=d (a21*a21 + a22*a22) / d a22

	M_sij_akl[2][2] = 2. * a31; //=d (a31^2+a32^2+a33^2)/ d a31
	M_sij_akl[2][3] = 2. * a32; //=d (a31^2+a32^2+a33^2)/ d a32
	M_sij_akl[2][4] = 2. * a33; //=d (a31^2+a32^2+a33^2)/ d a33

	M_sij_akl[3][0] = a31; //=d (a21*a31+ a22*a32)/ d a21
	M_sij_akl[3][1] = a32; //=d (a21*a31+ a22*a32)/ d a22
	M_sij_akl[3][2] = a21; //=d (a21*a31+ a22*a32)/ d a31
	M_sij_akl[3][3] = a22; //=d (a21*a31+ a22*a32)/ d a32

	M_sij_akl[4][0] = 0.; //=d a31/ d a21
	M_sij_akl[4][2] = 1.; //=d a31/ d a31

	M_sij_akl[5][0] = 1.;  //=d a21/ d a21

	m_S_covar = transform_sym_matrix(M_sij_akl, covar_Amat);
/*
	//m_S_covar = covar_Amat;  //0.0;
	cout<<"--cov --:"<<endl;
	for (int i=0; i<cov.size(); i++)
	{
		for (int j=0; j<cov.size(); j++) {
			cout<< cov(i,j) <<" ";
		}
		cout<<endl;
	}
	cout<<"--m_S_covar --:"<<endl;
	for (int i=0; i<m_S_covar.size(); i++)
	{
		for (int j=0; j<m_S_covar.size(); j++) {
			cout<< m_S_covar(i,j) <<" ";
		}
		cout<<endl;
	}
*/
}


Int4 EBSDIndexingModel::setParam(Double* param, ostringstream& outPutString)
{
//outPutString << "Parameters: " << endl;
//for (Int4 k=0; k<MainParaNum; k++) outPutString << param[k] << " ";
//outPutString << endl;

	if( param[2] <= -1.0 ) return 0;
	if( param[4] <= 0.0 || 1.0e+7 < param[3]*param[3]+param[4]*param[4] ) return 0;
	if( param[7] <= 0.0 || 1.0e+7 < param[5]*param[5]+param[6]*param[6]+param[7]*param[7] ) return 0;

	const Int4 ans = (param[11] < 0.0?2:1);
	if( param[11] < 0.0 )
	{
		param[11] = 0.0;
	}

	for(Int4 k=0; k<MainParaNum; k++)
	{
		m_param[k].value = param[k];
	}

	if(m_bravais_type.enumLatticeSymmetry() == D31d_rho) // Rhombohedral.
	{
		const double a = param[3];
		if(-0.5 >= a ) return 0;
		if( a >= 1.0 ) return 0;
		m_param[3].value = a;
		m_param[4].value = sqrt(1 - a*a);
		m_param[5].value = a;
		m_param[6].value = m_param[4].value*a/(1+a);
		m_param[7].value = sqrt((1-a)*(1+2*a)/(1+a));

		assert( (Int4)m_cmat[4].constr.begin()->index == m_indep_index[3] );
		assert( (Int4)m_cmat[6].constr.begin()->index == m_indep_index[3] );
		assert( (Int4)m_cmat[7].constr.begin()->index == m_indep_index[3] );

		// To do
		const double onePaPo3 = (1 + a) * (1 + a) * (1 + a) * (1 - a);
		m_cmat[4].constr.begin()->element = -a / sqrt (1 - a*a); // Set da22/da21.
		m_cmat[6].constr.begin()->element = (1 - a - a*a) / sqrt(onePaPo3); // Set da32/da21.
		m_cmat[7].constr.begin()->element = a*(2 + a) / sqrt (onePaPo3 * (1 + 2 * a)); // Set da33/da21.
	}

	NRMat<double> A_matrix(3, 3);
	A_matrix[0][0]= 1.;
	A_matrix[0][1]= 0.;
	A_matrix[0][2]= 0.;

	A_matrix[1][0]= m_param[3].value;
	A_matrix[1][1]= m_param[4].value;
	A_matrix[1][2]= 0.;

	A_matrix[2][0]= m_param[5].value;
	A_matrix[2][1]= m_param[6].value;
	A_matrix[2][2]= m_param[7].value;

	// Set the entries in m_reciprocal_lattice_basis by
	//                              (         1         0           0)
	// m_reciprocal_lattice_basis = (m_param[3] m_param[4]          0) * g(m_param[8], m_param[9], m_param[10])
	//                              (m_param[5] m_param[6] m_param[7])
	NRMat<double> ortho(3, 3);
	setOrthoMatrix(param[8], param[9], param[10], ortho);

	NRMat<double> AMultipg = mprod(A_matrix, ortho);

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++) m_reciprocal_lattice_basis[i][j].value = AMultipg[i][j];
	}


	for(int k=MainParaNum; k<ParaNum; k++)
	{
		m_cmat[k].constr.clear();
	}

	for(size_t i2=1, k2=3; i2<3; i2++)
	{
		for(size_t j2=0; j2<=i2; j2++, k2++)
		{
			NRMat<double> temp_dc_dp(3, 3, 0.);
			temp_dc_dp[i2][j2] = 1.; // temp_dc_dp = dA / d m_param[k2] (k2=3--7).

			const NRMat<double> dc_dA = mprod(temp_dc_dp, ortho);

			// Set d m_reciprocal_lattice_basis[i][j] / d m_param[k2] (k2=3--7) in dc_dp[i][j].
			for(int i=0, k=MainParaNum; i<3; i++)
			{
				for(int j=0; j<3; j++, k++)
				{
					if( m_cmat[k2].ID == _IDVary )
					{
						m_cmat[k].constr.push_back(index_set<double>(m_indep_index[k2], dc_dA[i][j]));
					}
					else if( m_cmat[k2].ID == _IDDepend )
					{
						for(Vec_DP_save::const_iterator it=m_cmat[k2].constr.begin(); it<m_cmat[k2].constr.end(); it++)
						{
//							assert( m_cmat[it->index].ID == _IDVary );
							m_cmat[k].constr.push_back(index_set<double>(it->index, it->element*dc_dA[i][j]));
						}
					}
				}
			}
		}
	}

	static const Double PI2 = M_PI*0.5;
	vector< NRMat<Double> > ortho_array(3, NRMat<double>(3, 3, 0.));
	vector< NRMat<Double> > dortho_array(3, NRMat<double>(3, 3, 0.));

	setOrthoMatrix(param[8], 0., 0.,  ortho_array[0]);
	setOrthoMatrix(0., param[9], 0.,  ortho_array[1]);
	setOrthoMatrix(0., 0., param[10], ortho_array[2]);

	setOrthoMatrix(param[8]+PI2, 0., 0.,  dortho_array[0]);
	setOrthoMatrix(0., param[9]+PI2, 0.,  dortho_array[1]);
	setOrthoMatrix(0., 0., param[10]+PI2, dortho_array[2]);
	dortho_array[0][2][2]= 0.;
	dortho_array[1][0][0]= 0.;
	dortho_array[2][2][2]= 0.;

	vector< NRMat<Double> > dc_dangle(3, NRMat<double>(3, 3));
	dc_dangle[0] = mprod(A_matrix, mprod(dortho_array[0], mprod( ortho_array[1],  ortho_array[2])));
	dc_dangle[1] = mprod(A_matrix, mprod( ortho_array[0], mprod(dortho_array[1],  ortho_array[2])));
	dc_dangle[2] = mprod(A_matrix, mprod( ortho_array[0], mprod( ortho_array[1], dortho_array[2])));

	for(size_t i2=0, k2=8; k2<11; i2++, k2++)
	{
		for(int i=0, k=MainParaNum; i<3; i++)
		{
			for(int j=0; j<3; j++, k++)
			{
				m_cmat[k].constr.push_back(index_set<double>(m_indep_index[k2], dc_dangle[i2][i][j]));
			}
		}
	}

	return ans;
}

bool EBSDIndexingModel::putFittingCoef_X_Z(const VecDat3<Int4>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	// Set the value X/Z = x*(x*Delta_x + y*Delta_y - (1-Delta_z)*z/(x^2+y^2)) in lclparam, where
    //
	// (x, y, z) = (h k l)*m_reciprocal_lattice_basis.
	//
	Double xValue = 0.;
	Double yValue = 0.;
	Double zValue = 0.;
	for (int i=0; i<3; i++)
	{
		xValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][0].value;
		yValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][1].value;
		zValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][2].value;
	}
	if( xValue*xValue + yValue*yValue <= 0.0 ) return false;

	const double DeltaZp1 = 1. - m_param[2].value;
	const double INV_XX_YY = 1.0/(xValue*xValue + yValue*yValue);
	const double lp0 = (m_param[0].value*xValue + m_param[1].value*yValue - DeltaZp1*zValue)*INV_XX_YY;
	lclparam = xValue*lp0;

	if( dlclparam_dother != NULL )
	{
		Int4 j = 0;
		{
			// Set d(X/Z) / d Delta_x = x^2/(x^2+y^2).
			if( m_cmat[0].ID == _IDVary ) dlclparam_dother[j++] = xValue*xValue*INV_XX_YY;
		}

		{
			// Set d(X/Z) / d Delta_y = x*y/(x^2+y^2).
			if( m_cmat[1].ID == _IDVary ) dlclparam_dother[j++] = xValue*yValue*INV_XX_YY;
		}

		{
			// Set d(X/Z) / d Delta_z = x*z/(x^2+y^2).
			if( m_cmat[2].ID == _IDVary ) dlclparam_dother[j++] = xValue*zValue*INV_XX_YY;
		}

		for(; j<m_num_indep; j++) dlclparam_dother[j] = 0.0;

		{
			// = d(X/Z)/dx = lp0 + x*( Delta_x - 2*x*(x*Delta_x + y*Delta_y - z*(1-Delta_z)/(x^2+y^2) )/(x^2+y^2)
			const Double dXZ_dx = lp0 + xValue*(m_param[0].value - 2.0*xValue*lp0)*INV_XX_YY;
			// = d(X/Z)/dy = x*( Delta_y - 2*y*(x*Delta_x + y*Delta_y - z*(1-Delta_z)/(x^2+y^2) )/(x^2+y^2)
			const Double dXZ_dy = xValue*(m_param[1].value - 2.0*yValue*lp0)*INV_XX_YY;
			// = d(X/Z)/dz = -x*(1-Delta_z)/(x^2+y^2)
			const Double dXZ_dz = -xValue * DeltaZp1 * INV_XX_YY;

			// Set d(X/Z)/d m_reciprocal_lattice_basis[i][j]
			//  = d(X/Z)/dx * dx/dm_r[i][j] + d(X/Z)/dy * dy/dm_r[i][j] + d(X/Z)/dz * dz/dm_r[i][j].
			NRMat<double> dXZ_dc(3, 3);
			for (int i=0; i<3; i++){
				dXZ_dc[i][0] = dXZ_dx * hkl[i];	// dx/dm_r[i][0] = hkl[i]. (x = \sum_i hkl[i]*dm_r[i][0])
				dXZ_dc[i][1] = dXZ_dy * hkl[i];	// dy/dm_r[i][1] = hkl[i].
				dXZ_dc[i][2] = dXZ_dz * hkl[i];	// dz/dm_r[i][2] = hkl[i].
			}
			for(int i=0, k=MainParaNum; i<3; i++)
			{
				for(int j=0; j<3; j++, k++)
				{
					for(Vec_DP_save::const_iterator it=m_cmat[k].constr.begin(); it<m_cmat[k].constr.end(); it++)
					{
						assert( (Int4)it->index < m_num_indep );
						dlclparam_dother[it->index] += it->element * dXZ_dc[i][j];
					}
				}
			}
		}
	}

	return true;
}

bool EBSDIndexingModel::putFittingCoef_Y_Z(const VecDat3<Int4>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	// Set the value Y/Z = y*(x*Delta_x + y*Delta_y - (1-Delta_z)*z/(x^2+y^2) in lclparam, where
    //
	// (x, y, z) = (h k l)*m_reciprocal_lattice_basis.
	//
	Double xValue = 0.;
	Double yValue = 0.;
	Double zValue = 0.;
	for (int i=0; i<3; i++)
	{
		xValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][0].value;
		yValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][1].value;
		zValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][2].value;
	}
	if( xValue*xValue + yValue*yValue <= 0.0) return false;

	const Double DeltaZp1 = 1. - m_param[2].value;
	const Double INV_XX_YY = 1.0/(xValue*xValue + yValue*yValue);
	const double lp0 = (m_param[0].value*xValue + m_param[1].value*yValue - DeltaZp1*zValue)*INV_XX_YY;
	lclparam = yValue*lp0;

	if( dlclparam_dother != NULL )
	{
		Int4 j = 0;
		{
			// Set d(Y/Z)/d Delta_x = x*y/(x^2+y^2).
			if( m_cmat[0].ID == _IDVary ) dlclparam_dother[j++] = xValue*yValue*INV_XX_YY;
		}
		{
			// Set d(Y/Z)/d Delta_y = y^2/(x^2+y^2).
			if( m_cmat[1].ID == _IDVary ) dlclparam_dother[j++] = yValue*yValue*INV_XX_YY;
		}
		{
			// Set d(Y/Z)/d Delta_z = y*z/(x^2+y^2).
			if( m_cmat[2].ID == _IDVary ) dlclparam_dother[j++] = yValue*zValue*INV_XX_YY;
		}

		for(; j<m_num_indep; j++) dlclparam_dother[j] = 0.0;

		{
			// = d(Y/Z)/dx = y*( Delta_x - 2*x*(x*Delta_x + y*Delta_y - z*(1-Delta_z)/(x^2+y^2) )/(x^2+y^2)
			const Double dYZ_dx = yValue*(m_param[0].value - 2.0*xValue*lp0)*INV_XX_YY;
			// = d(Y/Z)/dy = lp0 + y*( Delta_y - 2*y*(x*Delta_x + y*Delta_y - z*(1-Delta_z)/(x^2+y^2) )/(x^2+y^2)
			const Double dYZ_dy = lp0 + yValue*(m_param[1].value - 2.0*yValue*lp0)*INV_XX_YY;
 			// = d(Y/Z)/dz = -y*(1-Delta_z)/(x^2+y^2)
 			const Double dYZ_dz = -yValue * DeltaZp1 * INV_XX_YY;

			// Set d(Y/Z) / d m_reciprocal_lattice_basis[i][j]
			//  = dYZ_dx * dx/dm_r[i][j] + dYZ_dy * dy/dm_r[i][j] + dYZ_dz * dz/dm_r[i][j].
			NRMat<double> dYZ_dc(3, 3);
			for (int i=0; i<3; i++){
				dYZ_dc[i][0] = dYZ_dx * hkl[i];	// dx / d m_r[i][0] = hkl[i].
				dYZ_dc[i][1] = dYZ_dy * hkl[i];	// dy / d m_r[i][1] = hkl[i].
				dYZ_dc[i][2] = dYZ_dz * hkl[i];	// dz / d m_r[i][2] = hkl[i].
			}
			for(int i=0, k=MainParaNum; i<3; i++)
			{
				for(int j=0; j<3; j++, k++)
				{
					for(Vec_DP_save::const_iterator it=m_cmat[k].constr.begin(); it<m_cmat[k].constr.end(); it++)
					{
						assert( (Int4)it->index < m_num_indep );
						dlclparam_dother[it->index] += it->element * dYZ_dc[i][j];
					}
				}
			}
		}
	}

	return true;
}


bool EBSDIndexingModel::putFittingCoef_BandWidth(const VecDat3<Int4>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	// Set the following in lclparam:
	// BandWidth = (tan(sigma + theta) - tan(sigma - theta))*(1 - Delta_z),
	//  sigma = arctan(z/sqrt(x^2 + y^2)),
	//  theta = arcsin(wavelength*ScaleFactor*sqrt(x^2+y^2+z^2)/2), where
    //
	// (x, y, z) = (h k l)*m_reciprocal_lattice_basis.
	//
	Double xValue = 0.;
	Double yValue = 0.;
	Double zValue = 0.;
	for (int i=0; i<3; i++)
	{
		xValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][0].value;
		yValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][1].value;
		zValue  +=  hkl[i] * m_reciprocal_lattice_basis[i][2].value;
	}
	if( xValue*xValue + yValue*yValue <= 0.0) return false;
	if( zValue <= 0.0 ) return false;

	const Double d_star_cal = sqrt(xValue*xValue + yValue*yValue + zValue*zValue);
	if( fabs(m_wavelength*m_param[11].value*d_star_cal) > 2.0 ) return false;
	const Double theta = asin(0.5*m_wavelength*m_param[11].value*d_star_cal); assert( theta >= 0.0 );
	const Double INV_XX_YY = 1.0/(xValue*xValue + yValue*yValue);
	const Double sigma = atan(zValue*sqrt(INV_XX_YY)); assert( sigma > 0.0 );

	const Double tan_sigmab = tan(sigma-theta);
	const Double tan_sigmae = tan(sigma+theta);
	const Double inv_cos_theta = 1.0/cos(theta);

	// BandWidth = (tan(sigma + theta) - tan(sigma - theta))*(1 - Delta_z),
	const Double BandWidth0 = tan_sigmae - tan_sigmab;
	const Double DeltaZp1 = 1. - m_param[2].value;
	lclparam = BandWidth0*DeltaZp1;
	const Double dBW_dtheta = (2.0 + tan_sigmae*tan_sigmae + tan_sigmab*tan_sigmab)*DeltaZp1; // = 1.0/cos^2(sigmae) + 1.0/cos^2(sigmab).

	if( dlclparam_dother != NULL )
	{
		for(Int4 j=0; j<m_num_indep; j++) dlclparam_dother[j] = 0.0;

		{
			// Set d BandWidth / d Delta_z = tan_sigmae - tan_sigmab.
			if( m_cmat[2].ID == _IDVary ) dlclparam_dother[m_indep_index[2]] = -BandWidth0;
		}
		{
			// cos(theta)*dtheta = 0.5*m_wavelength*d_star_cal*dscale.
			const Double dtheta_dscale = 0.5*inv_cos_theta*m_wavelength*d_star_cal;
			if( m_cmat[11].ID == _IDVary ) dlclparam_dother[m_indep_index[11]] = dBW_dtheta*dtheta_dscale;
		}
/*
		{
			const double dtheta_dqstar_2 = 0.5*inv_cos_theta*m_wavelength*m_param[11].value/d_star_cal;

			// Herein, dBandWidth/dsigma = 0 (the influence of sigma is ignored).
			const Double dBW_dx = dBW_dtheta*dtheta_dqstar_2*xValue;
			const Double dBW_dy = dBW_dtheta*dtheta_dqstar_2*yValue;
			const Double dBW_dz = dBW_dtheta*dtheta_dqstar_2*zValue;

			// Set dBW / d m_reciprocal_lattice_basis[i][j]
			//  = dBW_dx * dx/dm_r[i][j] + dBW_dy * dy/dm_r[i][j] + dBW_dz * dz/dm_r[i][j].
			NRMat<double> dBW_dc(3, 3);
			for (int i=0; i<3; i++){
				dBW_dc[i][0] = dBW_dx * hkl[i];	// dx / d m_r[i][0] = hkl[i].
				dBW_dc[i][1] = dBW_dy * hkl[i];	// dz / d m_r[i][2] = hkl[i].
				dBW_dc[i][2] = dBW_dz * hkl[i];	// dz / d m_r[i][2] = hkl[i].
			}
			for(int i=0, k=MainParaNum; i<3; i++)
			{
				for(int j=0; j<3; j++, k++)
				{
					for(Vec_DP_save::const_iterator it=m_cmat[k].constr.begin(); it<m_cmat[k].constr.end(); it++)
					{
						assert( (Int4)it->index < m_num_indep );
						dlclparam_dother[it->index] += it->element * dBW_dc[i][j];
					}
				}
			}
		}
*/
	}

	return true;
}

bool EBSDIndexingModel::putFittingCoef(const pair< VecDat3<Int4>, char>& hkl, Double& lclparam, Double* dlclparam_dother, ostringstream& outPutString) const
{
	if( hkl.second == 0 )
	{
		return putFittingCoef_X_Z(hkl.first, lclparam, dlclparam_dother, outPutString);
	}
	else if( hkl.second == 1 )
	{
		return putFittingCoef_Y_Z(hkl.first, lclparam, dlclparam_dother, outPutString);
	}
	else
	{
		return putFittingCoef_BandWidth(hkl.first, lclparam, dlclparam_dother, outPutString);
	}

	return true;
}


// If wt != 2, lclsterr is used as weights in least square method. 
// If wt = 2, lclsterr is not used,  all the weights equal 1. 
pair<bool, ZErrorMessage> EBSDIndexingModel::setFittedParam(const vector< pair<VecDat3<Int4>, char> >& hkl,
		const vector<Double>& lclparam, const vector<Double>& lclvar,
		const vector<bool>& nxfit,
		const bool& output_view_flag, const Double& judge_conv,  const Int4& max_itnum,
		const vector<Double>& glbinit, Double& chisq_init, Double& chisq_all)
{
	Int4 itnum;
	
	LemarqMethod marq;
	marq.setParam(output_view_flag, judge_conv, max_itnum, 0.001);
	
	MarquardtFmodelBase< pair<VecDat3<Int4>, char> > LMFunc(*this);
	ZErrorMessage zerr = LMFunc.addHistogram(hkl, lclparam, lclvar, nxfit);
	if( zerr.putErrorType() != ZErrorNoError ) return pair<bool, ZErrorMessage>(false, zerr);
	
	Vec_DP chisq;
	Int4 reason_terminate_LM;
	return marq.execute< pair<VecDat3<Int4>, char> >(glbinit, LMFunc, itnum, chisq_init, chisq_all, chisq, reason_terminate_LM);
}
