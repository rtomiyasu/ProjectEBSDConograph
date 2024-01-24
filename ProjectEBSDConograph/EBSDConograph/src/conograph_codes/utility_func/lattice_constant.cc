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
#include "lattice_constant.hh"
#include "zmath.hh"


static const Double DegRad = 180.0 / M_PI;
static const Double RadDeg = M_PI /180.0; // = pi / 180.0;

static void calLatticeConstant(const SymMat<Double>& S, vector<Double>& lattice_consts,
VecDat3<Double>& cos_angle, VecDat3<Double>& sin_angle)
{
	assert( S.size() >= 3 );

	// a*, b*, c*
	const VecDat3<Double> co_length( sqrt( S(0,0) ), sqrt( S(1,1) ), sqrt( S(2,2) ) );
	// alpha*, beta*, gamma*
	const VecDat3<Double> cos_co_angle( S(1,2)/(co_length[1]*co_length[2]), S(0,2)/(co_length[0]*co_length[2]), S(0,1)/(co_length[0]*co_length[1]) );
	const VecDat3<Double> sin_co_angle( 	sqrt(1.0 - cos_co_angle[0]*cos_co_angle[0]),
											sqrt(1.0 - cos_co_angle[1]*cos_co_angle[1]),
											sqrt(1.0 - cos_co_angle[2]*cos_co_angle[2])	);
	const Double inv_vol = co_length[0]*co_length[1]*co_length[2]
							* sqrt( 1.0-cos_co_angle[0]*cos_co_angle[0]-cos_co_angle[1]*cos_co_angle[1]-cos_co_angle[2]*cos_co_angle[2]
																+ 2.0*cos_co_angle[0]*cos_co_angle[1]*cos_co_angle[2] );
	const Double vol = 1.0/inv_vol;

	// alpha, beta, gamma
	cos_angle[0] = (cos_co_angle[1]*cos_co_angle[2]-cos_co_angle[0])/(sin_co_angle[1]*sin_co_angle[2]);
	cos_angle[1] = (cos_co_angle[0]*cos_co_angle[2]-cos_co_angle[1])/(sin_co_angle[0]*sin_co_angle[2]);
	cos_angle[2] = (cos_co_angle[0]*cos_co_angle[1]-cos_co_angle[2])/(sin_co_angle[0]*sin_co_angle[1]);
	sin_angle[0] = sqrt(1.0 - cos_angle[0]*cos_angle[0]);
	sin_angle[1] = sqrt(1.0 - cos_angle[1]*cos_angle[1]);
	sin_angle[2] = sqrt(1.0 - cos_angle[2]*cos_angle[2]);

	lattice_consts.resize(6);
	lattice_consts[0] = co_length[1]*co_length[2]*sin_co_angle[0]*vol;
	lattice_consts[1] = co_length[0]*co_length[2]*sin_co_angle[1]*vol;
	lattice_consts[2] = co_length[0]*co_length[1]*sin_co_angle[2]*vol;
	lattice_consts[3] = atan2(sin_angle[0], cos_angle[0]) * DegRad;
	lattice_consts[4] = atan2(sin_angle[1], cos_angle[1]) * DegRad;
	lattice_consts[5] = atan2(sin_angle[2], cos_angle[2]) * DegRad;
}

void calLatticeConstant(const SymMat<Double>& S, vector<Double>& lattice_consts)
{
	// alpha, beta, gamma
	VecDat3<Double> cos_angle, sin_angle;
	calLatticeConstant(S, lattice_consts, cos_angle, sin_angle);
}

// S_covar is the covariant matrix on S(A*,B*,C*,D*,E*,F*).
void changeCovariantMatrixStoLatticeConstant(const SymMat<Double>& S_covar, const vector<Double>& lattice_consts,
const VecDat3<Double>& cos_angle, const VecDat3<Double>& sin_angle, SymMat<Double>& ans)
{
	static const Int4 pn_lat_const = 6;
	assert( S_covar.size() == pn_lat_const );

	const Double aa=0.5*lattice_consts[0]*lattice_consts[0];
	const Double bb=0.5*lattice_consts[1]*lattice_consts[1];
	const Double cc=0.5*lattice_consts[2]*lattice_consts[2];
	const Double ab=lattice_consts[0]*lattice_consts[1];
	const Double ac=lattice_consts[0]*lattice_consts[2];
	const Double bc=lattice_consts[1]*lattice_consts[2];

	const Double cos2alpha=cos_angle[0]*cos_angle[0];
	const Double cos2beta=cos_angle[1]*cos_angle[1];
	const Double cos2gamma=cos_angle[2]*cos_angle[2];
	const Double cos_alpha_beta=cos_angle[0]*cos_angle[1];
	const Double cos_alpha_gamma=cos_angle[0]*cos_angle[2];
	const Double cos_beta_gamma=cos_angle[1]*cos_angle[2];

	// Jacobian matrix from S to lattice constants.
	NRMat<Double> Jac(pn_lat_const, pn_lat_const);

	Jac[0][0] = lattice_consts[0] * aa;
	Jac[0][1] = lattice_consts[0] * bb * cos2gamma;
	Jac[0][2] = lattice_consts[0] * cc * cos2beta;
	Jac[0][3] = lattice_consts[0] * bc * cos_beta_gamma;
	Jac[0][4] = lattice_consts[0] * ac * cos_angle[1];
	Jac[0][5] = lattice_consts[0] * ab * cos_angle[2];

	Jac[1][0] = lattice_consts[1] * aa * cos2gamma;
	Jac[1][1] = lattice_consts[1] * bb;
	Jac[1][2] = lattice_consts[1] * cc * cos2alpha;
	Jac[1][3] = lattice_consts[1] * bc * cos_angle[0];
	Jac[1][4] = lattice_consts[1] * ac * cos_alpha_gamma;
	Jac[1][5] = lattice_consts[1] * ab * cos_angle[2];

	Jac[2][0] = lattice_consts[2] * aa * cos2beta;
	Jac[2][1] = lattice_consts[2] * bb * cos2alpha;
	Jac[2][2] = lattice_consts[2] * cc;
	Jac[2][3] = lattice_consts[2] * bc * cos_angle[0];
	Jac[2][4] = lattice_consts[2] * ac * cos_angle[1];
	Jac[2][5] = lattice_consts[2] * ab * cos_alpha_beta;

	Jac[3][0] = -DegRad / sin_angle[0] * aa * ( 2.0 * cos_beta_gamma - cos_angle[0] * (cos2beta + cos2gamma) );
	Jac[3][1] = -DegRad * sin_angle[0] * bb * cos_angle[0];
	Jac[3][2] = -DegRad * sin_angle[0] * cc * cos_angle[0];
	Jac[3][3] = -DegRad * sin_angle[0] * bc;
	Jac[3][4] = -DegRad * sin_angle[0] * ac * cos_angle[2];
	Jac[3][5] = -DegRad * sin_angle[0] * ab * cos_angle[1];

	Jac[4][0] = -DegRad * sin_angle[1] * aa * cos_angle[1];
	Jac[4][1] = -DegRad / sin_angle[1] * bb * ( 2.0 * cos_alpha_gamma - cos_angle[1] * (cos2alpha + cos2gamma) );
	Jac[4][2] = -DegRad * sin_angle[1] * cc * cos_angle[1];
	Jac[4][3] = -DegRad * sin_angle[1] * bc * cos_angle[2];
	Jac[4][4] = -DegRad * sin_angle[1] * ac;
	Jac[4][5] = -DegRad * sin_angle[1] * ab * cos_angle[0];

	Jac[5][0] = -DegRad * sin_angle[2] * aa * cos_angle[2];
	Jac[5][1] = -DegRad * sin_angle[2] * bb * cos_angle[2];
	Jac[5][2] = -DegRad / sin_angle[2] * cc * ( 2.0 * cos_alpha_beta - cos_angle[2] * (cos2alpha + cos2beta) );
	Jac[5][3] = -DegRad * sin_angle[2] * bc * cos_angle[1];
	Jac[5][4] = -DegRad * sin_angle[2] * ac * cos_angle[0];
	Jac[5][5] = -DegRad * sin_angle[2] * ab;

	ans = SymMat<Double>(pn_lat_const, 0.0);

	// ans = Jac * S_covar * transpose(Jac).
	for(Int4 j=0; j<pn_lat_const; j++)
	{
		for(Int4 k=0; k<=j; k++)
		{
			for(Int4 l=0; l<pn_lat_const; l++)
			{
				ans(j,k) += Jac[j][l] * Jac[k][l] * S_covar(l,l);
				for(Int4 m=0; m<l; m++) ans(j,k) += (Jac[j][l] * Jac[k][m] + Jac[k][l] * Jac[j][m]) * S_covar(l,m);
			}
		}
	}
}


void calLatticeConstant(const SymMat<Double>& S, const SymMat<Double>& S_covar, vector<Double>& lattice_consts,
SymMat<Double>& LatConst_covar)
{
	// alpha, beta, gamma
	VecDat3<Double> cos_angle, sin_angle;
	calLatticeConstant(S, lattice_consts, cos_angle, sin_angle);
	changeCovariantMatrixStoLatticeConstant(S_covar, lattice_consts, cos_angle, sin_angle, LatConst_covar);
}
