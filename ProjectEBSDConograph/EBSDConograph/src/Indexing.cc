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

#include <iomanip>
#include <fstream>
#include <cmath>
#include <set>
#include <map>
#include "Indexing.hh"
#include "marquardt_codes/model_function/EBSDIndexingModel.hh"
#include "marquardt_codes/levenberg_marquardt/SVdcmp.hh"
#include "conograph_codes/bravais_lattice_determination/BravaisLattice.hh"
#include "conograph_codes/bravais_lattice_determination/LatticeFigureOfMerit.hh"
#include "conograph_codes/bravais_lattice_determination/LatticeFigureOfMeritToCheckSymmetry.hh"
#include "conograph_codes/HC_algorithm/Determination_of_shortest_vectors.hh"
#include "conograph_codes/utility_data_structure/SymMat.hh"
#include "conograph_codes/utility_lattice_reduction/put_Selling_reduced_lattice.hh"
#include "conograph_codes/utility_func/lattice_constant.hh"
#include "conograph_codes/utility_func/zmath.hh"
#include "zlog/zlog.hh"


namespace {
	template<class T>
	inline T square(const T& arg)
	{
		return arg * arg;
	}

	inline NRMat<double> chToDouble(const NRMat<FittingParameter>& basis)
	{
		const int nr = basis.nrows();
		const int nc = basis.ncols();
		NRMat<double> ans(nr,nc);
		for(Int4 i2=0; i2<nr; i2++)
		{
			for(Int4 j2=0; j2<nc; j2++)
			{
				ans[i2][j2] = basis[i2][j2].value;
			}
		}
		return ans;
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

	NRMat<double> mprod(const NRMat<int>& lhs, const NRMat<double>& rhs)
	{
		const int irow = lhs.nrows();
		const int icol = rhs.ncols();
		const int n = lhs.ncols();
		assert( lhs.ncols()==rhs.nrows() );


		NRMat<double> ans( irow, icol, 0 );
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
	NRVec<T> right_act(const NRVec<int>&v, const NRMat<T> &m) {
		const int nc = m.ncols();
		const int nr = m.nrows();
		assert( nr == v.size() );
		NRVec<T> r(nc, 0);
		for(int i=0; i<nc; i++) {
			for(int j=0; j<nr; j++) {
				r[i] += v[j] * m[j][i];
			}
		}
		return r;
	}
}

namespace {
	inline bool equiv(const Indexing::OrientationVector& arg1, const Indexing::OrientationVector& arg2)
	{
		if(square(arg1.X0 - arg2.X0) <= square(arg1.err_X0) + square(arg2.err_X0)
				&& square(arg1.Y0 - arg2.Y0) <= square(arg1.err_Y0) + square(arg2.err_Y0)
				&& square(arg1.Z0 - arg2.Z0) <= square(arg1.err_Z0) + square(arg2.err_Z0) )
		{
			return true;
		}
		return false;
	}

	inline bool equiv(Indexing::OrientationVector arg1, const double& length1,
					  Indexing::OrientationVector arg2, const double& length2)
	{
		arg1.X0 *= length1;
		arg1.Y0 *= length1;
		arg1.Z0 *= length1;
		arg1.err_X0 *= length1;
		arg1.err_Y0 *= length1;
		arg1.err_Z0 *= length1;

		arg2.X0 *= length2;
		arg2.Y0 *= length2;
		arg2.Z0 *= length2;
		arg2.err_X0 *= length2;
		arg2.err_Y0 *= length2;
		arg2.err_Z0 *= length2;
		return equiv(arg1, arg2);
	}

	inline pair<Indexing::OrientationVector, double> put_vector_sum(const Indexing::OrientationVector& vec1, const double& len_vec1, const Indexing::OrientationVector& vec2, const double& len_vec2)
	{
		assert( fabs(vec1.X0*vec1.X0 + vec1.Y0*vec1.Y0 + vec1.Z0*vec1.Z0 - 1.0) <= 1.0e-14 );
		assert( fabs(vec2.X0*vec2.X0 + vec2.Y0*vec2.Y0 + vec2.Z0*vec2.Z0 - 1.0) <= 1.0e-14 );
		Indexing::OrientationVector vec_sum;
		vec_sum.X0 = vec1.X0*len_vec1 + vec2.X0*len_vec2;
		vec_sum.Y0 = vec1.Y0*len_vec1 + vec2.Y0*len_vec2;
		vec_sum.Z0 = vec1.Z0*len_vec1 + vec2.Z0*len_vec2;
		vec_sum.err_X0 = sqrt(square(vec1.err_X0*len_vec1) + square(vec2.err_X0*len_vec2));
		vec_sum.err_Y0 = sqrt(square(vec1.err_Y0*len_vec1) + square(vec2.err_Y0*len_vec2));
		vec_sum.err_Z0 = sqrt(square(vec1.err_Z0*len_vec1) + square(vec2.err_Z0*len_vec2));

		const double len = sqrt(vec_sum.X0*vec_sum.X0 + vec_sum.Y0*vec_sum.Y0 + vec_sum.Z0*vec_sum.Z0);

		vec_sum.X0 /= len;
		vec_sum.Y0 /= len;
		vec_sum.Z0 /= len;
		vec_sum.err_X0 /= len;
		vec_sum.err_Y0 /= len;
		vec_sum.err_Z0 /= len;
		return make_pair(vec_sum, len);
	}
}

namespace {
	void find_lambda(const vector<Indexing::DataSet2>& data, const size_t& i, const size_t& j,
			//sconst bool& isLambdaPositive,
			map<size_t, pair<double, double> >& tray_k_lambda)
	{
		tray_k_lambda.clear();

		NRMat<double> alpha(3,3);
		NRVec<double> d(3);

		const Indexing::OrientationVector veci = data[i].putOrientationVector();
		const Indexing::OrientationVector vecj = data[j].putOrientationVector();

		for(size_t k=0; k<data.size(); k++){
			if(i == k){
				continue;
			}
			if(j == k){
				continue;
			}

			const Indexing::OrientationVector veck = data[k].putOrientationVector();

/*
			//                      (data[i].X0 data[i].Y0 data[i].Z0)
			// alpha = V * V^T, V = (data[j].X0 data[j].Y0 data[j].Z0)
			//                      (data[k].X0 data[k].Y0 data[k].Z0)
			alpha[0][0] = 1.;
			alpha[1][1] = 1.;
			alpha[2][2] = 1.;
			alpha[0][1] = data[i].X0*data[j].X0 + data[i].Y0*data[j].Y0 + data[i].Z0*data[j].Z0;
			alpha[0][2] = data[i].X0*data[k].X0 + data[i].Y0*data[k].Y0 + data[i].Z0*data[k].Z0;
			alpha[1][2] = data[j].X0*data[k].X0 + data[j].Y0*data[k].Y0 + data[j].Z0*data[k].Z0;
			alpha[1][0] = alpha[0][1];
			alpha[2][0] = alpha[0][2];
			alpha[2][1] = alpha[1][2];

			SVdcmp(alpha, d, 0, 3);

			// | (alpha[0][col_index] alpha[1][col_index] alpha[2][col_index]) * V | is very small.
			int col_index = 0;
			if( d[1] < d[col_index] ) col_index = 1;
			if( d[2] < d[col_index] ) col_index = 2;

			const double lambda1 = -alpha[0][col_index] / alpha[2][col_index];
			const double lambda2 = -alpha[1][col_index] / alpha[2][col_index];
*/
			// (xi xj)*(lambda1)=(xk)
			// (yi yj) (lambda2) (yk)
			const double inv_det = 1.0 / (veci.X0*vecj.Y0 - veci.Y0*vecj.X0);
			const double lambda1 = (veck.X0*vecj.Y0 - veck.Y0*vecj.X0) * inv_det;
			const double lambda2 = (veck.Y0*veci.X0 - veck.X0*veci.Y0) * inv_det;
/*			if( isLambdaPositive )
			{
				if( lambda1 <= 0 )
				{
					continue;
				}
				if( lambda2 <= 0 )
				{
					continue;
				}
			}
*/

			const pair<Indexing::OrientationVector, double> vec_sum = put_vector_sum(data[i].putOrientationVector(), lambda1, data[j].putOrientationVector(), lambda2);

			if(equiv(vec_sum.first, data[k].putOrientationVector()))
			{
/*
if( isLambdaPositive )
{
cout << "The smallest eigenvalue of alpha = " << d[row_index] << endl;
cout << "data[" << i << "] = (" << data[i].x << "," << data[i].y << "," << data[i].z << ")\n";
cout << "data[" << j << "] = (" << data[j].x << "," << data[j].y << "," << data[j].z << ")\n";
cout << "data[" << k << "] = (" << data[k].x << "," << data[k].y << "," << data[k].z << ")\n";
cout << lambda1 << "*data[" << i << "] + " << lambda2 << "*data["<< j << "] = (" << vec_sum.x << "," << vec_sum.y << "," << vec_sum.z << ")\n\n";
}
*/
				tray_k_lambda[k] = make_pair(lambda1, lambda2);
			}
		}
	}

	int find(const vector<Indexing::DataSet2>& data, const map<size_t, pair<double, double> >& tray_k_lambda, const Indexing::OrientationVector& arg)
	{
		assert( fabs(arg.X0*arg.X0 + arg.Y0*arg.Y0 + arg.Z0*arg.Z0 - 1.0) <= 1.0e-14 );
		if (arg.Z0 <= 0) return -1;
		for(map<size_t, pair<double, double> >::const_iterator it = tray_k_lambda.begin(); it !=tray_k_lambda.end(); it++)
		{
			if(equiv(data[it->first].putOrientationVector(), arg)) return it->first;
		}
		return -1;
	}

	void make_pair_ij(const int& i, const Indexing::OrientationVector& data_i, const double& kappa1,
						const int& j, const Indexing::OrientationVector& data_j, const double& kappa2, const FittingParameter ScaleFactor,
						vector<Indexing::PairOfReciprocalVectors>& pairs_tray)  //
	{
//		assert( i >= 0 );
//		assert( j >= 0 );
		Indexing::PairOfReciprocalVectors pair_ij;
		pair_ij.vec1.index_in_data = i;
		pair_ij.vec1.orientation = data_i;
		pair_ij.vec1.length = kappa1;
		pair_ij.vec2.index_in_data = j;
		pair_ij.vec2.orientation = data_j;
		pair_ij.vec2.length = kappa2;
		pair_ij.ScaleFactor = ScaleFactor;

		pairs_tray.push_back(pair_ij);
		swap(pair_ij.vec1, pair_ij.vec2);
		pairs_tray.push_back(pair_ij);
	}

	void find_pair_of_reciprocal_vectors(const vector<Indexing::DataSet2>& data,
			const double& wavelength,
			const bool& get_unitcell_scale,
			// const bool& exhaustive_search,
			// const double& toleranceOnBandWidth,
			vector<Indexing::PairOfReciprocalVectors>& pairs_tray)
	{
		pairs_tray.clear();

		// 異なるm_data[i]=(xi, yi, zi), m_data[j]=(xj, yj, zj)に対して以下を行う。
		vector< FittingParameter> ScaleFactor(2, 0.0);
		for(size_t i=0; i<data.size(); i++)
		{
			for(size_t j=i+1; j<data.size(); j++)
			{
				// k != i,j, data[k]=(xk, yk, zk)で、(xk, yk, zk) = lambda1*(xi, yi, zi) + lambda2*(xj, yj, zj) となる　lambda1, lambda2。
				map<size_t, pair<double, double> > tray_k_lambda;
				find_lambda(data, i, j, tray_k_lambda);
/*
if( !tray_k_lambda.empty() )
{
	cout << "i=" << i << " j=" << j << " k=";
	for(map<size_t, pair<double, double> >::const_iterator it=tray_k_lambda.begin(); it!=tray_k_lambda.end(); it++)
	{
		const double& lambda1 = it->second.first;
		const double& lambda2 = it->second.second;
		if( lambda1 <= 0.0 || lambda2 <= 0.0 ) continue;
		cout << it->first << " ";
	}
	cout << endl;
}
*/
				for(map<size_t, pair<double, double> >::const_iterator it=tray_k_lambda.begin(); it!=tray_k_lambda.end(); it++)
				{
					const double& lambda1 = it->second.first;
					const double& lambda2 = it->second.second;
					if( lambda1 <= 0.0 || lambda2 <= 0.0 ) continue;
					// (i)
					// |a1 + a2|*(xk, yk, zk) = a1 + a2 を仮定した場合。
					// 長さの比 |a1| : |a2| : |a1 + a2| = lambda1 : lambda2 : 1.0;

					// (ii)
					// |2*a1 + a2|*(xk, yk, zk) = 2*a1 + a2 を仮定した場合。
					// 長さの比 |2*a1| : |a2| : |2*a1 + a2| = lambda1 : lambda2 : 1.0;

					// (iii)
					// |a1 + 2*a2|*(xk, yk, zk) = a1 + 2*a2 を仮定した場合。
					// 長さの比 |a1| : |2*a2| : |a1 + 2*a2| = lambda1 : lambda2 : 1.0;

//cout << lambda1 << " " << lambda2 << " " << lambda1+lambda2 << endl;
					// 正定値性を以下でチェックできる。
					assert( fabs(lambda1 - lambda2) < 1.0 );
					assert( 1.0 < lambda1 + lambda2 );

					if( get_unitcell_scale )
					{
						const double s1 = data[i].d_star_obs / lambda1;
						const double s2 = data[j].d_star_obs / lambda2;
						ScaleFactor[0] = FittingParameter(s1, data[i].err_d_star_obs / lambda1);
						ScaleFactor[1] = FittingParameter(s2, data[j].err_d_star_obs / lambda2);
					}

					const pair<Indexing::OrientationVector, double> a1pa2_ii  = put_vector_sum( data[i].putOrientationVector(), lambda1*0.5, data[j].putOrientationVector(), lambda2    );
					const pair<Indexing::OrientationVector, double> a1pa2_iii = put_vector_sum( data[i].putOrientationVector(), lambda1    , data[j].putOrientationVector(), lambda2*0.5);
					const int index_a1pa2[2] = { find(data, tray_k_lambda, a1pa2_ii.first),	// =  (a1 +2*a2)/2 in case (i), =a1 + a2 in case (ii)
										   	   	 find(data, tray_k_lambda, a1pa2_iii.first) };// = (2*a1 + a2)/2 in case (i), =a1 + a2 in case (iii)

					// In case (i), a1+a2 is observed. Save {a1, a1+a2}, {a2, a1+a2}.
					make_pair_ij(i, data[i].putOrientationVector(), lambda1, it->first, data[it->first].putOrientationVector(), 1.0, ScaleFactor[0], pairs_tray);
					make_pair_ij(j, data[j].putOrientationVector(), lambda2, it->first, data[it->first].putOrientationVector(), 1.0, ScaleFactor[1], pairs_tray);

					// In case (ii), 2*a1+a2 is observed. If a1+a2 is not observed, save {a1, a1+a2}, {a2, a1+a2}.
					if( index_a1pa2[0] < 0 )
					{
						make_pair_ij(i, data[i].putOrientationVector(), lambda1*0.5, -1, a1pa2_ii.first, a1pa2_ii.second, ScaleFactor[0], pairs_tray);
						make_pair_ij(j, data[j].putOrientationVector(), lambda2,     -1, a1pa2_ii.first, a1pa2_ii.second, ScaleFactor[1], pairs_tray);
					}

					// In case (iii), a1+2*a2 is observed. If a1+a2 is not observed, save {a1, a1+a2}, {a2, a1+a2}.
					if( index_a1pa2[1] < 0 )
					{
						make_pair_ij(i, data[i].putOrientationVector(), lambda1    , -1, a1pa2_iii.first, a1pa2_iii.second, ScaleFactor[0], pairs_tray);
						make_pair_ij(j, data[j].putOrientationVector(), lambda2*0.5, -1, a1pa2_iii.first, a1pa2_iii.second, ScaleFactor[1], pairs_tray);
					}
				}
			}
		}
	}

	void find_basis_of_reciprocal_lattice(const vector<Indexing::PairOfReciprocalVectors>& pairs_tray,
			const Indexing::conditions& cons, const vector<Indexing::DataSet2>& data, const bool& use_the_band_widths,
			vector< vector<Indexing::result_type> >& ans)
	{
		static const size_t NUM_BT = 14;
		ans.clear();
		ans.resize(NUM_BT);

		NRMat<int> trans_mat(3,3);
		// 3次元格子の基底を選出。
		FittingParameter ScaleFactor = 0.0;
		{
			NRMat<double> basis(3, 3);
			for(size_t i1=0; i1<pairs_tray.size(); i1++)
			{
				const Indexing::InfoOfReciprocalVector& v1 = pairs_tray[i1].vec1;
				const Indexing::InfoOfReciprocalVector& v2 = pairs_tray[i1].vec2;
				const double length_ratio12 = v2.length / v1.length;
				for(size_t i2=i1+1; i2<pairs_tray.size(); i2++)
				{
					// v1とpairs_tray[i2].vec1がスケールを無視したときに等しいかどうかチェック。
					if( v1.index_in_data != pairs_tray[i2].vec1.index_in_data ) continue;
					if( v1.index_in_data < 0 && !equiv(v1.orientation, pairs_tray[i2].vec1.orientation) ) continue;

					// vec2も同じ観測値から来ているときは、使わないことにする。
					const Indexing::InfoOfReciprocalVector& v3 = pairs_tray[i2].vec2;
					if( v2.index_in_data == v3.index_in_data )
					{
						if( v2.index_in_data >= 0 ) continue;
						if( equiv(v2.orientation, v3.orientation) ) continue;
					}

					// v1, v2, v3が線形独立になっているかチェック。
					vector<Indexing::DataSet2> data_i12(3);
					data_i12[0].setOrientationVector(v1.orientation);
					data_i12[1].setOrientationVector(v2.orientation);
					data_i12[2].setOrientationVector(v3.orientation);
					map<size_t, pair<double, double> > tray_k_lambda;
					find_lambda(data_i12, 0, 1, tray_k_lambda);
					if( !tray_k_lambda.empty() ) continue;

					const double length_ratio13 = v3.length / pairs_tray[i2].vec1.length;

					if( !cons.search_level )
					{
						// Check if v1 + v2 + v3 is observed.
						Indexing::OrientationVector vec_sum;
						vec_sum.X0 = v1.orientation.X0 + v2.orientation.X0*length_ratio12 + v3.orientation.X0*length_ratio13;
						vec_sum.Y0 = v1.orientation.Y0 + v2.orientation.Y0*length_ratio12 + v3.orientation.Y0*length_ratio13;
						vec_sum.Z0 = v1.orientation.Z0 + v2.orientation.Z0*length_ratio12 + v3.orientation.Z0*length_ratio13;
						vec_sum.err_X0 = sqrt(square(v1.orientation.err_X0) + square(v2.orientation.err_X0*length_ratio12) + square(v3.orientation.err_X0*length_ratio13));
						vec_sum.err_Y0 = sqrt(square(v1.orientation.err_Y0) + square(v2.orientation.err_Y0*length_ratio12) + square(v3.orientation.err_Y0*length_ratio13));
						vec_sum.err_Z0 = sqrt(square(v1.orientation.err_Z0) + square(v2.orientation.err_Z0*length_ratio12) + square(v3.orientation.err_Z0*length_ratio13));

						const double len = sqrt(vec_sum.X0*vec_sum.X0 + vec_sum.Y0*vec_sum.Y0 + vec_sum.Z0*vec_sum.Z0);
						vec_sum.X0 /= len;
						vec_sum.Y0 /= len;
						vec_sum.Z0 /= len;
						vec_sum.err_X0 /= len;
						vec_sum.err_Y0 /= len;
						vec_sum.err_Z0 /= len;

						vector<Indexing::DataSet2>::const_iterator it;
						for(it = data.begin(); it !=data.end(); it++)
						{
							if(equiv(it->putOrientationVector(), vec_sum)) break;
						}
						if( it >= data.end() ) continue;
					}

					if( use_the_band_widths )
					{
						const FittingParameter s1 = pairs_tray[i1].ScaleFactor * pairs_tray[i1].vec1.length;
						const FittingParameter s2 = pairs_tray[i2].ScaleFactor * pairs_tray[i2].vec1.length;

						if( cons.search_level || ( s1.value <= s2.value*cons.toleranceOnBandWidth && s2.value <= s1.value*cons.toleranceOnBandWidth ) )
						{
							// ( s1.value/s1.error ) = ( 1/s1.error ) * s
							// ( s2.value/s2.error )   ( 1/s2.error )
							const double s_var = 1.0 / ( 1.0/(s1.error*s1.error) + 1.0/(s2.error*s2.error) );
							const double s = s_var * ( s1.value / (s1.error*s1.error) + s2.value / (s2.error*s2.error) );
							ScaleFactor = FittingParameter(s, sqrt(s_var));
						}
						else continue;
					}

					basis[0][0] = v1.orientation.X0;
					basis[0][1] = v1.orientation.Y0;
					basis[0][2] = v1.orientation.Z0;

					basis[1][0] = v2.orientation.X0*length_ratio12;
					basis[1][1] = v2.orientation.Y0*length_ratio12;
					basis[1][2] = v2.orientation.Z0*length_ratio12;

					basis[2][0] = v3.orientation.X0*length_ratio13;
					basis[2][1] = v3.orientation.Y0*length_ratio13;
					basis[2][2] = v3.orientation.Z0*length_ratio13;

					SymMat<double> S_super(4);
					// 行列式 = ±1になるよう、スケールを調整。
					{
						const Double det = Determinant3(basis);
						const double scale = pow( fabs(det), -1.0/3.0 );
				    	basis *= scale;
						ScaleFactor *= 1.0/scale;

						// S_super = tmat*basis*transpose(tmat*basis);
						// 格子基底簡約。線形独立性が成立しているので、ここのチェックは通るはず。
				    	NRMat<int> tmat(4,3);
				    	if( !put_Selling_reduced_dim_3(putMetricTensor(basis), S_super, tmat, 1.0e-8) ) continue;
				    	moveSmallerDiagonalLeftUpper(S_super, tmat);

						// 格子基底簡約の結果を反映する (Delone簡約基底) => basis is Delone reduced.
				    	basis = mprod( put_transform_matrix_row4to3(tmat), basis );

						if( Determinant3(basis) < 0.0 ) basis *= -1.;
						assert( fabs(Determinant3(basis) - 1.0) < 1.0e-7 );
					}

					// ブラベー格子決定。
			    	BravaisLattice brav_lat;
			      	brav_lat.setParam(false, cons.resolutionForBLD*2.0);
			    	brav_lat.execute(S_super, cons.ABC_AXIS, cons.RH_AXIS);

			    	for(size_t i=0; i<NUM_BT; i++)
			    	{
				    	const vector<LatticeFigureOfMeritToCheckSymmetry>& result = brav_lat.putResult(BravaisType::Enum(i));
				    	for(size_t j=0; j<result.size(); j++)
				    	{
				    		LatticeFigureOfMerit LatFom(result[j].putBravaisType(), result[j].putReducedForm(), trans_mat);

				    		// The conventional cell is given by the row vectors in new_basis, where new_tmat*new_basis = basis.
				    		NRMat<int> new_tmat = mprod(mprod(result[j].putTransformToOriginalLattice(), put_transform_matrix_row4to3(result[j].putReducedForm().second)), trans_mat);

					    	// det(basis)>0のため、det(new_tmat) > 0 とする。
					    	if( Determinant3(new_tmat) < 0 ) new_tmat *= -1;
					    	assert( 1 <= Determinant3(new_tmat) && Determinant3(new_tmat) <= 4 );

					    	const FracMat<int> inv_new_tmat = FInverse3(FracMat<int>(new_tmat));
					    	const NRMat<double> new_basis_dbl( mprod(inv_new_tmat.mat, basis) * (1.0/inv_new_tmat.denom) );

					    	ans[i].resize( ans[i].size() + 1 );
							for(Int4 i2=0; i2<3; i2++)
							{
								for(Int4 j2=0; j2<3; j2++)
								{
									ans[i].back().result_of_abinitio.reciprocal_lattice_basis[i2][j2] = new_basis_dbl[i2][j2];
								}
							}
					    	ans[i].back().inv_tmat_to_prim = inv_new_tmat;
					    	ans[i].back().result_of_abinitio.ScaleFactor = ScaleFactor;
				    	}
					}
				}
			}
		}
	}
}


Indexing::Indexing()
{
}

Indexing::~Indexing()
{
}

void Indexing::set(const conditions& cons, const ProfileData& Data)
{
	// m_data.clear();
	m_data2.clear();

	// DataSet tray;
	DataSet2 tray2;
	OrientationVector vec;

	m_max_tan_sigma = 0.0;
	m_max_bwidth = 0.0;

	const vector<ProfileData::DataSet>& KikuchiData = Data.putData();
	m_get_unitcell_scale = Data.putFlagToGetUnitcellScale();
	m_wavelength = Data.putWavelength();
	for(size_t i=0; i<KikuchiData.size(); i++)
	{
		const double& phi = KikuchiData[i].phi;
		const double& sigma = KikuchiData[i].sigma;

		const double cos_phi = cos(phi);
		const double sin_phi = sin(phi);
		const double cos_sigma = cos(sigma);
		const double sin_sigma = sin(sigma);

		vec.X0 = -cos_sigma*cos_phi;
		vec.Y0 = -cos_sigma*sin_phi;
		vec.Z0 = sin_sigma;

		vec.err_X0 = sqrt(square(cos_sigma*sin_phi) + square(sin_sigma*cos_phi))*cons.ErrorForAngles_radian;
		vec.err_Y0 = sqrt(square(cos_sigma*cos_phi) + square(sin_sigma*sin_phi))*cons.ErrorForAngles_radian;
		vec.err_Z0 = fabs(cos_sigma)*cons.ErrorForAngles_radian;
		tray2.setOrientationVector(vec);

		const double tan_sigma = sin_sigma / cos_sigma;
		tray2.X = cos_phi * tan_sigma;
		tray2.Y = sin_phi * tan_sigma;
		tray2.err_X = sqrt(square(sin_phi*tan_sigma) + square(cos_phi/(cos_sigma*cos_sigma)))*cons.ErrorForAngles_radian;
		tray2.err_Y = sqrt(square(cos_phi*tan_sigma) + square(sin_phi/(cos_sigma*cos_sigma)))*cons.ErrorForAngles_radian;

		if( m_get_unitcell_scale )
		{
			const double sigmab = KikuchiData[i].sigma_begin;
			const double sigmae = KikuchiData[i].sigma_end;
			assert( sigmab < sigmae );
			const double btheta  = sigmae - sigmab;

			tray2.d_star_obs = sin(btheta*0.5)*2.0/m_wavelength;
			tray2.err_d_star_obs = sqrt(2.0) * cos(btheta*0.5) / m_wavelength * cons.ErrorForAngles_radian;
			assert( tray2.d_star_obs > 0.0 );

			const double tan_sigmab = tan(sigmab);
			const double tan_sigmae = tan(sigmae);
			const double inv_cos_sigmab_square = 1.+tan_sigmab*tan_sigmab;
			const double inv_cos_sigmae_square = 1.+tan_sigmae*tan_sigmae;

			tray2.BandWidth = tan_sigmae - tan_sigmab;
			tray2.err_BandWidth = sqrt(square(inv_cos_sigmab_square) + square(inv_cos_sigmae_square)) * cons.ErrorForAngles_radian;
		}
		else
		{
			tray2.d_star_obs = 0.0;
			tray2.err_d_star_obs = 0.0;
			tray2.BandWidth = 0.0;
			tray2.err_BandWidth = 0.0;
		}

		m_data2.push_back(tray2);
		const double tan_sigma_err = tan_sigma + cons.ErrorForAngles_radian/(cos_sigma*cos_sigma);
		const double bwidth_err = tray2.BandWidth + tray2.err_BandWidth;

		if ( m_max_tan_sigma < tan_sigma_err ) m_max_tan_sigma = tan_sigma_err;
		if ( m_max_bwidth < bwidth_err ) m_max_bwidth = bwidth_err;
	}

stringstream strstream;
strstream << "***                ***\n";
strstream << "*** Input data set ***\n";
strstream << "***                ***\n";
strstream << "# Number of data = " << m_data2.size() << "\n";
strstream << "# Kikuchi-line coordinates on the screen" << "\n";
strstream << "# phi, sigma, X, error, Y, error, band-width, error\n";
strstream.setf(ios::fixed);
strstream.precision(5);
for(size_t i=0; i<m_data2.size(); i++){
	strstream << setw(12) << KikuchiData[i].phi*(180./M_PI);
	strstream << setw(12) << KikuchiData[i].sigma*(180./M_PI);
	strstream << setw(12) << m_data2[i].X;
	strstream << setw(10) << m_data2[i].err_X;
	strstream << setw(12) << m_data2[i].Y;
	strstream << setw(10) << m_data2[i].err_Y;
	strstream << setw(12) << m_data2[i].BandWidth;
	strstream << setw(10) << m_data2[i].err_BandWidth << "\n";
}
strstream << "\n";
strstream << "# -cos(sigma)*cos(phi), error, -cos(sigma)*sin(phi), error, sin(sigma), error, d*(=1/d_spacing), error\n";
strstream.setf(ios::fixed);
strstream.precision(5);
for(size_t i=0; i<m_data2.size(); i++){
	const OrientationVector& data = m_data2[i].putOrientationVector();
	strstream << setw(12) << data.X0;
	strstream << setw(10) << data.err_X0;
	strstream << setw(12) << data.Y0;
	strstream << setw(10) << data.err_Y0;
	strstream << setw(12) << data.Z0;
	strstream << setw(10) << data.err_Z0;
	strstream << setw(12) << m_data2[i].d_star_obs;
	strstream << setw(10) << m_data2[i].err_d_star_obs << "\n";
}
strstream << "\n";
strstream << "Radius of Kikuchi pattern R = " << m_max_tan_sigma << " (used for computation of figure of merit)\n";
if ( m_get_unitcell_scale )
strstream << "Upper bound h of band-widths (on screen) = " << m_max_bwidth << "\n";
strstream << endl;
ZLOG_INFO( strstream.str() );

}

static bool equiv(const NRMat<FittingParameter>& data1, const NRMat<FittingParameter>& data2, const double resol)
{
	for(size_t i=0; i<3; i++){
		for (size_t j=0; j<3; j++){
			if( fabs(data1[i][j].value-data2[i][j].value) > resol*max(fabs(data1[i][j].value),fabs(data2[i][j].value)) )
			{
				return false;
			}
		}
	}
	return true;
}

void Indexing::run(const conditions& cons)
{
	// 2次元部分格子の基底候補を選出。
	vector<PairOfReciprocalVectors> pairs_tray;
	find_pair_of_reciprocal_vectors(m_data2, m_wavelength, m_get_unitcell_scale, pairs_tray);
/*
for(size_t i=0; i<pairs_tray.size(); i++)
{
	cout << pairs_tray[i].vec1.index_in_data << " ";
	cout << pairs_tray[i].vec2.index_in_data << " ";
	cout << pairs_tray[i].vec1.length << " ";
	cout << pairs_tray[i].vec2.length << ": ";
	cout << pairs_tray[i].ScaleFactor.value << "(";
	cout << pairs_tray[i].ScaleFactor.error << ")\n";
}
cout << "\n";
*/
ZLOG_INFO( "Number of detected pairs of reciprocal lattice vectors = " + num2str( pairs_tray.size() ) + "\n\n" );

	// 基底候補を選出＋格子基底簡約+ブラベータイプ決定。
	// 各格子のq値の小さなhklを準備しておく。
	find_basis_of_reciprocal_lattice(pairs_tray, cons, m_data2, m_get_unitcell_scale, m_result);

stringstream strstream;
strstream << "***                                            ***\n";
strstream << "*** Number of detected solutions (3D lattices) ***\n";
strstream << "***                                            ***\n";
	for(size_t i=0; i<m_result.size(); i++)
	{
strstream << "(" << i+1 << ") The number of candidates for " + BravaisType::ToString(BravaisType::Enum(i), cons.ABC_AXIS) + " : " << m_result[i].size() << "\n";
	}
strstream << "\n";
ZLOG_INFO( strstream.str() );

	this->refineParameters(cons);

	// ソート
	for(size_t i=0; i<m_result.size(); i++)
	{
		sort(m_result[i].begin(), m_result[i].end());
	}

	for(size_t i=0; i<m_result.size(); i++)
	{
		// m_result[i][j].result_of_refinement.reciprocal_lattice_basis[*][*].value (3*3 行列)の各成分が全部だいたい一緒なら、下位の m_result[i][j].output_flagをfalseとする。
		for(size_t j=0; j<m_result[i].size(); j++)
		{
			for(size_t j2=j+1; j2<m_result[i].size(); j2++)
			{
				if( equiv(m_result[i][j].result_of_refinement.reciprocal_lattice_basis,
							m_result[i][j2].result_of_refinement.reciprocal_lattice_basis, cons.resolutionForOutput) )
				{
					m_result[i][j2].output_flag = false;
				}
			}
		}
	}
}

void Indexing::refineParameters(const conditions& cons)
{
ZLOG_INFO( "Carrying out refinement of parameters ... (solutions with FOM under " + num2str(cons.lower_threshold_for_FOM) + " are removed)\n" );
	int index;
	for(size_t i=0; i<m_result.size(); i++)
	{
		index = 0;
		for(size_t j=0; j<m_result[i].size(); j++)
		{
			m_result[i][j].result_of_abinitio.PCShift[0] = 0.0;
			m_result[i][j].result_of_abinitio.PCShift[1] = 0.0;
			m_result[i][j].result_of_abinitio.PCShift[2] = 0.0;
			m_result[i][j].result_of_abinitio.S_covar = 0.0;

			// m_data(入力値)にhklをつける。
			m_result[i][j].result_of_abinitio.makeHKL_for_observed_data(m_data2, m_max_tan_sigma, m_get_unitcell_scale, m_max_bwidth, m_wavelength,
																		cons.num_hkl, cons.max_hkl, m_result[i][j].inv_tmat_to_prim);

			// All the fitting results are saved in result.
			// In m_result[i][j].result_of_abinitio, only the value of chi_square is modified.
			double chi_square_init;
			indexing_result_type result;
			bool succeed_in_refinement = m_result[i][j].result_of_abinitio.refineParameters(cons, BravaisType::Enum(i), m_data2, m_get_unitcell_scale, m_wavelength,
																							chi_square_init, result);
			m_result[i][j].result_of_abinitio.chi_square = chi_square_init;
			m_result[i][j].result_of_refinement = m_result[i][j].result_of_abinitio;

			if( succeed_in_refinement )
			{
				static const int max_itnum = 15;
				int itnum = 0;
				while( itnum++ < max_itnum )
				{
					// hklをつける。
					result.makeHKL_for_observed_data(m_data2, m_max_tan_sigma, m_get_unitcell_scale, m_max_bwidth, m_wavelength,
													 cons.num_hkl, cons.max_hkl, m_result[i][j].inv_tmat_to_prim);

					if( result.number_of_computed_lines <= 0 ) break;
					if( m_result[i][j].result_of_refinement.figure_of_merit >= result.figure_of_merit )
					{
						if( itnum <= 1 ) m_result[i][j].result_of_refinement = result;
						break;
					}
					m_result[i][j].result_of_refinement = result;

					if( !m_result[i][j].result_of_refinement.refineParameters(cons, BravaisType::Enum(i), m_data2, m_get_unitcell_scale, m_wavelength, chi_square_init, result) )
					{
						break;
					}
				}
			}

			// メモリ使用量を節約するため、しきい値より下回る結果は保存しない。
			if( m_result[i][j].result_of_refinement.figure_of_merit >= cons.lower_threshold_for_FOM )
			{
				m_result[i][index++] = m_result[i][j];
			}
		}
		ZLOG_INFO( "(" + num2str(i+1) + ") The number of candidates for " + BravaisType::ToString(BravaisType::Enum(i), cons.ABC_AXIS) + " : " + num2str( m_result[i].size() ) + " -> " + num2str( index ) + "\n" );
		m_result[i].resize(index);
	}
ZLOG_INFO( "Done.\n" );
}


void Indexing::print(const string& fileName, const conditions& cons) const
{
	static const size_t NUM_BT = 14;
	assert( m_result.size() == NUM_BT );

ZLOG_INFO( "Outputting: " + fileName + "\n\n" );
	ofstream ofs(fileName.c_str());
	if(!ofs) {
		ZLOG_ERROR( "[ERROR!] Cannot open file: " + fileName );
		return;
	}

	ofs << "###\n";
	ofs << "### Radius of Kikuchi pattern R = " << m_max_tan_sigma << " (used for computation of figure of merit)\n";
	if ( m_get_unitcell_scale )
	{
		ofs << "### Upper bound h of band-widths (on screen) = " << m_max_bwidth << endl;
	}
	ofs << "###\n";
	ofs << "###\n";
	ofs << "### Bravais-type: number of indexing solutions (best figure of merit)\n";
	ofs << "###\n\n";
	for(size_t i=0; i<NUM_BT; i++)
	{
		ofs << "(" << i+1 << ") " + BravaisType::ToString(BravaisType::Enum(i), cons.ABC_AXIS) + " : " << m_result[i].size();
		if( !m_result[i].empty() ) ofs << " (" << m_result[i][0].result_of_refinement.figure_of_merit << ")";
		ofs << "\n";
	}
	ofs << "\n";

	for(int i=NUM_BT-1; i>=0; i--)
	{
		ofs << "###                " << setw(20) << "###\n";
		ofs << "### Candidates for " << setw(15) << BravaisType::ToString(BravaisType::Enum(i), cons.ABC_AXIS) << setw(5) << "###\n";
		ofs << "###                " << setw(21) << "###\n\n";
		for(size_t j=0; j<m_result[i].size(); j++)
		{
			if( m_result[i][j].result_of_refinement.figure_of_merit < cons.lower_threshold_for_FOM ) break;
			if( !m_result[i][j].output_flag ) continue;
			ofs << "###" << setw(15) << " No." << j+1 << setw(20) << "###\n";
			ofs << "# Bravais type\n";
			ofs << setw(16) << BravaisType::ToString(BravaisType::Enum(i), cons.ABC_AXIS) << "\n";
			m_result[i][j].toText(ofs, cons, m_data2, m_wavelength, m_get_unitcell_scale);
		}
	}
}


bool Indexing::indexing_result_type::refineParameters(const conditions& cons, const BravaisType::Enum& bravais_type,
		const vector<Indexing::DataSet2>& Data, const bool& get_unitcell_scale,
		const double& wavelength, double& chi_square_init, indexing_result_type& result) const
{
	// Set all the observed data assigned to an hkl are fit.
	vector< pair< VecDat3<int>, char> > hkl_to_fit;
	vector<double> ydata, ydata_err;

	assert( Data.size() == hkls.size() );
	for(size_t k=0; k<hkls.size(); k++)
	{
		if( !hkls[k].isIndexed ) continue;

		const VecDat3<int> hkl(hkls[k].hkl[0], hkls[k].hkl[1], hkls[k].hkl[2]);
		hkl_to_fit.push_back( make_pair(hkl, 0) );
		ydata.push_back(Data[k].X);
		ydata_err.push_back(Data[k].err_X);

		hkl_to_fit.push_back( make_pair(hkl, 1) );
		ydata.push_back(Data[k].Y);
		ydata_err.push_back(Data[k].err_Y);

		if( get_unitcell_scale )
		{
			hkl_to_fit.push_back( make_pair(hkl, 2) );
			ydata.push_back(Data[k].BandWidth);
			ydata_err.push_back(Data[k].err_BandWidth);
		}
	}

	const Vec_BOOL nxfit( hkl_to_fit.size(), true); // Fit by using all the observed data.

	// Fitting of projection centers, lattice parameters, direction of the unitcell.
	EBSDIndexingModel EIModel;
	EIModel.setWaveLength(wavelength);
	EIModel.setConstraint(BravaisType(bravais_type, cons.ABC_AXIS, cons.RH_AXIS), cons.fitProjectionCenterShift, get_unitcell_scale);

	// Parameters are not fit in this case.
	if( (int)hkl_to_fit.size() <= EIModel.putNumberOfIndependentParam() ) return false;
//	cout << "Number of observed points: "<<  hkl_to_fit.size() << endl;
//	cout << "Number of independent parameters: " <<EIModel.putNumberOfIndependentParam()<<endl;

	const NRMat<double> cp_initial_basis = chToDouble(this->reciprocal_lattice_basis);
	vector<double> init_param;
	if( get_unitcell_scale )
	{
		if( !EBSDIndexingModel::putInitialParameters(BravaisType(bravais_type, cons.ABC_AXIS, cons.RH_AXIS), PCShift, ScaleFactor.value, cp_initial_basis, init_param) ) return false;
	}
	else if( !EBSDIndexingModel::putInitialParameters(BravaisType(bravais_type, cons.ABC_AXIS, cons.RH_AXIS), PCShift, 0.0, cp_initial_basis, init_param) ) return false;

	static const bool output_view_flag = false;
	result = *this;
	pair<bool, ZErrorMessage> ans = EIModel.setFittedParam(hkl_to_fit, ydata, ydata_err, nxfit, output_view_flag, 0.0, 15, init_param, chi_square_init, result.chi_square);
	if( !ans.first ) return false;

//	cout<<"Chi-square: " << chisq_result <<endl;
	EIModel.putResult(result.PCShift, result.EulerAngles, result.ScaleFactor, result.reciprocal_lattice_basis, result.S_covar);

	// 行列式=det(reciprocal_lattice_basis)が初期値と一致するよう、スケールを調整。
	NRMat<double> cp_refined_basis = chToDouble( result.reciprocal_lattice_basis );
	const double det_init = Determinant3(cp_initial_basis);
	const double det_refn = Determinant3(cp_refined_basis);
	assert( det_refn != 0.0 );
	const double scale = pow( fabs(det_init/det_refn), 1.0/3.0);
	assert( fabs(det_init - det_refn*scale*scale*scale) <= fabs(det_init*1.0e-8) );

	for(Int4 i=0; i<3; i++)
	{
		for(Int4 j=0; j<3; j++)
		{
			result.reciprocal_lattice_basis[i][j] *= scale;
		}
	}
	result.S_covar *= scale*scale*scale*scale;
	result.ScaleFactor *= 1.0/scale;
/*
	if( use_the_band_widths )
	{
		assert( hkls.size() == Data_d_star_obs.size() );

		double numer_sum = 0.;
		double denom_sum = 0.;
		for(size_t k=0; k<hkls.size(); k++)
		{
			if( !hkls[k].isIndexed ) continue;

			const NRVec<double> V( right_act(hkls[k].hkl, cp_refined_basis) );
			const double d_star_cal = sqrt( inner_product(V, V) );

			numer_sum += d_star_cal * Data_d_star_obs[k].value / square( Data_d_star_obs[k].error );
			denom_sum += square( d_star_cal / Data_d_star_obs[k].error );
		}

		result.ScaleFactor.value = numer_sum/denom_sum;
		result.ScaleFactor.error = 1.0/sqrt(denom_sum);
	}
*/
	return true;
}

void Indexing::indexing_result_type::putBandParameters(const NRVec<int>& hkl, const double& wavelength, double& Xcal, double& Ycal, double& BandWidth) const
{
	const NRMat<double> cp_basis_dbl = chToDouble(reciprocal_lattice_basis);

	NRVec<double> V( right_act(hkl, cp_basis_dbl) );

	const double& x = V[0];
	const double& y = V[1];
	const double& z = V[2];

	const double xyz_len2 = x*x + y*y;
	const double xyz_len3 = xyz_len2 + z*z;
	const double INV_XX_YY = 1.0/xyz_len2;
	const double DeltaZp1 = 1. - PCShift[2].value;

	// X/Z = x*( (x*Delta_x + y*Delta_y - z*(1-Delta_z) )/(x^2+y^2),
	// Y/Z = y*( (x*Delta_x + y*Delta_y - z*(1-Delta_z) )/(x^2+y^2),
	const double Xcal0 = -x*z*INV_XX_YY;
	const double Ycal0 = -y*z*INV_XX_YY;
	const double lp0 = (x*PCShift[0].value + y*PCShift[1].value - z*DeltaZp1)*INV_XX_YY;
	Xcal = x*lp0;
	Ycal = y*lp0;

	const Double d_star_cal = ScaleFactor.value*sqrt(xyz_len3);

	// band-width = (tan(sigma + theta) - tan(sigma - theta))*(1 - Delta_z)
	//            = 2*theta*(1 - Delta_z)/cos^2(sigma) = 2*theta*(1 - Delta_z)*(1 + Xcal^2 + Ycal^2),
	//      sigma = arctan(sqrt(Xcal^2 + Ycal^2)),
	//    2 theta = 2 sin theta = wavelength * d_star_cal.
	const double tan_sigma_square = square(Xcal0) + square(Ycal0);
	BandWidth = wavelength*d_star_cal*(1.0-PCShift[2].value)*(1.0+tan_sigma_square);
}

// Returns true, if all bands are indexed. Otherwise, returns false.
bool Indexing::indexing_result_type::computeFigureOfMerit(const vector<Indexing::hkl_cal>& hkl_cal_data, const size_t& num_hkl,
								const vector<Indexing::DataSet2>& data_set, const double& max_tan_sigma,
								const bool& use_the_band_widths, const double& max_bwidth, const double& wavelength,
								const int& max_hkl, indexing_result_type& result)
{
	assert( hkl_cal_data.size() >= num_hkl );

	vector<hkl_for_data>& hkls = result.hkls;
	int& number_of_computed_lines = result.number_of_computed_lines;
	double& figure_of_merit = result.figure_of_merit;

	vector<double> Xcals(num_hkl);
	vector<double> Ycals(num_hkl);
	vector<double> BandWidths_cal(num_hkl);
	{
		for(size_t n=0; n<num_hkl; n++)
		{
			result.putBandParameters(hkl_cal_data[n].hkl, wavelength, Xcals[n], Ycals[n], BandWidths_cal[n]);
		}
	}

	hkls.clear();
	hkls.resize(data_set.size());
	double sum_dist = 0.0;
	number_of_computed_lines = 0;
	for(size_t k=0; k<data_set.size(); k++)
	{
		const Indexing::DataSet2& data2 = data_set[k];

		double min_dist = max_tan_sigma;	// = Radius of the EBSD pattern.
		int min_dist_index = num_hkl;
		for(size_t n=0; n<num_hkl; n++)
		{
			const double d1 = data2.X - Xcals[n];
			const double d2 = data2.Y - Ycals[n];

			double dist = square(d1) + square(d2);
			int multi_pos = 1;
			if(use_the_band_widths)
			{
				const int multi = iround_half_up<int>( data2.BandWidth / BandWidths_cal[n] );
				multi_pos = (multi>1?min(multi, max_hkl):1);// multi_pos <= 0 is not allowed.

				const double d3 = data2.BandWidth - BandWidths_cal[n];	// multi_pos = 1のhklを優先
				dist += square(d3);
			}

			if( min_dist > dist )
			{
				min_dist = dist;
				min_dist_index = n;
				hkls[k].hkl[0] = hkl_cal_data[n].hkl[0]*multi_pos;
				hkls[k].hkl[1] = hkl_cal_data[n].hkl[1]*multi_pos;
				hkls[k].hkl[2] = hkl_cal_data[n].hkl[2]*multi_pos;
			}

			if( d1*d1 + d2*d2 <= square(data2.err_X) + square(data2.err_Y) )
			{
				hkls[k].isIndexed = true;	// data[k] is given a index.
			}
		}
		sum_dist += sqrt(min_dist);
		if( hkls[k].isIndexed ) number_of_computed_lines = max(min_dist_index, number_of_computed_lines);
		else number_of_computed_lines = num_hkl;
	}
	sum_dist /= data_set.size();

	bool ans = false;
	if( number_of_computed_lines < (int)num_hkl )
	{
		number_of_computed_lines++;
		ans = true;
	}

	if ( use_the_band_widths )
	{
		static const double gamma_1_3 = 2.6789385347077476337;  // = Gamma(1/3) https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function
		const double epsilon_Q = gamma_1_3 *pow(max_tan_sigma*max_tan_sigma* max_bwidth/(36.*number_of_computed_lines), 1./3.);
		figure_of_merit = epsilon_Q/sum_dist;
	}
	else{
		const double epsilon_Q = 0.5 * max_tan_sigma * sqrt(M_PI/number_of_computed_lines);
		figure_of_merit = epsilon_Q/sum_dist;
	}
	return ans;
}


void Indexing::indexing_result_type::makeHKL_for_observed_data(const vector<Indexing::DataSet2>& data_set, const double& max_tan_sigma,
								const bool& use_the_band_widths, const double& max_two_sin_theta, const double& wavelength,
								const size_t& num_hkl, const int& max_hkl, const FracMat<int>& inv_tmat_to_prim)
{
	const NRMat<double> cp_basis_dbl = chToDouble(reciprocal_lattice_basis);

	vector< pair<NRVec<int>, double> > shortest_vec_a;
	if( use_the_band_widths )
	{
		// band-widthから推定されるd*-value(=1/d: d-spacing)がどの観測値よりも大きいhklはindexingの候補から除く。
    	// 今、単位胞の体積が1なので、多くて約973(=2/3*PI*60^{3/2})個のhklが出るはず。
		get_short_vectors(putMetricTensor(cp_basis_dbl), min(60.0, square( max_two_sin_theta / (ScaleFactor.value * wavelength) )), shortest_vec_a);
	}
	else
	{
		// 今、単位胞の体積が1なので、約973(=2/3*PI*60^{3/2})個のhklが出るはず。
		get_short_vectors(putMetricTensor(cp_basis_dbl), 60.0, shortest_vec_a);
	}

	vector<Indexing::hkl_cal> hkl_cal_data;
	for(size_t i=0; i<shortest_vec_a.size(); i++)
	{
		const NRVec<int>& hkl = shortest_vec_a[i].first;
		if ( abs(hkl[0]) > max_hkl || abs(hkl[1]) > max_hkl || abs(hkl[2]) > max_hkl ) continue;

		const NRVec<int> hkl2 = right_act(hkl, inv_tmat_to_prim.mat);
		if( abs(hkl2[0]) % abs(inv_tmat_to_prim.denom) != 0
				|| abs(hkl2[1]) % abs(inv_tmat_to_prim.denom) != 0
				|| abs(hkl2[2]) % abs(inv_tmat_to_prim.denom) != 0 ) continue;
		if ( GCD_Euclid(GCD_Euclid(hkl2[0], hkl2[1]), hkl2[2]) != abs(inv_tmat_to_prim.denom) ) continue;

		const NRVec<double> V( right_act(hkl, cp_basis_dbl) );

		if( (V[0]*V[0] + V[1]*V[1]) * max_tan_sigma * max_tan_sigma < V[2]*V[2] ) continue; // |tan(sigma)| がどの観測値よりも大きいhklはindexingの候補から除く。

		hkl_cal_data.resize( hkl_cal_data.size() + 1 );
		hkl_cal_data.back().hkl = (V[2]>0.0?hkl:hkl*(-1));
		hkl_cal_data.back().q = inner_product(V, V);
	}
	sort( hkl_cal_data.begin(), hkl_cal_data.end() );

	const size_t max_num_hkl = min(hkl_cal_data.size(), num_hkl);
	computeFigureOfMerit(hkl_cal_data, max_num_hkl, data_set, max_tan_sigma,
							use_the_band_widths, max_two_sin_theta, wavelength, max_hkl, *this);
}

// Returns error of s1*a1+s2*a2 when cov11, cov12, cov13 are covariance between a1, a2.
inline double put_propagation_error(const double& cov11, const double& cov12, const double& cov22,
const double& s1, const double& s2)
{
	assert( square(s1)*cov11 + square(s2)*cov22 + 2.0*s1*s2*cov12 >= 0.0 );
	return sqrt( square(s1)*cov11 + square(s2)*cov22 + 2.0*s1*s2*cov12 );
}

void Indexing::result_type::toText(ofstream& ofs, const conditions& cons, const vector<DataSet2>& Data_XY,
		const double& wavelength, const bool& use_the_band_width) const
{
	static const double DegRad = 180.0/M_PI;
	const double ErrCoef = 1.0/(cons.ErrorForAngles_radian*DegRad);

	ofs << "# a : b : c  alpha  beta  gamma (degree) scale_factor (before refinement)\n";

	vector<Double> lattice_consts;
	SymMat<double> lattice_consts_covar(6);
	const NRMat<double> cp_basis_dbl_ab_initio = chToDouble(result_of_abinitio.reciprocal_lattice_basis);
	const NRMat<double> cp_basis_dbl_refined   = chToDouble(result_of_refinement.reciprocal_lattice_basis);
	calLatticeConstant(putMetricTensor(cp_basis_dbl_ab_initio), SymMat<double>(6, 0.0), lattice_consts, lattice_consts_covar);

	ofs.setf(ios::fixed);
	ofs.setf(ios::right);
	ofs.precision(4);

	ofs.width(10); ofs << lattice_consts[0];
	ofs.width(10); ofs << lattice_consts[1];
	ofs.width(10); ofs << lattice_consts[2];
	ofs.width(10); ofs << lattice_consts[3];
	ofs.width(10); ofs << lattice_consts[4];
	ofs.width(10); ofs << lattice_consts[5];
	ofs.width(12); ofs << result_of_abinitio.ScaleFactor.value << "\n";

	if( use_the_band_width )
	{
		ofs << "# a : b : c  alpha  beta  gamma (degree) scale_factor, a/c, b/c, a, b, c (after refinement)\n";
	}
	else{
		ofs << "# a : b : c  alpha  beta  gamma (degree) scale_factor, a/c, b/c (after refinement)\n";
	}
	calLatticeConstant(putMetricTensor(cp_basis_dbl_refined), result_of_refinement.S_covar, lattice_consts, lattice_consts_covar);

	ofs.width(10); ofs << lattice_consts[0];
	ofs.width(10); ofs << lattice_consts[1];
	ofs.width(10); ofs << lattice_consts[2];
	ofs.width(10); ofs << lattice_consts[3];
	ofs.width(10); ofs << lattice_consts[4];
	ofs.width(10); ofs << lattice_consts[5];

	ofs.width(12); ofs << result_of_refinement.ScaleFactor.value;
	ofs.width(15); ofs << lattice_consts[0]/lattice_consts[2];
	ofs.width(10); ofs << lattice_consts[1]/lattice_consts[2];
	if( use_the_band_width )
	{
		ofs.width(10); ofs << lattice_consts[0]/result_of_refinement.ScaleFactor.value;
		ofs.width(10); ofs << lattice_consts[1]/result_of_refinement.ScaleFactor.value;
		ofs.width(10); ofs << lattice_consts[2]/result_of_refinement.ScaleFactor.value;
	}
	ofs << "\n";

	ofs << "# propagation errors when the errors of the input angles are assumed to be within 1 deg.\n";
	ofs.width(10); ofs << sqrt( lattice_consts_covar(0,0) ) * ErrCoef;
	ofs.width(10); ofs << sqrt( lattice_consts_covar(1,1) ) * ErrCoef;
	ofs.width(10); ofs << sqrt( lattice_consts_covar(2,2) ) * ErrCoef;
	ofs.width(10); ofs << sqrt( lattice_consts_covar(3,3) ) * ErrCoef;
	ofs.width(10); ofs << sqrt( lattice_consts_covar(4,4) ) * ErrCoef;
	ofs.width(10); ofs << sqrt( lattice_consts_covar(5,5) ) * ErrCoef;
	ofs.width(12); ofs << result_of_refinement.ScaleFactor.error * ErrCoef;

	// Error of a/c
	ofs.width(15); ofs << put_propagation_error(lattice_consts_covar(0,0), lattice_consts_covar(0,2), lattice_consts_covar(2,2),
		1.0/lattice_consts[2], -lattice_consts[0]/square(lattice_consts[2]) ) * ErrCoef;
	// Error of b/c
	ofs.width(10); ofs << put_propagation_error(lattice_consts_covar(1,1), lattice_consts_covar(1,2), lattice_consts_covar(2,2),
		1.0/lattice_consts[2], -lattice_consts[1]/square(lattice_consts[2]) ) * ErrCoef;
	if( use_the_band_width )
	{
		// Error of a/scale
		ofs.width(10); ofs << put_propagation_error(lattice_consts_covar(0,0), 0.0, square(result_of_refinement.ScaleFactor.error),
			1.0/result_of_refinement.ScaleFactor.value, -lattice_consts[0]/square(result_of_refinement.ScaleFactor.value) ) * ErrCoef;
		// Error of b/scale
		ofs.width(10); ofs << put_propagation_error(lattice_consts_covar(1,1), 0.0, square(result_of_refinement.ScaleFactor.error),
			1.0/result_of_refinement.ScaleFactor.value, -lattice_consts[1]/square(result_of_refinement.ScaleFactor.value) ) * ErrCoef;
		// Error of c/scale
		ofs.width(10); ofs << put_propagation_error(lattice_consts_covar(2,2), 0.0, square(result_of_refinement.ScaleFactor.error),
			1.0/result_of_refinement.ScaleFactor.value, -lattice_consts[2]/square(result_of_refinement.ScaleFactor.value) ) * ErrCoef;
	}
	ofs << "\n";

	ofs << "# Buerger-reduced reciprocal_lattice basis (before refinement)\n";
	ofs.width(7); ofs << "a1 = (";
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[0][0].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[0][1].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[0][2].value << ")\n";
	ofs.width(7); ofs << "a2 = (";
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[1][0].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[1][1].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[1][2].value << ")\n";
	ofs.width(7); ofs << "a3 = (";
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[2][0].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[2][1].value;
	ofs.width(10); ofs << result_of_abinitio.reciprocal_lattice_basis[2][2].value << ")\n";

	ofs << "# Buerger-reduced reciprocal_lattice basis, propagation errors  (after refinement)\n";
	ofs.width(7); ofs << "a1 = (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][0].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][1].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][2].value << ") (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][0].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][1].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[0][2].error * ErrCoef << ")\n";
	ofs.width(7); ofs << "a2 = (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][0].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][1].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][2].value << ") (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][0].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][1].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[1][2].error * ErrCoef << ")\n";
	ofs.width(7); ofs << "a3 = (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][0].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][1].value;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][2].value << ") (";
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][0].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][1].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.reciprocal_lattice_basis[2][2].error * ErrCoef << ")\n";

	ofs << "# Euler angles: theta1, theta2, theta3, Error_theta1, Error_theta2, Error_theta3 (deg) (after refinement)\n";
	ofs.width(10); ofs << result_of_refinement.EulerAngles[0].value * DegRad;
	ofs.width(10); ofs << result_of_refinement.EulerAngles[1].value * DegRad;
	ofs.width(10); ofs << result_of_refinement.EulerAngles[2].value * DegRad;
	ofs.width(10); ofs << result_of_refinement.EulerAngles[0].error * DegRad * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.EulerAngles[1].error * DegRad * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.EulerAngles[2].error * DegRad * ErrCoef << "\n";

	ofs << "# Projection center shifts: Delta_x, Delta_y, Delta_z, Error_Delta_x, Error_Delta_y, Error_Delta_z, \n";
	ofs.width(10); ofs << result_of_refinement.PCShift[0].value;
	ofs.width(10); ofs << result_of_refinement.PCShift[1].value;
	ofs.width(10); ofs << result_of_refinement.PCShift[2].value;
	ofs.width(12); ofs << result_of_refinement.PCShift[0].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.PCShift[1].error * ErrCoef;
	ofs.width(10); ofs << result_of_refinement.PCShift[2].error * ErrCoef << "\n\n";

	ofs << "# Number of computed bands\n";
	ofs.width(15); ofs << result_of_abinitio.number_of_computed_lines;
	ofs.width(15); ofs << result_of_refinement.number_of_computed_lines << endl;

	ofs << "# Figure of merit at the beginning and the end of the refinement\n";
	ofs.width(15); ofs << result_of_abinitio.figure_of_merit;
	ofs.width(15); ofs << result_of_refinement.figure_of_merit << endl;

	ofs << "# Chi-squares at the beginning and the end of the refinement\n";
	ofs.width(10); ofs << result_of_abinitio.chi_square / square(ErrCoef);
	ofs.width(10); ofs << result_of_refinement.chi_square / square(ErrCoef) << "\n";

//	ofs << "# Number of observed data that was not assigned with a Miller index at the beginning and the end of the refinement:\n";
//	ofs.width(10); ofs << count_is_not_indexed(hkls_for_data);
//	ofs.width(10); ofs << count_is_not_indexed(hkls_for_data_refined) << "\n";

	ofs << "# Indexing with the parameters before refinement\n";
	if (use_the_band_width)
	{
		ofs << "# Band_number, Miller_index, X_cal, Y_cal, X_obs, Y_obs, distance, good_fit?, band_width_cal, band_width_obs\n";
	}else
	{
		ofs << "# Band_number, Miller_index, X_cal, Y_cal, X_obs, Y_obs, distance, good_fit?\n";
	}

	assert( result_of_abinitio.hkls.size() == Data_XY.size() );

	NRVec<int> hkl(3);
	double Xcal, Ycal, BandWidth;
	for(size_t k=0; k<result_of_abinitio.hkls.size(); k++)
	{
		ofs.width(5); ofs << k+1 << (cons.output_tex_format?" &":"");

		hkl[0] = result_of_abinitio.hkls[k].hkl[0];
		hkl[1] = result_of_abinitio.hkls[k].hkl[1];
		hkl[2] = result_of_abinitio.hkls[k].hkl[2];
		ofs.width(5); ofs << hkl[0] << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << hkl[1] << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << hkl[2] << (cons.output_tex_format?" &":"");

		result_of_abinitio.putBandParameters(hkl, wavelength, Xcal, Ycal, BandWidth);

		ofs.width(10); ofs << Xcal << (cons.output_tex_format?" &":"");
		ofs.width(10); ofs << Ycal << (cons.output_tex_format?" &":"");

		ofs.width(10); ofs << Data_XY[k].X << (cons.output_tex_format?" &":"");
		ofs.width(10); ofs << Data_XY[k].Y << (cons.output_tex_format?" &":"");

		ofs.width(10); ofs << sqrt(square(Xcal - Data_XY[k].X) + square(Ycal - Data_XY[k].Y)) << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << result_of_abinitio.hkls[k].isIndexed << (cons.output_tex_format?" &":"");
		if (use_the_band_width)
		{
			ofs.width(10); ofs << BandWidth << (cons.output_tex_format?" &":"");
			ofs.width(10); ofs << Data_XY[k].BandWidth;
		}
		ofs  << (cons.output_tex_format?" \\\\":"") <<"\n";
	}
	ofs << "\n";

	ofs << "# Indexing with the parameters after refinement\n";
	if (use_the_band_width)
	{
		ofs << "# Band_number, Miller_index, X_cal, Y_cal, X_obs, Y_obs, distance, good_fit?, band_width_cal, band_width_obs\n";
	}
	else
	{
		ofs << "# Band_number, Miller_index, X_cal, Y_cal, X_obs, Y_obs, distance, good_fit?\n";
	}

//	ofs << "# line_number, phi_obs(rad), sigma_obs(rad), X_obs, Y_obs, Miller index, phi_cal(rad), sigma_cal(rad), phi_diff, sigma_diff, X_cal, Y_cal, X_cal_diff, Y_cal_diff\n";

	assert( result_of_refinement.hkls.size() == Data_XY.size() );
	for(size_t k=0; k<result_of_refinement.hkls.size(); k++)
	{
		ofs.width(5); ofs << k+1 << (cons.output_tex_format?" &":"");

		hkl[0] = result_of_refinement.hkls[k].hkl[0];
		hkl[1] = result_of_refinement.hkls[k].hkl[1];
		hkl[2] = result_of_refinement.hkls[k].hkl[2];
		ofs.width(5); ofs << hkl[0] << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << hkl[1] << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << hkl[2] << (cons.output_tex_format?" &":"");

		NRVec<double> V(3);
		V[0] = hkl[0]*result_of_refinement.reciprocal_lattice_basis[0][0].value + hkl[1]*result_of_refinement.reciprocal_lattice_basis[1][0].value + hkl[2]*result_of_refinement.reciprocal_lattice_basis[2][0].value;
		V[1] = hkl[0]*result_of_refinement.reciprocal_lattice_basis[0][1].value + hkl[1]*result_of_refinement.reciprocal_lattice_basis[1][1].value + hkl[2]*result_of_refinement.reciprocal_lattice_basis[2][1].value;
		V[2] = hkl[0]*result_of_refinement.reciprocal_lattice_basis[0][2].value + hkl[1]*result_of_refinement.reciprocal_lattice_basis[1][2].value + hkl[2]*result_of_refinement.reciprocal_lattice_basis[2][2].value;
		assert( V[2] >= 0.0 );

		result_of_refinement.putBandParameters(hkl, wavelength, Xcal, Ycal, BandWidth);

		ofs.width(10); ofs << Xcal << (cons.output_tex_format?" &":"");
		ofs.width(10); ofs << Ycal << (cons.output_tex_format?" &":"");

		ofs.width(10); ofs << Data_XY[k].X << (cons.output_tex_format?" &":"");
		ofs.width(10); ofs << Data_XY[k].Y << (cons.output_tex_format?" &":"");

		ofs.width(10); ofs << sqrt(square(Xcal - Data_XY[k].X) + square(Ycal - Data_XY[k].Y)) << (cons.output_tex_format?" &":"");
		ofs.width(5); ofs << result_of_refinement.hkls[k].isIndexed << (cons.output_tex_format?" &":"");
		if (use_the_band_width)
		{
			ofs.width(10); ofs << (cons.output_tex_format?" &":"") << BandWidth;
			ofs.width(10); ofs << (cons.output_tex_format?" &":"") << Data_XY[k].BandWidth;

		}
		ofs  << (cons.output_tex_format?" \\\\":"") <<"\n";
	}
	ofs << "\n";
}
