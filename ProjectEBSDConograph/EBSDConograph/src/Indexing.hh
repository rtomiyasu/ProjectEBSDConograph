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

#ifndef INDEXING_HH_
#define INDEXING_HH_

#include <vector>
#include "ProfileData.hh"
#include "conograph_codes/bravais_type/BravaisType.hh"
#include "conograph_codes/bravais_type/enumAxis.hh"
#include "conograph_codes/lattice_symmetry/HKL_Q.hh"
#include "conograph_codes/utility_data_structure/nrutil_nr.hh"
#include "conograph_codes/utility_data_structure/SymMat.hh"
#include "marquardt_codes/zparam/FittingParameter.hh"

using namespace std;

class Indexing
{
public:
	/*!
	 * 計算に必要なパラメータ
	 */
	struct conditions {
		// ブラベー格子決定後の軸を選択するパラメータ。
		eABCaxis ABC_AXIS;
		eRHaxis RH_AXIS;

		// 0: quick search, 1: exhaustive search.
		int search_level;
		bool fitProjectionCenterShift[3];

		double ErrorForAngles_radian;
		size_t num_hkl; // 指数付けのため生成するミラー指数の数。
		int max_hkl;  // 指数付けのため生成するhklに対し、許容される絶対値の最大値。
	    double resolutionForBLD;    // 格子ベクトルの長さ|l|に対し、|l|*resolの差異を許す。Used for Bravais-lattice determination.
	    double resolutionForOutput;    // 格子ベクトルの各成分(l1, l2, l3)に対し、|li|*resolの差異を許す。Used for selection of output solutions.
		double toleranceOnBandWidth; // when s1, s2 are the unit-cell scale from distinct bands, (1/toleranceOnBandWidth)*s2 <= s1 <= toleranceOnBandWidth*s2 is allowed.
		double lower_threshold_for_FOM;

		bool output_tex_format;
	};

	class OrientationVector {
	public:
		double X0;
		double Y0;
		double Z0;

		double err_X0;
		double err_Y0;
		double err_Z0;

	};

	class DataSet2 {
	private:
		// (-cos \sigma cos \phi, -cos \sigma sin \phi, sin \sigma)
		OrientationVector m_orientation; // vector of length 1.

	public:
		// (cos \phi tan \sigma, sin \phi tan \sigma)
		double X;
		double Y;
		double err_X;
		double err_Y;

		double BandWidth;
		double err_BandWidth;
		double d_star_obs;
		double err_d_star_obs;
		
		inline void setOrientationVector(const OrientationVector& arg) { m_orientation = arg; };
		inline const OrientationVector& putOrientationVector() const { return m_orientation; };
	};

	class InfoOfReciprocalVector
	{
	public:
		int index_in_data;
		OrientationVector orientation;
		double length;
	};

	class PairOfReciprocalVectors
	{
	public:
		// Pair of reciprocal lattice vectors used to construct a 3D lattice.
		InfoOfReciprocalVector vec1; 	// vec1.first: index in m_data, vec1.second: length of the corresponding reciprocal lattice vector.
		InfoOfReciprocalVector vec2; 	// vec2.first: index in m_data, vec2.second: length of the corresponding reciprocal lattice vector.
		FittingParameter ScaleFactor;  // Set by using band-widths.

//		bool operator<(const PairOfReciprocalVectors& rhs) const
//		{
//			if( vec1.index_in_data < rhs.vec1.index_in_data ) return true;
//			if( vec1.index_in_data > rhs.vec1.index_in_data ) return false;
//			if( vec2.index_in_data < rhs.vec2.index_in_data ) return true;
//			if( vec2.index_in_data > rhs.vec2.index_in_data ) return false;
//			if( vec1.length < rhs.vec1.length ) return true;
//			if( vec1.length > rhs.vec1.length ) return false;
//			if( vec2.length < rhs.vec2.length ) return true;
//			return false;
//		}
	};

	class hkl_cal
	{
	public:
		NRVec<int> hkl; // Miller index.
		double q;	// q-values.

		hkl_cal(): hkl(3,0), q(0.0){};
		bool operator<(const hkl_cal& rhs) const { return this->q < rhs.q; };
	};

	class hkl_for_data
	{
	private:
		enum{ dim = 3 };
	public:
		int hkl[dim]; // Miller index.
		bool isIndexed;

		hkl_for_data(): isIndexed(false){};
	};

	class indexing_result_type
	{
	private:
		enum{ dim = 3 };
	public:
		int number_of_computed_lines;
		double chi_square;
		double figure_of_merit;

    	// reciprocal_bases:
		NRMat<FittingParameter> reciprocal_lattice_basis;
    	SymMat<double> S_covar;

    	FittingParameter ScaleFactor;
    	FittingParameter PCShift[dim];

    	// ( cos(theta) sin(theta)  0)   (1          0          0 )   ( cos(phi) sin(phi)  0)
    	// (-sin(theta) cos(theta)  0) * (0  cos(sigma) sin(sigma)) * (-sin(phi) cos(phi)  0)
    	// (        0           0   1)   (0 -sin(sigma) cos(sigma))   (        0       0   1)
    	FittingParameter EulerAngles[dim]; // =(theta, sigma, phi)
		vector<hkl_for_data> hkls;

		indexing_result_type(): chi_square(0.0), figure_of_merit(0.0), reciprocal_lattice_basis(3,3), S_covar(6, 0.0), ScaleFactor(0.0){ PCShift[0]=0.0; PCShift[1]=0.0; PCShift[2]=0.0; };

		void putBandParameters(const NRVec<int>& hkl, const double& wavelength, double& Xcal, double& Ycal, double& BandWidth) const;

		// Returns true, if all bands are indexed. Otherwise, returns false.
		static bool computeFigureOfMerit(const vector<Indexing::hkl_cal>& hkl_cal_data, const size_t& num_hkl,
										const vector<Indexing::DataSet2>& data_set, const double& max_tan_sigma,
										const bool& get_unitcell_scale, const double& max_two_sin_theta, const double& wavelength,
										const int& max_hkl, indexing_result_type& result);

		bool refineParameters(const conditions& cons, const BravaisType::Enum& bravais_type,
								const vector<Indexing::DataSet2>& Data, const bool& get_unitcell_scale,
								const double& wavelength, double& chi_square_init, indexing_result_type& result) const;

		void makeHKL_for_observed_data(const vector<Indexing::DataSet2>& data_set, const double& max_tan_sigma,
								const bool& get_unitcell_scale, const double& max_sin_beta_2, const double& wavelength,
								const size_t& num_hkl, const int& max_hkl, const FracMat<int>& inv_tmat_to_prim);
	};

	class result_type
	{
	public:
//    	NRMat<int> tmat_to_prim;
    	FracMat<int> inv_tmat_to_prim;

    	class indexing_result_type result_of_abinitio;	// Parameters before refinement;
    	class indexing_result_type result_of_refinement;

		bool output_flag;

    	result_type(): inv_tmat_to_prim(NRMat<int>(3,3)), output_flag(true){};
		void toText(ofstream& ofs, const conditions& cons, const vector<DataSet2>& Data_XY,
				    const double& wavelength, const bool& get_unitcell_scale) const;

		bool operator<(const result_type& rhs) const { return this->result_of_refinement.figure_of_merit > rhs.result_of_refinement.figure_of_merit; };
	};

	Indexing();
	virtual ~Indexing();

	void set(const conditions&, const ProfileData& Data);
	void run(const conditions&);
	void refineParameters(const conditions& cons);
	void print(const string& fname, const conditions& cons) const;
	inline const vector< vector<result_type> >& putResult() const { return m_result; };

private:
	// Data on (-cos \sigma cos \phi, -cos \sigma sin \phi, sin \sigma)
	// vector<DataSet> m_data;

	// Data on (cos \phi tan \sigma, sin \phi tan \sigma)
	vector<DataSet2> m_data2;

	bool m_get_unitcell_scale;
	double m_wavelength;

	double m_max_tan_sigma;
	double m_max_bwidth;

	// [0]: Triclinic, [1]: Monoclinic_P, [2]: Monoclinic_B
	// [3]: Orthorhombic_P, [4]: Orthorhombic_C, [5]: Orthorhombic_I, [6]: Orthorhombic_F,
	// [7]: Tetragonal_P, [8]: Tetragonal_I, [9]: Rhombohedral, [10]: Hexagonal,
	// [11]: Cubic_P, [12]: Cubic_I, [13]: Cubic_F.
	vector< vector<result_type> > m_result;
};

#endif /* INDEXING_HH_ */
