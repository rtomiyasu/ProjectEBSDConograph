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

#include <cmath>
#include "../../conograph_codes/utility_func/zstring.hh"
#include "linear_constraint.hh"

// On output, cmat has a form as follows;
// 1 0 * 0 *
// 0 1 * 0 *
// 0 0 0 1 *
void gauss_elimination_row(vector< vector<double> >& cmat, const double& thred)
{
	const size_t irow = cmat.size();
	if(irow == 0) return;
	const size_t icol = cmat[0].size();

	// Gauss elimination by row.
	size_t row_start=0;
	for(size_t k=0; k<icol; k++)
	{
		// non-zero element among fabs(cmat[j][k]) s.t. j>=k. 
		size_t index = row_start;
		while(index<irow && fabs(cmat[index][k])<thred) index++;
		if(index==irow) continue;

		if(row_start != index)
			for(size_t i=k; i<icol; i++) swap(cmat[row_start][i], cmat[index][i]);

		// Set cmat[row_start][k]=1.0;
		double t = 1.0/cmat[row_start][k];
		for(size_t i=0; i<icol; i++) cmat[row_start][i]*=t;

		for(size_t j=0; j<irow; j++){
			 if(j!=row_start){
				for(size_t i=k+1; i<icol; i++) cmat[j][i] -= cmat[row_start][i]*cmat[j][k];
				 cmat[j][k] = 0.0;
			 }
		}
		row_start++;
	}
	cmat.resize(row_start);	// Eliminate rows with all elements zero.
}

// On output, cmat_implicit is arranged to the least size.
void putExplicitConstraints(const size_t& pn, vector< vector<double> >& cmat_implicit,
		Mat_DP_constr& cmat_explicit, const double& thred)
{
	cmat_explicit.clear();
	cmat_explicit.resize(pn);
	for(size_t k=0; k<pn; k++) cmat_explicit[k].ID = _IDVary;

	if( cmat_implicit.empty() ) return;

	gauss_elimination_row(cmat_implicit, thred);
	const size_t irow2 = cmat_implicit.size();

	index_set<double> tray;
	double t;
	for(size_t k=0, jrow=0; k<pn && jrow<irow2; k++){
		if( fabs(cmat_implicit[jrow][k]) >thred)
		{
			t = -1.0/cmat_implicit[jrow][k];
			for(size_t i=k+1; i<pn; i++){
				tray.element = cmat_implicit[jrow][i]*t;
				if( fabs(tray.element)<thred ) continue;
				tray.index = i;
				cmat_explicit[k].constr.push_back(tray);
			}
			if( cmat_explicit[k].constr.empty() ) cmat_explicit[k].ID = _IDFixed;
			else cmat_explicit[k].ID = _IDDepend;
			jrow++;
		}
	}
}



bool checkIfFirstConstraintsSatisfySecond(const vector< vector<double> >& constr1,
		const vector< vector<double> >& constr2, const double& thred)
{
	const size_t isize = constr1.size();
	vector< vector<double> > constr1_copy = constr1;
	
	constr1_copy.insert( constr1_copy.end(), constr2.begin(), constr2.end() );
	gauss_elimination_row(constr1_copy, thred);

	return constr1_copy.size() == isize;
}




static string chTermToString(const double& coef, const string& par_label)
{
	if( coef == 1.0 ) return par_label;
	return num2str<double>(coef) + " * " + par_label;
}

string chEquationToString(const int index, const Vec_DP_save& cont,
const vector<string>& par_label)
{
	string ans = par_label[index] + " = ";

	if( cont.empty() ) return ans + "0";
	if( cont[0].element < 0.0 ) ans += "- ";
	
	ans += chTermToString( abs( cont[0].element), par_label[cont[0].index] );

	for(size_t l=1; l<cont.size(); l++)
	{
		if( cont[l].element < 0.0 ) ans += " - ";
		else ans += " + ";

		ans += chTermToString( abs( cont[l].element ), par_label[cont[l].index] );
	}
	return ans;
}
