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

#ifndef _ProfileData_h_
#define _ProfileData_h_

// ProfileData.hh
#include <vector>
#include <string>

using namespace std;

class ProfileData
{
public:
	struct DataSet {
		 // sigma(rad), phi(rad).
		double sigma;
		double phi;
		double sigma_begin;
		double sigma_end;
	};

    ProfileData();
    ~ProfileData();

    string readFile(const string& fname);

    const vector<DataSet>& putData() const { return m_data; };
    const bool& putFlagToGetUnitcellScale()  const { return get_unitcell_scale; };
    const double& putWavelength() const { return wavelength; };

private:
    vector<DataSet> m_data;     // Profile data (entire data)
	bool get_unitcell_scale;
	double wavelength;
};

#endif
