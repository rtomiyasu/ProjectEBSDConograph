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
#ifndef enumSymmetricOperation_HH_
#define enumSymmetricOperation_HH_

enum eZSymmetricOperation{
		Id, Inv, C2X, C2Y, C2Z, // C2hex, 
		SigmaX, SigmaY, SigmaZ, // Sigmahex, 

		C4Xplus, C4Yplus, C4Zplus, 
		C4Xminus, C4Yminus, C4Zminus, 
		S4Xminus, S4Yminus, S4Zminus, 
		S4Xplus, S4Yplus, S4Zplus, 

		C31plus_rho, C32plus_rho, C33plus_rho, C34plus_rho,
		C31minus_rho, C32minus_rho, C33minus_rho, C34minus_rho,
		S61minus_rho, S62minus_rho, S63minus_rho, S64minus_rho,
		S61plus_rho, S62plus_rho, S63plus_rho, S64plus_rho,
		
		C6plus_hex, C3plus_hex, C3minus_hex, C6minus_hex,
		S3minus_hex, S6minus_hex, S6plus_hex, S3plus_hex,

		C2primeX0_rho, C2primeX1_rho, C2primeY0_rho, C2primeY1_rho, C2primeZ0, C2primeZ1, // C2f, C2d, C2e, C2c, C2b(=C"23), C2a(=C"23)
		SigmaDX0_rho, SigmaDX1_rho, SigmaDY0_rho, SigmaDY1_rho, SigmaDZ0, SigmaDZ1, // sigma_f, sigma_d, sigma_e, sigma_c, sigma_b(sigma_d3), sigma_a(sigma_v3)

		C2primeX0_hex, C2primeX1_hex, C2primeY0_hex, C2primeY1_hex, // C2primeZ0_hex, // C'21, C'22, C'23, C"21, C"22
		SigmaDX0_hex, SigmaDX1_hex, SigmaDY0_hex, SigmaDY1_hex, // SigmaDZ0_hex, // sigma_d1, sigma_v1, sigma_d2, sigma_v2
		
		D2dummy
	};

#endif /*enumSymmetricOperation_HH*/
