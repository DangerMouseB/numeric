//************************************************************************************************************************************************
//
//    Copyright (c) 2011 David Briant - see https://github.com/DangerMouseB
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//************************************************************************************************************************************************


// Provide a C style interface to make calling from Python and VB/VBA easy.
// We allocate and free memory though the client is responsible to initiate freeing of objects.
//
// See other Sobol implementations:
//		Numerical Recipes
// 	Kuo & Joe - http://web.maths.unsw.edu.au/~fkuo/sobol/index.html
//		Jungman - http://www.koders.com/c/fid205FACC04ABC4B63199A5C306F544FAB17D52A01.aspx
//		Burkardt - http://people.sc.fsu.edu/~jburkardt/py_src/sobol/sobol.html
//
// sobol_newSI - creates a SobolInitialiser with a file from Kuo & Joe - new-joe-kuo-6.21201.txt is recommened by Joe & Kuo
// sobol_freeSI
// sobol_newSG
// sobol_freeSG
// sobol_next
// sobol_skip
//
// TODO:
// file
//		error handling - testing the fileformat error detection code
// speed
// 	can speed up sobol_next and sobol_skipWhileDuplicates
// interface
//		sobol_hasDuplicates?
//		sobol_blockNext?


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>		// getline global fn
#include <sstream>	// stringstream
#include "math.h"
#include <set>


#define INVALID_FILE 1
#define SI_FILE_SYNTAX_ERROR 2
#define SEQUENCE_EXHAUSTED 3
#define ALLOCATION_ERROR 4
#define INVALID_ARGUMENTS 5
#define HAS_DUPLICATES 6


typedef unsigned int uint32;

using namespace std;

struct IndexError{};

struct SobolInitialiser {
	uint32 max_d;
	uint32 max_s;
	uint32 data[1];

	inline uint32 *ps(uint32 d) {
		if (max_d == 0 | max_s ==0 | d > max_d) throw IndexError();
		return &data[(d - 1) * (max_s + 2)];
	}

	inline uint32 *pa(uint32 d) {
		if (max_d == 0 | max_s ==0 | d > max_d) throw IndexError();
		return &data[(d - 1) * (max_s + 2) + 1];
	}

	inline uint32 *pm(uint32 d) {
		if (max_d == 0 | max_s ==0 | d > max_d) throw IndexError();
		return &data[(d - 1) * (max_s + 2) + 2];
	}

	inline uint32& m(uint32 d, uint32 s) {
		if (max_d == 0 | max_s ==0 | d > max_d | s > max_s) throw IndexError();
		return data[(d - 1) * (max_s + 2) + 2 + (s - 1)];
	}

	inline uint32 m(uint32 d, uint32 s) const {
		if (max_d == 0 | max_s ==0 | d > max_d | s > max_s) throw IndexError();
		return data[(d - 1) * (max_s + 2) + 2 + (s - 1)];
	}

	inline uint32 s(uint32 d) const {return data[(d - 1) * (max_s + 2)];}
	inline uint32 a(uint32 d) const {return data[(d - 1) * (max_s + 2) + 1];}
};


static double denominator;

struct SobolGenerator {
	// sub-arrays of SG (pX, pS & pV) are indexed rather than offsetted, i.e. base 1 rather than base 0, because we think of dimension 1 and direction number 1
	uint32 D;							// number of dimensions
	uint32 N;							// number of points
	uint32 *pG;						// pG[N] - Gray coding - see NR 20.2 - fn implmented as an array as an optimisation
	uint32 *pX;						// pX[D] - current numerator for each dimension
	uint32 *pS;						// pS[D] - sequence tracker - this allows each dimension to be in a different part of the sequence
	uint32 L;							// max number of bits needed
	uint32 *pV;						// V[L][D] - direction numbers - accessed V(d, l)
	bool exhausted;					// true if any S has reached N

	inline uint32 G(uint32 n) const {return pG[n];}
	inline uint32& X(uint32 d) {return pX[(d-1)];}
	inline uint32 X(uint32 d) const {return pX[(d-1)];}
	inline uint32& S(uint32 d) {return pS[(d-1)];}
	inline uint32 S(uint32 d) const {return pS[(d-1)];}
	inline uint32& V(uint32 d, uint32 l) {return pV[D * (l-1) + (d-1)];}
	inline uint32 V(uint32 d, uint32 l) const {return pV[D * (l-1) + (d-1)];}

};


extern "C" int __stdcall sobol_newSI(char const * const filename, SobolInitialiser ** const o_ppsi) {
	ifstream initfile;
	uint32 d, s, a, max_d, max_s, i;
	string lineString;
	uint32 *m;
	SobolInitialiser *psi;

	denominator = 1 / pow(2.0, 32);

	*o_ppsi = NULL;							// clear the ouput variable

	initfile.open(filename, ios::in);
	if (! initfile.is_open()) return INVALID_FILE;


	// find max d and max s
	initfile.seekg (0, ios::end);
	if (initfile.tellg() < 1024) {
		initfile.seekg (0, ios::beg);
	} else {
		initfile.seekg (-1024, ios::end);
	}

	getline(initfile, lineString);					// skip partial line or header line
	while(getline(initfile, lineString)) {
		std::stringstream line(lineString);
		d = 0;
		line >> d;
		if (line.good() && d != 0) {
			max_d = d;
			line >> max_s;
			if (max_s == 0) return SI_FILE_SYNTAX_ERROR;
		}
	}


	// over-allocate storage
	psi = (SobolInitialiser *)malloc(sizeof(SobolInitialiser) - sizeof(int) + max_d * (max_s + 2) * sizeof(int));
		// -sizeof(int) to account for the int[1] already in the struct def.
		// + 2 for each s and a at the beginning of each line.
	psi->max_d = max_d;
	psi->max_s = max_s;


	// read the data in
	initfile.clear();
	initfile.seekg (0, ios::beg);
	getline(initfile, lineString);					// skip the header line

	try {
		while(getline(initfile, lineString)) {
			std::stringstream line(lineString);
			d = 0;
			line >> d;
			if (line.good() && d != 0 && d <= max_d) {
				line >> s;
				*(psi->ps(d)) = s;
				if (s == 0 | s > max_s) {
					free(psi);
					return SI_FILE_SYNTAX_ERROR;
				}
				line >> a;
				*(psi->pa(d)) = a;
				if (d > 2 && a == 0) {
					free(psi);
					return SI_FILE_SYNTAX_ERROR;
				}
				m = psi->pm(d);
				for (i=0; i<s; i++) {
					line >> m[i];
					if (m[i] == 0) {
						free(psi);
						return SI_FILE_SYNTAX_ERROR;
					}
				}
			}
		}
	} catch (IndexError) {
		free(psi);
		return SI_FILE_SYNTAX_ERROR;
	}

	// clean up
	initfile.close();

	*o_ppsi = psi;
	return 0;
}


extern "C" int __stdcall sobol_freeSI(SobolInitialiser * const psi) {
	free(psi);
	return 0;
}


extern "C" int __stdcall sobol_newSG(uint32 D, uint32 N, SobolInitialiser const * const psi, SobolGenerator ** const o_ppsg) {
	uint32 d, l, n, s, a, k, value, L;

	if (psi == NULL | D == 0 | D > psi->max_d | N == 0) return INVALID_ARGUMENTS;

	SobolInitialiser const &si = *psi;
	L = (uint32)ceil(log((double)N)/log(2.0));
	// allocate storage
	SobolGenerator * psg = (SobolGenerator*)malloc(sizeof(SobolGenerator));
	if (psg == NULL) return ALLOCATION_ERROR;
	SobolGenerator &sg = *psg;
	sg.D = D;
	sg.N = N;
	sg.pG = (uint32*) malloc(sizeof(uint32) * N);
	sg.pX = (uint32*) calloc(D, sizeof(uint32));
	sg.pS = (uint32*) calloc(D, sizeof(uint32));
	sg.L = L;
	sg.pV = (uint32*) calloc(L * D, sizeof(uint32));
	sg.exhausted = false;

	if (sg.pG == NULL | sg.pX == NULL | sg.pS == NULL | sg.pV== NULL) {
		free(sg.pG);
		free(sg.pX);
		free(sg.pS);
		free(sg.pV);
		return ALLOCATION_ERROR;
	}

	// Gray coding - gets the bit that needs to change to get to the next one (hence is a fn not an index)
	sg.pG[0] = 1;
	for (n = 1; n < N; n++) {
		sg.pG[n] = 1;
		value = n;
		while (value & 1) {
			value >>= 1;
			sg.pG[n]++;
		}
	}
//	if (false) {
//		cout << "\nG [";
//		for (n = 0; n < N; n++) {
//			if (n > 0) cout << ",";
//			cout <<sg.pG[n];
//		}
//		cout << "]\n";
//	}

	// calculate the direction numbers
	d = 1;
//	cout << "\n[\n [1,0,0,";
	for (l=1; l<=L; l++) sg.V(d, l) = 1 << (32-l);
//	for (l=1; l<=L; l++) {
//		if (l > 1) cout << ",";
//		cout << sg.V(d, l);
//	}
//	cout << "],\n";
	for (d=2; d<=D; d++) {
		s = si.s(d);
		a = si.a(d);
//		cout <<" [" << d << "," << s <<"," << a << ",";
		if (s >= L) {
			for (l=1; l<=L; l++) {
				//cout << "1(" << d << "," << l << "),";
				sg.V(d, l) = si.m(d, l) << (32-l);
			}
		} else {
			for (l=1; l<=s; l++) {
				//cout << "2(" << d << "," << l << "),";
				sg.V(d, l) = si.m(d, l) << (32-l);
			}
			for (l=s+1; l<=L; l++) {
				//cout << "3(" << d << "," << l << "),";
				sg.V(d,l) = sg.V(d,l-s) ^ (sg.V(d,l-s) >> s);
				for (k=1; k<=s-1; k++) {
					sg.V(d,l) ^= (((a >> (s-1-k)) & 1) * sg.V(d,l-k));
				}
			}
		}
//		for (l=1; l<=L; l++) {
//			if (l > 1) cout << ",";
//			cout << sg.V(d, l);
//		}
//		cout << "]";
//		if (d<D) cout << ",";
//		cout << endl;
	}
//	cout << "]" << endl;


	*o_ppsg = (SobolGenerator *)&sg;
	return 0;
}


extern "C" int __stdcall sobol_freeSG(SobolGenerator * const psg) {
	free(psg->pG);
	free(psg->pX);
	free(psg->pS);
	free(psg->pV);
	free(psg);
	return 0;
}


inline void sobolStep(SobolGenerator &sg, uint32 const &d) {
	sg.X(d) ^= sg.V(d, sg.G((++sg.S(d))-1));			// update X by xoring it with appropiate direction number
}


extern "C" int __stdcall sobol_next(SobolGenerator * const psg, double *oVector) {
	if (psg == NULL | oVector == NULL) return INVALID_ARGUMENTS;
	SobolGenerator &sg = *psg;
	if (sg.exhausted) return SEQUENCE_EXHAUSTED;

	/*
	cout << "\n   S: ";
	for  (uint32 d=1; d<=sg.D; d++) {
		cout << sg.S(d) << ", ";
	}
	cout << endl;

	cout << "   G: ";
	for  (uint32 d=1; d<=sg.D; d++) {
		cout << sg.G((sg.S(d))) << ", ";
	}
	cout << endl;

	cout << "   V: ";
	for  (uint32 d=1; d<=sg.D; d++) {
		cout << sg.V(d, sg.G(sg.S(d))) << ", ";
	}
	cout << endl;

	cout << "oldX: ";
	for  (uint32 d=1; d<=sg.D; d++) {
		cout << sg.X(d) << ", ";
	}
	cout << endl;

	cout << "newX: ";
	for  (uint32 d=1; d<=sg.D; d++) {
		cout << (sg.X(d) ^ sg.V(d, sg.G(sg.S(d)))) << ", ";
	}
	cout << endl;
	*/

	for (uint32 d=1; d<=sg.D; d++) {
		sobolStep(sg, d);
		if (sg.S(d) >= sg.N) sg.exhausted = true;
		oVector[d-1] = sg.X(d) * denominator;
	}

	return 0;
}


extern "C" int __stdcall sobol_skip(SobolGenerator * const psg, uint32 const d, uint32 const n) {
	// allows the user to easily "warm-up" the sequence by skipping n
	if (psg == NULL | d == 0 | d > psg->D | n == 0) return INVALID_ARGUMENTS;
	SobolGenerator &sg = *psg;
	if (sg.exhausted) return SEQUENCE_EXHAUSTED;
	uint32 N = sg.N;
	uint32 d_1 = d - 1;
	uint32 D = sg.D;
	uint32 X = sg.X(d);
	uint32 S = sg.S(d);
	for (uint32 i=0; i<n; i++) {
		//if (sg.S(d) < N) sobolStep(sg, d);									// 2.8817985844141845 to skip 70k numbers in 1000 dimensions
		//sobolStep(sg, d);															// 2.434695197376186
		//uint32 l_1 = sg.pG[(++sg.pS[d_1])-1]-1;						// 1.1444231154179754
		//sg.pX[d_1] ^= sg.pV[D * l_1 + d_1];
		X ^= sg.pV[D * (sg.pG[S]-1) + d_1];								// 0.5769590054970718 with no exhaustion protection
		S++;
		if (S >= N) break;															// 0.6525645302852215 with exhaustion protection
	}
	sg.X(d) = X;
	sg.S(d) = S;
	if (sg.S(d) >= N) sg.exhausted = true;
	return 0;
}


/*
extern "C" int __stdcall sobol_skip(SobolGenerator * const psg, uint32 const d, uint32 const n) {
	// allows the user to easily "warm-up" the sequence by skipping n
	if (psg == NULL | d == 0 | d > psg->D | n == 0) return INVALID_ARGUMENTS;
	SobolGenerator &sg = *psg;
	if (sg.exhausted) return SEQUENCE_EXHAUSTED;
	for (uint32 i=0; i<n; i++) {
		if (sg.S(d) < sg.N) sobolStep(sg, d);									// 2.8817985844141845
	}
	if (sg.S(d) >= sg.N) sg.exhausted = true;
	return 0;
}
*/

bool hasDuplicates(SobolGenerator const &sg, set<uint32> s) {
	s.clear();
	uint32 D = sg.D;
	for (uint32 d=D; d>=1; d--) {
		s.insert(sg.X(d));
		if (s.size() != (D-d+1)) return true;
	}
	return false;
}

extern "C" int __stdcall sobol_skipWhileDuplicates(SobolGenerator * const psg, uint32 const max) {
	uint32 n, i, d;
	set<uint32> s;
	SobolGenerator &sg = *psg;
	if (!hasDuplicates(sg, s)) return 0;
	if (sg.exhausted) return HAS_DUPLICATES;
	for (n=0; n<max; n++) {
		for (uint32 d=1; d<=sg.D; d++) {
			if (sg.S(d) < sg.N) {
				sobolStep(sg, d);
			} else {
				sg.exhausted = true;
			}
		}
		if (!hasDuplicates(sg, s)) return 0;
		if (sg.exhausted) return HAS_DUPLICATES;
	}
	return HAS_DUPLICATES;
}