/* Laruence Brown
   SMU Mathematics
   Math 3316
   21 November 2016 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "matrix.h"
#include "fcn.hpp"

 /*
	 Create an array of 400 evenly-spaced T values over the interval [800, 1200] K. Output this
	to disk as the file Temp.txt.
	• Create an array of 600 evenly-spaced t values from t = 1 second up to t = 48 hours.
	Output this to disk as the file time.txt.
	• Create a 400 × 600 array that contains C(0.002, t, T). Output this to disk as the file
	C2mm.txt.
	• Create a 400 × 600 array that contains C(0.004, t, T). Output this to disk as the file
	C4mm.txt.
*/

using namespace std;


// function prototypes
double carbon(const double x, const double t, const double T,	
				const double rtol, const double atol);

// main routine
int main(int argc, char* argv[]) {

	//settings
	double rtol = -pow(10,11);
	double atol = -pow(10,15);

	//test vectors
	vector<double> Ts = Linspace(800, 1200, 400);
	vector<double> ts = Linspace(1, 172800, 600);

	Matrix C2mm(400,600);
	Matrix C4mm(400,600);
	for (int i=0; i<400; i++){
		for (int j=0; j<600; j++){
			C2mm.M[j][i]=carbon(0.002, ts[j], Ts[i], rtol, atol);
			C4mm.M[j][i]=carbon(0.004, ts[j], Ts[i], rtol, atol);
		}
	}

	vector<double> C2mm_800K;
	vector<double> C2mm_900K;
	vector<double> C2mm_1000K;
	vector<double> C2mm_1100K;
	vector<double> C2mm_1200K;

	vector<double> C4mm_800K;
	vector<double> C4mm_900K;
	vector<double> C4mm_1000K;
	vector<double> C4mm_1100K;
	vector<double> C4mm_1200K;

	for (int j=0; j<600; j++){
		C2mm_800K.push_back(carbon(0.004, ts[j], 800, rtol, atol));
		C2mm_900K.push_back(carbon(0.004, ts[j], 900, rtol, atol));
		C2mm_1000K.push_back(carbon(0.004, ts[j], 1000, rtol, atol));
		C2mm_1100K.push_back(carbon(0.004, ts[j], 1100, rtol, atol));
		C2mm_1200K.push_back(carbon(0.004, ts[j], 1200, rtol, atol));

		C4mm_800K.push_back(carbon(0.004, ts[j], 800, rtol, atol));
		C4mm_900K.push_back(carbon(0.004, ts[j], 900, rtol, atol));
		C4mm_1000K.push_back(carbon(0.004, ts[j], 1000, rtol, atol));
		C4mm_1100K.push_back(carbon(0.004, ts[j], 1100, rtol, atol));
		C4mm_1200K.push_back(carbon(0.004, ts[j], 1200, rtol, atol));
	}

	//write data
	write(ts, "time.txt");
	write(Ts, "Temp.txt");

	C2mm.write("C2mm.txt");
	C4mm.write("C4mm.txt");

	write(C2mm_800K, "C2mm_800K.txt");
	write(C2mm_900K, "C2mm_900K.txt");
	write(C2mm_1000K, "C2mm_1000K.txt");
	write(C2mm_1100K, "C2mm_1100K.txt");
	write(C2mm_1200K, "C2mm_1200K.txt");

	write(C4mm_800K, "C4mm_800K.txt");
	write(C4mm_900K, "C4mm_900K.txt");
	write(C4mm_1000K, "C4mm_1000K.txt");
	write(C4mm_1100K, "C4mm_1100K.txt");
	write(C4mm_1200K, "C4mm_1200K.txt");


	return 0;
}