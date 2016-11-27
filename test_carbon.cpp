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

/*
//create a new vector of linearly space data
vector<double> Linspace(double a, double b, size_t n){
  if (n<2) std::cerr << "Linspace::length must be > 1\n";
  std::vector<double> v(n);
  
  double h = (b-a)/(n-1);
  for (size_t i=0; i<n; i++)
    v[i] = a + i*h;

  return v;
}

//write a vector to a file
int write(std::vector<double> v, const char *outfile){ 

  // open output file
  FILE *fptr = NULL;
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    cerr << "Write:: error, unable to open " << outfile << " for writing\n";
    return 1;
  }

  // print data to file
  for (size_t i=0; i<v.size(); i++) {
      fprintf(fptr, "  %.16g", v[i]);
  }

  // close output file and return
  fclose(fptr);
  return 0;
}
*/

// main routine
int main(int argc, char* argv[]) {

	//settings
	double rtol = -pow(10,11);
	double atol = -pow(10,15);

	//test vectors
	vector<double> Ts = Linspace(800, 1200, 400);
	vector<double> ts = Linspace(1, 172800, 600);

	//vector<double> C2mm;
	//vector<double> C4mm;

	//vector<double> C_002_t_800;

	Matrix C2mm(600,400);
	Matrix C4mm(600,400);
	for (int i=0; i<600; i++){
		for (int j=0; j<400; j++){
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