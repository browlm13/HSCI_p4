/* Laruence Brown
   SMU Mathematics
   Math 3316
   21 November 2016 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "fcn.hpp"

// function prototypes
int adaptive_int(Fcn& f, const double a, const double b, const double rtol,
              const double atol, double& R, int& n, int& Ntot);

// Integrand for erf
class Ierf : public Fcn {
public:
  double operator()(double z) {   // function evaluation
    return (-exp(-pow(z,2)));
  }
};

// D(T), teperature
//6.2 * 10^-7 exp(- (8*10^4)/(8.31T))
double D(double T){return 6.2 * (1/pow(10,7)) * exp(-(8*pow(10,4))/(8.31*T));}

//calls adaptive int for solving integral
double erf(Fcn& ierf, const double y, const double rtol, const double atol){
	double R;
	int n, Ntot;
	R = (2/M_PI) * adaptive_int(ierf, 0, y, rtol, atol, R, n, Ntot);
	return R;
}

double carbon(const double x, const double t, const double T,	const double rtol, const double atol){
	double C_0 = 0.001;
	double C_s = 0.02;

	double parm = x / sqrt(4*t*(D(T)));
	return  C_s - (C_s - C_0) * erf(parm);
}