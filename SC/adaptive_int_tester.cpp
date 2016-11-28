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

using namespace std;

// function prototypes
double composite_int(Fcn& f, const double a, 
			const double b, const int n);

/*
	-a and b are interval endpoints
	-Ntot = total number of intervals tried
	-R = Integral estimation
	-rtol and atol are used in modling equation

	model equation with...
	|Rn+k(f) − Rn(f)| < rtol |Rn+k(f)| + atol.

	trick is picking smart n and k values
*/

int integrate(Fcn& f, Fcn& g, double a, double b, int m, const int k, const double rtol,
							const double atol, double& R, int& n, int& Ntot){

	double tol = rtol * fabs(composite_int(f,a,b,m+k)) + atol;
	double err = fabs(composite_int(f,a,b,m+k) - composite_int(f,a,b,m));

	//break value
	if( m > 1000){
		Ntot = 1000000;
		return 0;
	}

	if( err > tol ){
		Ntot += m+k;
		m *= int(g(m));
		integrate(f, g, a, b, m, k, rtol, atol, R, n, Ntot);
	}

	//printf("\ntol:%f, \terr:%f\n", tol, err);
	R = composite_int(f,a,b,m);
	Ntot += m+k;
	n = m;
	return 1;




}

void adaptive_int_tester(Fcn& f, Fcn& g, int m, int k, const double a, const double b, const double rtol,
									const double atol, double& R, int& n, int& Ntot){
	//number of intervals
	Ntot = 0;
	integrate(f, g, a, b, m, k, rtol, atol, R, n, Ntot);
	
} // end of function

