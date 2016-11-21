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
double composite_int_gen(Fcn& f, const double a, 
			const double b, const int n, int order);

/*
	-a and b are interval endpoints
	-Ntot = total number of intervals tried
	-R = Integral estimation
	-rtol and atol are used in modling equation

	model equation with...
	|Rn+k(f) − Rn(f)| < rtol |Rn+k(f)| + atol.

	trick is picking smart n and k values
*/

void integrate(Fcn& f, double a, double b, int m, const double rtol,
							const double atol, double& R, int& n, int& Ntot){
	m *= 4;
	int k = m*8;
	double tol = rtol * fabs(composite_int_gen(f,a,b,m,6));
	double err = fabs(composite_int_gen(f,a,b,m+k,6)) - fabs(composite_int_gen(f,a,b,m,6));

	if( err > tol ){
		integrate(f, a, b, m, rtol, atol, R, n, Ntot);
		Ntot += m+k;
	}
	else{
		R = composite_int_gen(f,a,b,m,6);
		n = m;
	}

}

void adaptive_int(Fcn& f, const double a, const double b, const double rtol,
							const double atol, double& R, int& n, int& Ntot){
	//number of intervals
	int m = 1;
	integrate(f, a, b, m, rtol, atol, R, n, Ntot);
	
} // end of function

