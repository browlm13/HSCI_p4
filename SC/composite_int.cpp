/* 
   Laruence Brown
   SMU Mathematics
   Math 3316
   21 November 2016 
*/

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "fcn.hpp"

using namespace std;

// This routine numerically approximates the integral
//    int_a^b fun(x) dx
// using the composite Gaussian quadrature rule with 6 points per 
// subinterval, over n subintervals.  We 
// require that fun have the calling syntax 
//    y = fun(x)
// where y is a double and x is a const double.
//
// Usage: F = composite_int(fun, a, b, n);
//
// inputs:   f       integrand (Fcn object)
//           a       lower limit of integration
//           b       upper limit of integration
//           n       number of subintervals
//
// outputs:  F       value of numerical integral
// 

double composite_int(Fcn& f, const double a, 
			const double b, const int n) {

  // check input arguments
  if (b < a) {
    cerr << "error: illegal interval, b < a\n";
    return 0.0;
  }
  if (n < 1) {
    cerr << "error: illegal number of subintervals, n < 1\n";
    return 0.0;
  }

  // set subinterval width
  double h = (b-a)/n;

  // set nodes/weights defining the quadrature method
  double xi_6[6] = {-0.93246951420315,-0.66120938646626,
      -0.23861918608320, 0.23861918608320, 0.66120938646626,
      0.93246951420315};
  double w_6[6] = {0.17132449237917, 0.36076157304814,
      0.46791393457269, 0.46791393457269, 0.36076157304814,
      0.17132449237917};

  // initialize result
  double F = 0.0;

  // loop over subintervals, accumulating result
  double xmid;
  for (int i=0; i<n; i++) {
    xmid  = a + (i+0.5)*h;
    for (int j=0; j<6; j++)
      F += w_6[j]*f(xmid + 0.5*h*xi_6[j]);
   

  } // end loop

  // return final result
  return (0.5*h*F);

} // end of function