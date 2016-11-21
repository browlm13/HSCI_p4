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

// This routine numerically approximates the integral
//    int_a^b fun(x) dx
// using the composite Gaussian quadrature rule with 4 points per 
// subinterval, over n subintervals.  We 
// require that fun have the calling syntax 
//    y = fun(x)
// where y is a double and x is a const double.
//
// Usage: F = composite_Gauss2(fun, a, b, n);
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
  double x1 = -0.86113631159405;
  double x2 = -0.33998104358486;
  double x3 = 0.33998104358486;
  double x4 = 0.86113631159405;
  double w1 = 0.34785484513745;
  double w2 = 0.65214515486254;
  double w3 = 0.65214515486255;
  double w4 = 0.34785484513745;

  // initialize result
  double F = 0.0;

  // loop over subintervals, accumulating result
  double xmid, node1, node2, node3, node4;
  for (int i=0; i<n; i++) {
   
    // determine evaluation points within subinterval
    xmid  = a + (i+0.5)*h;
    node1 = xmid + 0.5*h*x1;
    node2 = xmid + 0.5*h*x2;
    node3 = xmid + 0.5*h*x3;
    node4 = xmid + 0.5*h*x4;


    // add Gauss2 approximation on this subinterval to result
    F += w1*f(node1) + w2*f(node2) + w3*f(node3) + w4*f(node4);

  } // end loop

  // return final result
  return (0.5*h*F);

} // end of function