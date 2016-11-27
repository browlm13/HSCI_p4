/* 
   Laruence Brown
   SMU Mathematics
   Math 3316
   1 December 2016 
*/

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include "fcn.hpp"

using namespace std;

// This routine numerically approximates the integral
//    int_a^b fun(x) dx
// using the composite Gaussian quadrature rule with x points per 
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
//           x       number of nodes per subinterval
//
// outputs:  F       value of numerical integral
// 

//modify weights and points for given x
void get_xw(int x, vector<double> &xi, vector<double> &w){

  double xi_n[9];
  double w_n[9];

  //weights and points
  double xi_1[1] = {0};
  double w_1[1] = {2};

  double xi_2[2] = {-0.57735026918963, 0.57735026918963};
  double w_2[2] = {1,1};

  double xi_3[3] = {-0.77459666924148, 0 ,0.77459666924148};
  double w_3[3] = {0.55555555555555, 0.88888888888889, 0.55555555555555};

  double xi_4[4] = {-0.86113631159405,-0.33998104358486,
      0.33998104358486,0.86113631159405};
  double w_4[4] = {0.34785484513745,0.65214515486254,
      0.65214515486255,0.34785484513745};

  double xi_5[5] = {-0.90617984593866,-0.53846931010568,
      0.00000000000000,0.53846931010568,0.90617984593866};
  double w_5[5] = {0.23692688505619,0.47862867049937,
      0.56888888888889,0.47862867049937,0.23692688505619};

  double xi_6[6] = {-0.93246951420315,-0.66120938646626,
      -0.23861918608320, 0.23861918608320, 0.66120938646626,
      0.93246951420315};
  double w_6[6] = {0.17132449237917, 0.36076157304814,
      0.46791393457269, 0.46791393457269, 0.36076157304814,
      0.17132449237917};

  if (x == 1){
    copy(begin(xi_1), end(xi_1), begin(xi_n));
    copy(begin(w_1), end(w_1), begin(w_n));
  }
  else if (x == 2){
    copy(begin(xi_2), end(xi_2), begin(xi_n));
    copy(begin(w_2), end(w_2), begin(w_n));
  }
  else if (x == 3){
    copy(begin(xi_3), end(xi_3), begin(xi_n));
    copy(begin(w_3), end(w_3), begin(w_n));
  }
  else if (x == 4){
    copy(begin(xi_4), end(xi_4), begin(xi_n));
    copy(begin(w_4), end(w_4), begin(w_n));
  }
  else if (x == 5){
    copy(begin(xi_5), end(xi_5), begin(xi_n));
    copy(begin(w_5), end(w_5), begin(w_n));
  }
  else if (x == 6){
    copy(begin(xi_6), end(xi_6), begin(xi_n));
    copy(begin(w_6), end(w_6), begin(w_n));
  }
  else {
    cout << "Illegal x value..." << endl;
  }

  //update passed vectors
  xi.insert(xi.begin() , xi_n , xi_n + x) ; 
  w.insert(w.begin() , w_n , w_n + x) ; 

}

double composite_int_gen(Fcn& f, const double a, 
			const double b, const int n, int x) {

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
  vector<double> xi;
  vector<double> w;
  get_xw(x, xi, w);

  // initialize result
  double F = 0.0;

  // loop over subintervals, accumulating result
  //double xmid, node1, node2, node3, node4;
  double xmid;

  double nodes[x];
  for (int i=0; i<n; i++) {
   
    // determine evaluation points within subinterval
    xmid  = a + (i+0.5)*h;

    for (int i=0; i<x; i++){
      nodes[i] = xmid + 0.5*h*xi[i];
    }

    // add Gauss-n approximation on this subinterval to result
    for (int i=0; i<x; i++){
      F += w[i]*f(nodes[i]);
    }

  } // end loop

  // return final result
  return (0.5*h*F);

} // end of function