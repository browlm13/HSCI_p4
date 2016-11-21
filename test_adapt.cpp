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
#include "fcn.hpp"

using namespace std;

// function prototypes
int adaptive_int(Fcn& f, const double a, const double b, const double rtol,
              const double atol, double& R, int& n, int& Ntot);

// Integrand
class fcn : public Fcn {
public:
  double c, d;
  double operator()(double x) {   // function evaluation
    return (exp(c*x) + sin(d*x));
  }
  double antiderivative(double x) { // function evaluation
    return (exp(c*x)/c - cos(d*x)/d);
  }
};


// This routine tests the Gauss-4 method on a simple integral
int main(int argc, char* argv[]) {

  // limits of integration
  double a = -3.0;
  double b = 5.0;

  // integrand
  fcn f;
  f.c = 0.5;
  f.d = 25.0;

  // true integral value
  double Itrue = f.antiderivative(b) - f.antiderivative(a);
  printf("\n True Integral = %22.16e\n", Itrue);


  // test the Gauss-4 rule
  cout << "\n adaptive solver:\n";
  cout << "     n             R(f)            relerr      Ntot\n";
  cout << "  ---------------------------------------------------\n";
  //vector<int> n = {20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};
  vector<double> errors;
  //vector<double> hvals(n.size());

  // iterate over n values, computing approximations, error, convergence rate
  double Iapprox; //integral estimation
  int Ntot; //number of intervals tried
  int n;    //final number of intervals used

  double rtol = 1/(10*10*10*10);
  double atol = rtol/ (1000);

  //call adaptive solver
  adaptive_int(f, a, b, rtol, atol, Iapprox, n, Ntot);

  errors.push_back(fabs(Itrue-Iapprox)/fabs(Itrue));

  printf("\t%d\t%22.16e\t%7.1e\t%d\n", n, Iapprox, errors[0], Ntot);


  cout << "  ---------------------------------------------------\n";

}
