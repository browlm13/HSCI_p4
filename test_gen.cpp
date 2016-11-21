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
double composite_int_gen(Fcn& f, const double a, 
			const double b, const int n, int order);

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

//get order of convergece
double magnitude_change(double arg1, double arg2){
 return (arg2 == 0) ? 0 : (int)(1 + log10(fabs(arg1)) ) - (int)(1 + log10(fabs(arg2)) );
}

//mean
double get_mean(vector<double> sample){
  double mean = 0;
  for (int i=0; i<sample.size(); i++){
    if (isinf(sample[i]))
      mean += (mean/(sample.size()-1));         //[todo] fix model
    else
      mean += sample[i];
  }
  return mean /= sample.size();
}

//standard devation
double get_sd(vector<double> sample){
  double mean = get_mean(sample);

  //sd = sqrt((E (x-mean)^2)/n-1)
  double variance = 0;
  for (int i=0; i<sample.size(); i++){
    if (isinf(sample[i]))                       //[todo] fix model
      sample[i] = (mean/(sample.size()-1));
    variance += pow((sample[i] - mean),2);
  }
  variance /= sample.size() - 1;
  return sqrt(variance);

}

// This routine tests the Gauss-5 method on a simple integral
int main(int argc, char* argv[]) {

  // limits of integration
  double a = -3.0;
  double b = 5.0;

  // integrand
  fcn f;
  f.c = 0.5;
  f.d = 25.0;

  //set order for guass for all tests
  for (int order=1; order<7; order++){
    // true integral value
    double Itrue = f.antiderivative(b) - f.antiderivative(a);
    printf("\n True Integral = %22.16e\n", Itrue);


    // test the Gauss-4 rule
    printf("\n Gauss-%d approximation:\n", order);
    cout << "     n             R(f)            relerr    conv rate\n";
    cout << "  ---------------------------------------------------\n";
    vector<int> n = {10, 20, 40, 80, 160, 320};
    //vector<int> n = {5, 10, 15, 20, 25, 30, 35, 40, 80, 160, 320};
    vector<double> errors(n.size());
    vector<double> hvals(n.size());
    vector<double> conv_rates;
    vector<double> error_order_imprv;

    // iterate over n values, computing approximations, error, convergence rate
    double Iapprox;
    for (int i=0; i<n.size(); i++) {

      printf("   %6i", n[i]);

      Iapprox = composite_int_gen(f, a, b, n[i], order);
      errors[i] = fabs(Itrue-Iapprox)/fabs(Itrue);
      hvals[i] = (b-a)/n[i];
      if (i == 0) 
        printf("  %22.16e  %7.1e     ----\n", Iapprox, errors[i]);
      else{
        conv_rates.push_back( (log(errors[i-1]) - log(errors[i])) / (log(hvals[i-1]) - log(hvals[i])) );
        error_order_imprv.push_back(magnitude_change( errors[i-1], errors[i]));
        printf("  %22.16e  %7.1e   %f\n", Iapprox, errors[i], conv_rates[i-1]);
      }
    }

    //display converges rate anylisis
    //mean convergence rate
    double mean = get_mean(conv_rates);
    //sd
    double sd = get_sd(conv_rates);

    //mean relitive error change in magnitude
    //double mean_error_change = get_mean(error_order_imprv);
    //double sd_error_change = get_sd(error_order_imprv);

    printf("\n\tmean convergence: %f, sd: %f", mean, sd);
    //printf("\n\tmean error exp change: %f, sd: %f\n", mean_error_change, sd_error_change);

    cout << "  ---------------------------------------------------\n";

  }

}


