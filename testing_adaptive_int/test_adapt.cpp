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
int adaptive_int(Fcn& f, Fcn& g, const double a, const double b, const double rtol,
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

// Interval Updating Function
class Fupd : public Fcn {
public:
  double c;
  double operator()(double n) {   // function evaluation
    return (pow(n,c));
  }
};

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
int write(std::vector<double> v, const char *outfile) { 

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


// This routine tests the Gauss-4 method on a simple integral
int main(int argc, char* argv[]) {

  // limits of integration
  double a = -3.0;
  double b = 5.0;

  // integrand
  fcn f;
  f.c = 0.5;
  f.d = 25.0;

  //
  // true integral value
  //
  double Itrue = f.antiderivative(b) - f.antiderivative(a);


    //nupdating function x^c, c values (NUM TOL for now)
  const int NUM_FS = 10;
  vector<double> n_up_f_consts = Linspace(0.4, 1.0,NUM_FS);


  //
  //    run trials
  //
  for (int j=0; j<NUM_FS; j++){
        //g() interval number update function
    Fupd g;
    g.c = n_up_f_consts[j];
    cout << "\n\n\n\n Interval updating function x^n n=%f:" << n_up_f_consts[j] << endl << endl;


        //ugly caluclating tolerence with - exp list
    const int NUM_TOL = 15;
    vector<double>rtol_exps = Linspace(1,14,NUM_TOL);

    //
    //    data vectors
    //
    vector<double> rtols;
    for (int i=0; i<NUM_TOL; i++){
      rtols.push_back(1/pow(10,rtol_exps[i]));
    }
    vector<double> errors;
    vector<int> Ntots;
    vector<int> ns;


    for (int i=0; i<NUM_TOL; i++){

      //excepted tolerance
      double rtol = rtols[i];
      double atol = rtol/(1000);
      double tolerance = rtol* fabs(Itrue) + atol;

      printf("\n\tTrue Integral = %22.16e\n", Itrue);
      printf("\tTolerence = %7.1e\n", tolerance);


      // test the Gauss-4 rule
      cout << "\n adaptive solver:\n";
      cout << "     n             R(f)            relerr      Ntot\n";
      cout << "  ---------------------------------------------------\n";
      //vector<int> n = {20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};

      // iterate over n values, computing approximations, error, convergence rate
      double Iapprox; //integral estimation
      int Ntot; //number of intervals tried
      int n;    //final number of intervals used


      //
      //    call adaptive solver
      //
      adaptive_int(f, g, a, b, rtol, atol, Iapprox, n, Ntot);



      //update data
      errors.push_back(fabs(Itrue-Iapprox)/fabs(Itrue));
      Ntots.push_back(Ntot);
      ns.push_back(n);

      printf("    %d\t  %22.12e  %7.1e   %d\n", ns[i], Iapprox, errors[i], Ntots[i]);

      cout << "  ---------------------------------------------------\n";
    }
  }

  //
  //    write data vectors to files
  //
  //write(rtols, "tols.txt");
  //write(errors, "errors.txt");
  //write(Ntots, "Ntots.txt");
  //write(ns, "ns.txt");
}
