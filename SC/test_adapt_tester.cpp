/* 
   Laruence Brown
   SMU Mathematics
   Math 3316
   21 November 2016 
*/

/* 
        Finding corret m and k values, and interval size update function
*/


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "matrix.h"
#include "fcn.hpp"

using namespace std;

// function prototypes
int adaptive_int_tester(Fcn& f, Fcn& g, int m, int k, const double a, const double b, const double rtol,
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



//sum vector
int sum_vector(vector<int> v){
  int sum =0;
  for (int i=0; i<v.size(); i++)
    sum += v[i];
  return sum;
}

//get lowest value from vector
int low_element(vector<int> v){
  int low_val = v[0];
  for (int i=1; i<v.size(); i++){
    if (low_val > v[i])
      low_val = v[i];
  }
  return low_val;
}

//get index of lowest value from vector
int low_element_index(vector<int> v){
  int low_index = 0;
  for (int i=1; i<v.size(); i++){
    if (v[low_index] > v[i])
      low_index = i;
  }
  return low_index;
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

  //(m,k)
  vector<double> ms = {10,20,40,80,160,320};
  vector<double> ks = {10,20,40,80,160,320};

  Matrix Ntots_all(6,6);          //matrix to fill with summed total of operation for best coef
  Matrix Ntots_all_coefs(6,6);    //matrix to fill with coef that lead to that sum

  /*
        Trials for coef of (m,k) pair.
  */

  int m, k;
  for (int z=0; z<6; z++){
    for (int y=0; y<6; y++){

      m = ms[z];
      k = ks[y] + m;


      cout << "n:" << m << ", k:" << k << "\n";

      const int NUM_FS = 15;
      vector<double> coefs = Linspace(0.4, 1.0,NUM_FS);
      vector<int> ntots_all_mk_trial;


      for (int j=0; j<NUM_FS; j++){
            //g() interval number update function
        Fupd g;
        g.c = coefs[j];

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
        vector<int> ns;
        vector<int> Ntots;

        for (int i=0; i<NUM_TOL; i++){

          //excepted tolerance
          double rtol = rtols[i];
          double atol = rtol/(1000);
          double tolerance = rtol* fabs(Itrue) + atol;

          // iterate over n values, computing approximations, error, convergence rate
          double Iapprox; //integral estimation
          int Ntot; //number of intervals tried
          int n;    //final number of intervals used


          //
          //    call adaptive solver
          //
          adaptive_int_tester(f, g, m, k, a, b, rtol, atol, Iapprox, n, Ntot);



          //update data
          errors.push_back(fabs(Itrue-Iapprox)/fabs(Itrue));
          Ntots.push_back(Ntot);
          ns.push_back(n);


        }//end tol trial
          //update most resent all ntots for mk piar
          ntots_all_mk_trial.push_back(sum_vector(Ntots));

      }//end coef trial

      //add lowest ntots to (m,k) osition and lowest coef
      Ntots_all.M[z][y] = low_element(ntots_all_mk_trial);
      Ntots_all_coefs.M[z][y] = coefs[low_element_index(ntots_all_mk_trial)];

      printf(" Lowest total number of operations %f,  corresponding coef: %f\n", Ntots_all.M[z][y], Ntots_all_coefs.M[z][y]);

    }//end k
  } //end m

  //write findings to files
  Ntots_all.write("Total_operations.txt");
  Ntots_all_coefs.write("Optimal_coef.txt");

  //write ns, and ksvalues to file
  write(ms, "ns.txt");
  write(ks, "ks.txt");

  //find smallest value in matrix with associated n and k value
  double min_val = Ntots_all.M[0][0];
  int min_index[] = {0,0};
  for (int i=0; i<6; i++){
    for(int j=0; j<6; j++){
       if (min_val > Ntots_all.M[i][j]){
          min_val = Ntots_all.M[i][j];
          min_index[0] = i;
          min_index[1] = j;
        }
       }
  }

  int optimal_n = ms[min_index[0]];
  int optimal_k = ks[min_index[1]] + optimal_n;
  double optimal_coef = Ntots_all_coefs.M[min_index[0]][min_index[1]];

  printf("\n\nOptimal settings:\n\tn:%d, k:%d, coef: %f\n", optimal_n, optimal_k, optimal_coef);
}
