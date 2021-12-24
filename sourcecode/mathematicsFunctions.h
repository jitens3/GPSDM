#ifndef _MATHFUNS_H
#define _MATHFUNS_H
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
// #include <gsl/gsl_fft_real.h>
#include <random>
#include <thread>
// #include <gsl/gsl_rng.h>

// //XXX ???
// //This has got to do with the parallelisable random number generator
// // (see MFS_randDouble)
// //I have tested with g++. I think _MSC_VER is windows??
// #if defined (_MSC_VER)  // Visual studio
//     #define thread_local __declspec( thread )
// #elif defined (__GCC__) // GCC
//     #define thread_local __thread
// #endif

const double PI = M_PI; //?

///////////////////////////////////////////////////////////////////////////////


bool MFS_polywfit(int obs,int degree, std::vector<double> dx,
              std::vector<double> dy, std::vector<double> dw,
              std::vector<double> &coef, double &chisq, double &R2);


double MFS_randGausVal(double x0, double sig);
double MFS_randDouble(double a, double b);


int MFS_invertMatrix(std::vector< std::vector<double> > inmat,
                        std::vector< std::vector<double> > &outmat);
int MFS_invertMatrix(std::vector< std::vector<double> > &inmat);
int MFS_invertMatrix(std::vector< std::vector<float> > inmat,
                        std::vector< std::vector<float> > &outmat);
int MFS_invertMatrix(std::vector< std::vector<float> > &inmat);

double MFS_calcDeterminant(std::vector< std::vector<double> > inmat);
double MFS_calcDeterminant(std::vector< std::vector<float> > inmat);

int MFS_linsolve(std::vector< std::vector<double> > inmat,
                 std::vector<double> invec, std::vector<double> &outvec);


double MFS_inverseErf(double x);
double MFS_fastInvErf(double x);
double MFS_fastErf(double x);
double MFS_fastExp(double x);



#endif
