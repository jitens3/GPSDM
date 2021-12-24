#ifndef _NOISEPFUNS_H
#define _NOISEPFUNS_H
#include <cmath>
#include <omp.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

int allanVariance(double d[], int M, double av[], int N, bool para);

int autoCorrelation(double d[], int M, double acf[], int N, int jmin,
                    bool para);

int powerSpectrumFFT(double d[], int M, double psd[], int N, int dif,
                  gsl_fft_real_workspace * work,
                  gsl_fft_real_wavetable * real_wt
                  );

int histogram(double d[], int M, double hist[], int N, double w, int jmin);


#endif
