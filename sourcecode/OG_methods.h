#ifndef _FEVNT_H
#define _FEVNT_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <algorithm> //for sort
#include <time.h>
#include <omp.h>

#include "JplGpsDataClass.h"
#include "mathematicsFunctions.h"
#include "miscFunctions.h"
#include "NumericCdfInverseClass.h"

const int MAXITS=2048; //250 //??? XXX
const int MAXJW=61;
const int MAXCLOCKS=65; //??
const int MAXEPOCHS=2880;

int MClikelihoods(
  std::string model,
  JplGpsData data,
  int iCov, bool use_ccc, bool sub_w_mean,
  std::vector< std::vector<double> > &odds,
  int jbeg, int jend, int JW,
  int NmMC,
  double dh,
  NumericCdfInverse v_prior, 
  NumericCdfInverse psi_prior,
  double *Einv,
  int exact
);

int MCsnr(std::string model, JplGpsData data, int icov, 
	bool use_ccc, bool sub_w_mean, double *snr,
	int jbeg, int jend, int jw,	int NMC, double h0,
	NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  int exact, double * Einv
); 


int likeParameterEstimation(
  std::string model,
  JplGpsData data,
  int iCov, bool use_ccc, bool sub_w_mean,
  int itn_or_day,
  int j0, int JW,
  int NMC,
  double h0,
  NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  std::vector<double> &best_params,
  double *Einv,
  int exact
);

int snrParameterEstimation(
  std::string model,
  JplGpsData data,
  int icov, bool use_ccc, bool sub_w_mean,
  int itn_or_day,
  int j0, int JW,
  int NMC,
  double h0,
  NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  int exact, double * Einv, int invest
);

double likelihoodLogNorm(int iCov, int Nclk, int JW, std::vector<float> csd,
                         double * Hijl);

void prepareData(int hJW, int j0, int &jmin,
                 double dD[MAXCLOCKS][MAXJW],
                 double rsat[MAXCLOCKS][3], double rref[3],
                 JplGpsData data);

void prepareDataSearch(int hJW, int j0, int &jmin,
                 double dD[MAXCLOCKS][MAXJW],
                 double rsat[MAXCLOCKS][3], double rref[3],
                 JplGpsData data);

void randomParameters(NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
                        std::string model, int j0,
                        double &t0, double &v, double n[3],
                        double &d, double &R, double &a);

double hAnalytic(double ds, double ss, double h0);

int weightedMeanSignal(double s[MAXCLOCKS][MAXJW], std::vector<float> sdev,
                       double W, int num_clocks, int JW);

int genHijl(
  JplGpsData data,
  double *hessian, //output Hessian matrix
  int JW,
  int max_cov,
  double acf_cut=0 //ACF/Hessian cut-off threshold
  );

int genHijlCholesky(
  JplGpsData data,
  double Hijl[MAXCLOCKS][MAXJW][MAXJW], //output Hessian matrix
  int JW,
  int max_cov,
  double acf_cut=0 //ACF/Hessian cut-off threshold
  );

int calcEinv( double * Einv, JplGpsData & data, int jw, int icov,
      double max_cov, double acf_cut, int day, int week, double ref_sig, 
      int by_hand_ccc, int sanity);


int calcEinvWhiteNoise( double * Einv, JplGpsData & data, int jw, 
    double sigma, double sigma_x);

double calcChiExact(double x[MAXCLOCKS][MAXJW], 
  double y[MAXCLOCKS][MAXJW], double * Einv, 
  int numClks, int jw);

double calcChi_BLAS(double *E_inverse, double *X, 
  double *Y, int num_clocks, int jw);

double calcChi(
  double d1[MAXCLOCKS][MAXJW], //input "first" vector (d or s)
  double d2[MAXCLOCKS][MAXJW], //input "second" vector (d or s)
  double *hessian, //Input Hessian matrix
  //double csd[MAXCLOCKS], //input clock standard deviation
  std::vector<float> csd, //input clock standard deviation
  //int jmin, int jmax,
  int JW,
  int iClocks, int iCov);


double calcChiW(
  double d1[MAXCLOCKS][MAXJW], //input "first" vector (d or s)
  double d2[MAXCLOCKS][MAXJW], //input "second" vector (d or s)
  double b0, //Input cross-clock correlation term
  //double csd[MAXCLOCKS], //input clock standard deviation
  std::vector<float> csd, //input clock standard deviation
  int JW,
  int iClocks);

// this is a new function 11/26/19 that finds the vector v
double calcdE(
              double d1[MAXCLOCKS][MAXJW],       // this is the data vector
              double *m_Einv_e,     // this is the hessian matrix
              double dEij[MAXCLOCKS*MAXJW],         // this is the saved vector variable
              int JW,
              int num_clocks);

// this is a new function 11/26/19 that finds the chi_ds term for exact case
double calcdEs(
              double dEij[MAXCLOCKS*MAXJW],       // this is the saved vector
              double d2[MAXCLOCKS][MAXJW],         // this is the saved vector variable
              int JW,
               int num_clocks);


// *** Signal Generator ***
int generateSignal(std::string model, double s[MAXCLOCKS][MAXJW],
                   double rsat[MAXCLOCKS][3], double rref[3],
                   int num_clocks, int jmin, int JW,
                   std::vector<int> dif, std::vector<float> keff,
                   double t0, double v, double eDM[3],
                   double d=0, double R=0, double a=0);
//overloaded:
int generateSignal(std::string model, double s[MAXCLOCKS][MAXJW],
                   double rsat[MAXCLOCKS][3], double rref[3],
                   int num_clocks, int jmin, int JW, int id,
                   std::vector<float> keff,
                   double t0, double v, double eDM[3],
                   double d=0, double R=0, double a=0);


int thinWallSignal(double s[MAXCLOCKS][MAXJW], double rsat[MAXCLOCKS][3],
                   int num_clocks, int jmin, int JW,
                   std::vector<int> dif, std::vector<float> keff,
                   double t0, double v, double eDM[3],
                   double tr, double rel_hr);

double expqd2(std::string model, double rx, double ry, double rz, double xi,
              double eDM[3], double d, double a, double R);
double projectClock(double rx, double ry, double rz, double eDM[3]);
double perpendicularClock(double rx, double ry, double rz, double xi);
double cosGamma(double rx, double ry, double rz, double xi, double eDM[3],
                double a);

double firstOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r);
double secondOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r);
double zerothOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r);

int randJumpSignal(double t0, int iClocks, int jmin, int JW,
                   double s[MAXCLOCKS][MAXJW]);
int refOnlySignal(double t0, double v, double rref[3], double eDM[3],
                  int iClocks, int jmin, int JW, double s[MAXCLOCKS][MAXJW]);
int refOnlySignal(double tr, int iClocks, int jmin, int JW,
                  double s[MAXCLOCKS][MAXJW]);


#endif
