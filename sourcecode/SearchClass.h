#ifndef _CSEARCH_H
#define _CSEARCH_H
#include <iostream>
#include <fstream>
#include <string>

#include <cstdlib>

#include "JplGpsDataClass.h"
#include "NumericCdfInverseClass.h"
#include "methods.h"


class SearchClass : public JplGpsData{
  //friend std::ofstream & operator<< (std::ofstream & os, SearchClass & s) //what is this and why are we using it.
  
  public:

    //********************************//
    //  Constrcutors and Destructors  //
    //********************************//

    //We need a default constructor/destructor
	SearchClass(std::string path_to_data = "", int NMC = 0, 
				int jplweek = 0, int jplday=0, 
				bool detrend = false, bool difference = 
        false, double o = 0, double s = 0, double h = 0);

  ~SearchClass();
    //******************//
    // Public variables //
    //******************//

    //Any Jpl Data information
    //JplGpsData m_data; we dont need this as this class inherits all elements from JplGpsData
	
    //A storage for the SNR for each epoch, each day.
    double * m_snr;
    double * m_snr_temp;
    double * tday_snr;

   // storage of the velocities for all NMC and for all best matches
    double * m_vel;
    double * vel_tmp;

    //A storage for the odds and log odds for each epoch, each day.
    double * m_odds;

    //A storage for the likelihood for each epoch, each day.
    //each epoch has 2000 elements: a sHs and a sHd for each monte carlo iteration.
    double * m_likelihood; //we need to print max size here every day so 
									 //we know that we can store all the necessary vals
	  double * dHd;
    
    // storage for the chi squared values
    double * m_chi;
    double * m_chi_temp;

    //A storage for the likelihood for each epoch, each day.
    //each epoch has two elements: a mean and a standard dev.
    std::vector<double> m_posterior;

    //A storage for the h values for each epoch, each day.
    //each epoch has two elements: an h* and a standard dev.
    std::vector<double> m_h;

    //A storage of the inverse covariance matrix, need to pass pointer to inverse functions!!!
    double* m_Einv_e; 

    double* m_Einv_p;

    //cut-offs for odds and SNR, hmin for analytic h-integral
    double odds_cut;
    double SNR_cut;
	  double hmin;

    int m_investigate;
    int m_estNMC;

    //******************//
    // Public functions //
    //******************//

	
//actually should be able to reduce some of this since a lot of these are part of JplGpsData
int MCSearchStats(
  std::string model,
  SearchClass &obj,
  int icov, int JW,
  int jbeg, int jend,
  int NMC,
  NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  int exact,
  int week,
  int day,
  int ilimit_cores,
  int sanity,
  int event
);

    // this is for chi_squared filter
    int chi_thresh;

void serialize(int week, int day, int sanity, int NMC, int jbeg, int jend, int chi_thresh);


};

#endif
