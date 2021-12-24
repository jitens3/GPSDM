#ifndef _JPLGPSDATA_H
#define _JPLGPSDATA_H
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
//#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm> //used for "sort"
#include <omp.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "mathematicsFunctions.h"
#include "miscFunctions.h"

  //some constants. I often want these outside of the class object!
  const double RGPS=26600; //GPS radius!, in km
  const double REARTH=6371; //earth radius!, in km

//NOTE: Bias data _must_ be doubles, due to huge 'offset' in raw bias data
typedef double bias_data_t;

//******************************************************************************
class JplGpsData {


  //----------------------------------------------------------------------------
  public:

    //Constructor: (is called automatically upon creation)
    JplGpsData(std::string path="", int in_jplweek=0, int in_jplday=0,
                       bool detrend=false, bool difference=false);

    //******************//
    // Public variables //
    //******************//

    //Clock bias and formal error data:
    std::vector< std::vector<bias_data_t> > bias;
    std::vector< std::vector<float> > ferr;

    //Clock and reference positions arrays:
    std::vector< std::vector< std::vector<float> > > pos;
    std::vector< std::vector<float> > refpos;

    //r^2 value for polynomial fit (ever actually needed?)
    std::vector<float> dR2;

    //standard deviation for each clock:
    std::vector<float> sdev;

    // Array to hold the ACF for each clock.
    std::vector< std::vector<float> > acf;

    //Stores the cross-clock correlation.
    //Note: same for each clock pair!
    // Array, in case need non-zero lag. typically, just use 0th element.
    std::vector<float> b0;

    // stores the average network CCC value 03/11/20
    std::vector<double> CCC;
    
    // stores each network end date 03/11/20
    std::vector<double> netEnd;

    // k_effective: "relative" coupling strengths between clocks:
    //nb: reference clock value held at the very end!
    std::vector<float> keff;

    //Number of clocks of each type
    int num_receivers;  //base stations
    int num_satellites; //num_satellites
    int num_clocks;     //clocks? num_clocks?

    //Number of epochs:
    int num_epochs; //epochs

    //If verbose=false, won't print (except some important msgs).
    bool verbose; //False by dflt

    //Stores level of differening.
    //Vector, that isn't re-sized untill diff is called.
    //Therefore, vdiff.size()=0 until it's called!
    std::vector<int> vdiff;

    //Clock name arrays:
    std::vector<std::string> prn;
    std::vector<std::string> svn;
    std::vector<std::string> clk;
    std::vector<std::string> blk;
    //reference clock name:
    std::string refprn;
    std::string refsvn;
    std::string refclk;
    std::string refblk;
    std::string refname;
    
    std::vector<int> tmpsvn;        // added on 09/30/20 tday for Cs elim
    std::vector<std::string> tmpsvnord;        // added on 09/30/20 tday for Cs elim

    //location of JPL files
    std::string jpl_file_dir;

    //dates. For now, string. Make int versions?
    std::string date; //human-readable date ("yyyy-mm-dd")
    std::string week; //JPL week
    std::string day;  //JPL day of week (0=sunday)

    //******************//
    // Public functions //
    //******************//

    //XXX new: reads in 1s data. NOT TESTED
    int readJplBiasData_1s(std::string path, int in_day=0);

    //Read in the JPL clock data files:
    int readJplBiasData(std::string path, int in_jplweek, int in_jplday,
                        bool download=false);

    //Checks to see if JPL files exists. If not, attemps to download them
    int checkJPLfiles(std::string path, bool download=false);

    //reads/writes data+positions to binary file
    int binaryJpl30s(std::string path, int in_jplweek, int in_jplday,
      bool pos=false, bool write=false);

    //Read in the ECI position files (Geoff's files)
    int readEciPos(std::string which="both");

    //Difference the data (and "zero" it)
    //(by default, will just "center/zero" the data)
    int differenceData(int iDif=0);

    //Remove a polynomial
    int polynomialDetrend(int poly_order=2, int weighted=1);

    //Subtract mean:
    int subtractWeightedMean(int dif); //??

    //Calculates the s.d. of each clock:
    double calculateStdDev(void);

    //Read in the ACF functions, form ACF array for each clock
    int formAutoCorrelation(std::string path_to_acf="calculateACF",
                            std::string acf_label="");

    //Calculates the cross-clock correlations.
    //Stores in b0 vector.
    int calculateCrossClockCorrelation(int max_lag=0);

    // replaces the daily CCC value to the network averaged CCC
    // 03/11/20
    int avgNetCCC( int day, int week, std::string path_to_avgccc);

    // reads in the network averaged CCC values
    // 03/11/20
    int avgNetCCCread( int num_networks, std::string path_to_avgccc);

    // Fills array (keff) with the relative coupling strengths.
    int formRelativeCouplings(std::string interaction="");

    //Swaps/changes the reference clock to one of the GPS satellites
    bool swapReference(std::string clock="Rb", std::string block="all",
                       bool by_svn=false);

    int gpsSimulator(std::string pathToPSD, std::string inLabel,
                     std::vector<int> &numClks, double sig=0.05,
                     int in_receivers=0, double sta_sd=0.007,
                     std::string sRef="USN3", int in_epochs=2880,
                     bool useRandSVN=true,
                     bool swap_white_ref=false, double ref_sig=0.0075,
                     bool print=false);

    //A function that copies ALL public members in to a new JPlGpsData object
    int makeACopyOf(JplGpsData inObject);

    //Read in the PRN_GPS file, and map the PRN and SVN etc values:
    int mapSVNPRN(std::string path, bool gps_dm=true);
    int mapStations(std::string path);

  //----------------------------------------------------------------------------
  private:

    //*******************//
    // Private Variables //
    //*******************//

    //Position array indices. Can be modified to include velocities
    const int iXYZ=3;     //3 for pos only, 6 for pos+velocity!??
    const int XX=0, YY=1, ZZ=2;

    static const int JPLEPOCHS=2880; //num epochs in JPL files
    const int MAXPSD=1440; //Maximum length of (half) PSD

    //Number of points of the ACF to take. Typical 5. Set on construction.
    int acf_points;

    //*******************//
    // Private Functions //
    //*******************//

    //XXX this is ineficient! should not have to do this each time!!
    int readSVNPRN(std::string path,
                   std::vector< std::vector<std::string> > &sPRNmap,
                   bool gps_dm);
    //XXX this is ineficient! should not have to do this each time!!
    int readStaMap(std::string path,
                   std::vector< std::vector<std::string> > &sPRNmap);

    //Read in the ECI SATellite and STAtion positions:
    int readEciSatPos(std::string path);
    int readEciStaPos(std::string path);

    // XXX here: binary read/write
    int binaryReadWriteData(std::string fname, bool write=false);
    int binaryReadWritePositions(std::string fname, bool write=false);

    //Reads in a single ACF file
    //int readAcfFile(std::string acf_fname, int ic, int id);   // commented out 09/30/20 by tday
    int readAcfFile(std::string acf_fname, int ic, int id, int nt);
    //Calculates ACF from the data for given day/clock
    //int calculateAcfFromData(int ic); 	// commented out by tday 03/19/21 to update for gpsdays before 12800
    int calculateAcfFromData(int ic, int nt);

    //Functions to re-size the std::vectors. works for 1,2,3D vectors.
    template <class T>
    int reSizeVec(std::vector< std::vector<T> > &invec, int rows, int cols);
    //--overloaded:
    template <class T>
    int reSizeVec(std::vector<T> &invec, int rows);
    //--overloaded:
    template <class T>
    int reSizeVec(std::vector< std::vector< std::vector<T> > > &invec,
                  int rows, int cols, int slcs);

    //GPS simulator:
    int simulatedData(std::vector<double> power_spectrum,
                      std::vector<double> &out_z0, int temp_epochs,
                      int dif_order,
                      gsl_fft_real_wavetable * real_wt,
                      gsl_fft_halfcomplex_wavetable * hc_wt,
                      gsl_fft_real_workspace * work
                      );
    int whiteData(std::vector<double> &outZ0, double sig);
    int simulatedPositions(bool bFull_Random=false);


};
//******************************************************************************











#endif
