/*

-------------------------------------------------------------------------
====== CHANGE LOG ======


*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#include "JplGpsDataClass.h"
#include <gsl/gsl_fft_real.h>


//******************************************************************************
int formPSD(std::vector<double> &d)
/*
This functions takes in a data vector, and calculates the power-spectrum
NOTE: For the factor, assumes sampling rate is 1s.
*/
{


  const int MAXDATASIZE = 2048;
  double data[MAXDATASIZE];

  size_t n = d.size();

  if(n>MAXDATASIZE){
    std::cout<<"\nWARNING 21 in formPSD: Data too long. Using first "
    <<MAXDATASIZE<<" points only!!\n";
    n=MAXDATASIZE;
  }

  //Transfer data from input vector into array
  for(size_t i=0; i<n; i++)
    data[i]=d[i];

  size_t stride=1;
  gsl_fft_real_radix2_transform(data,stride,n);

  // The output is a half-complex sequence, which is stored in-place.
  // The arrangement of the half-complex terms uses the following scheme:
  // for k < n/2 the real part of the k-th term is stored in location k,
  // and the corresponding imaginary part is stored in location n-k.
  // The terms for k=0 and k=n/2 are both purely real,
  // and count as a special case.
  // Their real parts are stored in locations  0 and n/2 respectively,
  // while their imaginary parts which are zero are not stored.

  double t0 = 1.; //Assume sample-time is 1 unit.
  double A = t0/n;

  d.clear();
  d.push_back(A*pow(data[0],2)); //DC (k=0) term [pure real]
  for(size_t k=1; k<n/2; k++){
    double re = data[k];  //real part of FT
    double im = data[n-k];    //real part of FT
    double psd_k = A*(pow(re,2)+pow(im,2));
    d.push_back(psd_k);
  }
  d.push_back(A*pow(data[n/2],2)); //Niquist term [pure real]

  //std::cout<<n<<" "<<n/2+1<<" "<<d.size()<<"\n";


  return 0;

}




//******************************************************************************
//******************************************************************************
int main (void)
/*
*/
{

  int T=128; //XXX MUST be a power of 2!

  //Acceptable fraction of missing points:
  double acc_missing = 0;

  std::string prefix = "../gps1sdata/unr";
  std::string suffix = "_completeDataSet.txt";

  std::vector<int> prev_days = {18621,19622,19623};
  std::vector<int> post_days = {19624,19625};

  //This will hold a list of clock names (PRN-SVN)
  std::vector<std::string> clock_list;

  //XXX OUTPUT:
  // a) Each clock sepperately
  // b) Straight average
  // c) weighted average
  //NOTE: Already weightd, because of P/P_av??


  //Step 1: Just get a list of PRN/SVNs we need
  for(size_t i=0; i<post_days.size(); i++){
    std::string sday = std::to_string(post_days[i]);
    JplGpsData data;
    data.readJplBiasData_1s(prefix+sday+suffix);
    for(int c=0; c<data.num_clocks; c++){
      //form unique clock "name"
      std::string name = data.svn[c]+data.prn[c]+data.clk[c];
      bool newclock = true;
      for(size_t j=0; j<clock_list.size(); j++){
        if(name==clock_list[j]){
          newclock = false;
          break;
        }
      }
      if(newclock) clock_list.push_back(name);
    }
  }
  int total_clocks = clock_list.size();

  //Create list of frequencies (these are the 'x' axis of Fourier transform):
  //There are (T/2+1) frequencies, including DC term
  //Assumes 1s sampling rate!
  std::vector<double> flist;
  int kmax = T/2+1;
  for(int k=0; k<kmax; k++) flist.push_back(double(k)/T);

  //Array to store the averaged PSD for each clock
  std::vector< std::vector<double> > Pav(total_clocks,std::vector<double>(kmax));
  std::vector<int> Ngood(total_clocks);

  int idif = 1;

  //Loop over all non-event (preceeding) days, to avg P
  for(size_t id=0; id<prev_days.size(); id++){
    std::string sday = std::to_string(prev_days[id]);
    std::cout<<sday<<":\n";
    JplGpsData data;
    data.readJplBiasData_1s(prefix+sday+suffix);
    data.differenceData(idif);
    int tmin = idif;
    int tmax = data.bias[0].size()-T-1;

    for(int i=0; i<data.num_clocks; i++){

      //Identify clock index:
      int ic=-1;
      std::string this_name = data.svn[i]+data.prn[i]+data.clk[i];
      bool found = false;
      for(size_t c=0; c<clock_list.size(); c++){
        if(this_name == clock_list[c]){
          found=true;
          ic = c;
          break;
        }
      }
      if(!found) continue;

      for(int j=tmin; j<tmax; j++){
        //Count number of 'bad epochs'. if too many, skip
        int bad_epochs=0;
        for(int t=j; t<j+T; t++){
          if(data.ferr[i][j]==0) bad_epochs++;
        }
        if (double(bad_epochs)/T > acc_missing) continue;
        //Do FT/PSD for chunk of data (moving window)
        std::vector<double> tempd(data.bias[i].begin()+j, data.bias[i].begin()+j+T);
        formPSD(tempd);
        for(size_t j2=0; j2<tempd.size(); j2++) Pav[ic][j2] += tempd[j2];
        Ngood[ic]++;
      }
    }
  }

  //Devide by total number of "good" windows - to get the average!
  for(size_t c=0; c<Pav.size(); c++){
    for(size_t j2=0; j2<Pav[c].size(); j2++){
      Pav[c][j2] /= Ngood[c];
    }
  }


  //Loop over all post-event days
  for(size_t id=0; id<post_days.size(); id++){
    std::string sday = std::to_string(post_days[id]);
    std::cout<<sday<<":\n";
    JplGpsData data;
    data.readJplBiasData_1s(prefix+sday+suffix);
    data.differenceData(idif);
    // int tmin = idif;
    // int tmax = data.bias[0].size()-T-1;

    std::cout<<data.bias[0].size()<<"\n";

    std::ofstream of;
    of.open("P-"+sday+"_"+std::to_string(T)+".out");

    for(int i=0; i<data.num_clocks; i++){

      //Identify clock index:
      int ic=-1;
      std::string this_name = data.svn[i]+data.prn[i]+data.clk[i];
      bool found = false;
      for(size_t c=0; c<clock_list.size(); c++){
        if(this_name == clock_list[c]){
          found=true;
          ic = c;
          break;
        }
      }
      if(!found) continue;

      of<<this_name<<" \n"; //for header title

      for(int j=idif; j<data.num_epochs-T; j++){
        std::vector<double> tempd(data.bias[i].begin()+j, data.bias[i].begin()+j+T);
        formPSD(tempd);
        int t = j+T/2;
        for(size_t k=0; k<tempd.size(); k++){
          of<<t<<" "<<flist[k]<<" "<<(tempd[k]/Pav[ic][k])<<"\n";
        }
      }
      of<<"\n";
    }
    of.close();
  }



  return 1;

  // JplGpsData data;
  // data.readJplBiasData_1s("../gps1sdata/unr19624_completeDataSet.txt");
  //
  // data.differenceData(1);
  // data.calculateStdDev();
  //
  // //Length of data to FT for each time.
  //
  //
  // //Create list of frequencies (these are the 'x' axis of Fourier transform):
  // //There are (T/2+1) frequencies, including DC term
  // //Assumes 1s sampling rate!
  // // std::vector<double> flist;
  // // for(int k=0; k<T/2+1; k++)
  // //   flist.push_back(double(k)/T);
  //
  // //for now, don't loop over clocks. Just choose one.
  // int i=20;
  //
  // //This part just for testing
  // /*
  // //(injects sinusoidal signal that peaks at t=125, and lasts for ~75s)
  // for(int j=0; j<data.num_epochs; j++)
  // {
  //   double ti = 4000;
  //   double st = 500;
  //   double fj = 0.18;
  //   if(j>ti) fj = 0.18*(1.+0.7*(j-ti)/st);
  //   data.bias[i][j] += 2*data.sdev[i]*sin(2.*M_PI*fj*j)*exp(-0.5*pow(ti-j,2)/pow(st,2));
  // }
  // */
  //
  // //Calculate PSD, and write to file (function of t and k)
  // std::ofstream of;
  // of.open("psdx.out");
  // //for(int i=0; i<data.num_clocks; i++)
  // {
  //   //of<<"# "<<data.prn[i]<<"-"<<data.svn[i]<<" "<<data.sdev[i]<<"\n";
  //   for(int j=2; j<data.num_epochs-T; j++)
  //   //for(int j=2; j<500; j++) //to use fewer points (speed)
  //   {
  //     int t = j+T/2;
  //     //create temporary small vector
  //     //NOTE: there are more efficient ways of doint this..but probs fine.
  //     std::vector<double> tempd(data.bias[i].begin()+j, data.bias[i].begin()+j+T);
  //     formPSD(tempd);
  //     for(size_t k=0; k<tempd.size(); k++){
  //       of<<t<<" "<<flist[k]<<" "<<log10(tempd[k])<<"\n";
  //     }
  //   }
  //   //of<<"\n";
  // }
  //
  // of.close();
  //
  //
  //
  // return 1;



}
//******************************************************************************
//******************************************************************************
