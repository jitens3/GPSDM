#include "JplGpsDataClass.h"
#include "noiseProfileFunctions.h"
/*
170904.
New program, replaces older version of same name.

Calculates power spectrums, autocorrelation functions, histograms,
Allan variance, and standard deviation (with uncertianty), for each SVN, clock,
and reference clock combination.
Then, averages those functions, to find the above for each block (averaged over)
SVNs, and averaged over all reference clocks.
Also, it can "swap" the reference clock to one of the other satellite clocks.

Note: The list of possible reference clocks is hard-coded in.
It is easy to update.
Probably, this could be automated.
Also hard-coded in, is which SVNs belong to which satellite block.
When block III satellites are launched, this needs to be updated [svnToBlock()].
If not, the newer satellites will be interpreted as block IIF.

All output files are placed into a directory specified by the user.
The output files are names like:
**ClkBlk-svn-REF-outlabel.ext**.
e.g.,
  * RbIIR-61-USN3-test1.acf

hold the autocorrelation function for the Rb IIR clock with svn=61, for days
when USN3 reference was used

Other information about the run is written to the header lines of each file.
The header lines are marked with an '#'.
Generally, the data is given as a table, each row corresponding to an epoch/lag/
frequency component, and each collumn for 0th, 1st or 2nd order differecing.

The standard deviations are written to the ACF files.
The summary of all standard deviations (including the uncertianty in those
standard deviations) are written to another file:
**sdWunc-outlabel.txt**

  XXX - Sepperate RbIIR and RbIIRM ??

=====================
CHANGE LOG:
170905- Added if(jplok!=0)continue; to fix errors when JPL data not exists.
170912- Added ability to use this code to test the GPS simulator!
170921- Added GSL FFT to the power spectrum calculation. ~17x faster!
170928- Starting to add Cross-Clock Correlations
*/


//******************************************************************************
//******************************************************************************
//function declaration:
std::string svnToBlock(int in_svn);
//******************************************************************************
int main (void){

  time_t start,end;
  time (&start);

  //define input parameters:
  std::string path;     //directory of input files
  int ibin;
  std::string odir;     //directory for output files
  int minweek,maxweek;  //first/last week of data to consider
  int weekday;          //which day to consider? 0=sunday. 7=all
  std::string sref;
  std::string out_label;


  //read in the input file:
  std::string junk;//holds the "junk" after the input variables.
  std::ifstream fInput;
  fInput.open ("stationNoise.dat");
    fInput >> ibin >> path;                           getline(fInput,junk);
    fInput >> minweek >> maxweek;                     getline(fInput,junk);
    fInput >> sref;                                   getline(fInput,junk);
  fInput.close();

  bool binary=false;
  if(ibin==1) binary=true;



  //Formats input, & reduces liklihood of errors
  if(minweek<1060||minweek>2000) minweek=1060;
  if(maxweek>3000) maxweek=2000;
  if(maxweek==0) maxweek=minweek;
  if(maxweek<minweek) maxweek=minweek;
  if(weekday<0||weekday>7) weekday=7;



  // //Make the output directory (if it doesn't exist already)!
  // MSC_execute("mkdir -p "+odir);
  // if(!MSC_direxists(odir)){
  //   std::cout<<"ERROR 59: Directory "<<odir<<" does not exist. "
  //            <<"Please create it, and then try again.\n";
  //   return 1;
  // }


  std::cout<<"\n          ###########################################";
  std::cout<<"\n                    Process Station Noise            ";
  std::cout<<"\n          ###########################################\n\n";

  //----------------------------------------------------------------------------


  //std::string sref = "USN7";

  std::vector<std::string> stations;
  std::vector<int> n_sta;
  std::vector<double> sd;

  double ref_cc=0;
  int n_days=0;


  int init_day=minweek*7;
  int max_day=(maxweek+1)*7;
  for(int days=init_day; days<max_day; days++){
    MSC_progressBar("",40,days-init_day,max_day-init_day);
    int iw=(int)floor(days/7);
    int id=days - iw*7;

    JplGpsData data;
    int jplok;
    if(binary) jplok = data.binaryJpl30s(path,iw,id);
    else       jplok = data.readJplBiasData(path,iw,id);

    if(jplok!=0) continue;

    if(data.refname != sref) continue;

    data.differenceData(1);

    data.calculateStdDev(); //little wasteful..

    double tmp=0;
    double tmp_norm=0;
    int n_pairs=0;

    //Calculate standard deviations
    for(int i=0; i<data.num_receivers; i++){

      if(data.clk[i]!="H" && data.clk[i]!="Rb" && data.clk[i]!="Cs") continue;
      if(data.sdev[i]>1.0) continue;

      std::string sta_name = data.clk[i]+"-"+data.prn[i];

      bool new_sta = true;
      int index=-1;
      for(size_t k=0; k<stations.size(); k++){
        if(sta_name==stations[k]){
          new_sta=false;
          index = k;
          break;
        }
      }

      if(new_sta){
        stations.push_back(sta_name);
        n_sta.push_back(1);
        sd.push_back(data.sdev[i]);
      }else{
        n_sta[index] ++;
        sd[index] += data.sdev[i];
      }

      if(data.clk[i]!="H") continue;
      if(data.sdev[i]>0.5) continue;

      for(int k=i+1; k<data.num_receivers; k++){
        if(data.clk[k]!="H") continue;
        if(data.sdev[k]>0.5) continue;

        double tmp_ik=0;
        int good_epochs=0;
        for(int j=0; j<data.num_epochs; j++){
          if(data.ferr[i][j]*data.ferr[k][j]==0) continue;
          good_epochs++;
          tmp_ik += data.bias[i][j]*data.bias[k][j];
        }
        if(good_epochs!=0){
          tmp += tmp_ik/(good_epochs*pow(data.sdev[i]*data.sdev[k],2));
          tmp_norm += 1./(pow(data.sdev[i]*data.sdev[k],2));
          n_pairs++;
        }
      }


    }

    if(n_pairs!=0){
      //std::cout<<iw<<id<<" "<<tmp/tmp_norm<<" "<<n_pairs<<"\n";
      ref_cc += tmp/tmp_norm;
      n_days++;
    }

  }//END loop over days
  MSC_progressBar("",40);

  for(size_t i=0; i<stations.size(); i++) sd[i] /= n_sta[i];
  ref_cc /= n_days;

  std::cout<<"Station standard deviations (with "<<sref<<" as reference)\n";
  for(size_t i=0; i<stations.size(); i++){
    if(stations[i].substr(0,1)=="C"){
      std::cout<<stations[i]<<" "<<n_sta[i]<<" "<<sd[i]<<"\n";
    }
  }
  std::cout<<"\n";
  for(size_t i=0; i<stations.size(); i++){
    if(stations[i].substr(0,1)=="R"){
      std::cout<<stations[i]<<" "<<n_sta[i]<<" "<<sd[i]<<"\n";
    }
  }
  std::cout<<"\n";
  for(size_t i=0; i<stations.size(); i++){
    if(stations[i].substr(0,1)=="H"){
      std::cout<<stations[i]<<" "<<n_sta[i]<<" "<<sd[i]<<"\n";
    }
  }
  std::cout<<"\n";
  std::cout<<"Reference clock s.d. (calculated from cross-correlations)\n";
  std::cout<<sref<<" "<<n_days<<" "<<sqrt(ref_cc)<<"\n";


  //change to H:mm:ss
  time (&end);
  double dif = difftime (end,start);
  printf ("Total elasped time is %.0f seconds.\n", dif );

  return 0;
}
