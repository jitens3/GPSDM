/*

-------------------------------------------------------------------------
====== CHANGE LOG ======


*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//using namespace std;

#include "JplGpsDataClass.h"


int calc_k3(JplGpsData data, int i, double &k3, double &dk3){

  int good_epochs=0;
  k3=0;
  //#pragma omp parallel for
  for(int j=0; j<data.num_epochs; j++){
    if(data.ferr[i][j]==0) continue;
    k3 += pow(data.bias[i][j],3);
    good_epochs++;
  }
  k3 /= (good_epochs*pow(data.sdev[i],3));
  dk3 = sqrt(6./good_epochs);

  return 0;
}


int calc_k4(JplGpsData data, int i, double &k4, double &dk4){

  int good_epochs=0;
  k4=0;
  //#pragma omp parallel for
  for(int j=0; j<data.num_epochs; j++){
    if(data.ferr[i][j]==0) continue;
    k4 += pow(data.bias[i][j],4);
    good_epochs++;
  }
  k4 /= (good_epochs*pow(data.sdev[i],4));
  k4 -= 3; //subtract the "Normal" 3
  dk4 = sqrt(24./good_epochs);

  return 0;
}


//******************************************************************************
//******************************************************************************
int main (void)
/*
*/
{

  //Input Parameters:
  std::string path,psd_path,psd_label;
  int minweek,maxweek;
  std::string ref_clock;
  int i_bin,i_k3,isim,i_dif,i_poly;
  std::string label;

  std::string junk;
  std::ifstream fInput;
  fInput.open ("skewness.dat");
    fInput >> i_bin >> path;                     getline(fInput,junk);
    fInput >> minweek >> maxweek;                getline(fInput,junk);
    fInput >> i_k3 >> ref_clock;                 getline(fInput,junk);
    fInput >> i_dif >> i_poly;                   getline(fInput,junk);
    fInput >> isim >> psd_path >> psd_label;     getline(fInput,junk);
    fInput >> label;                             getline(fInput,junk);
  fInput.close();

  bool use_binary = false;
  if(i_bin==1) use_binary=true;

  bool remove_poly = true;
  if(i_poly==0) remove_poly = false;

  std::string file_pref="k3";

  bool do_skew=true;
  if(i_k3 == 1){
    do_skew=false;
    file_pref="k4";
  }

  // Formats input, & reduces liklihood of errors:

  if(minweek<1060||minweek>2000) minweek=1060;
  if(maxweek>2000||maxweek<1060) maxweek=2000;

  bool simulate = false;
  if(isim==1) simulate = true;

  std::string sim_ref_clock = ref_clock; //which reference to simulate
  if(ref_clock=="all" || ref_clock=="any"){
    ref_clock="all";
    sim_ref_clock="USN3";
  }

  std::cout<<"\n          ###########################################"
           <<"\n          Looking for asymmetry/skew + An. modulation"
           <<"\n          ###########################################";
  std::cout<<"\n\n";

  //Output params to screen:
  printf("Searching the 30s clock data, for weeks %i -> %i\n",minweek,maxweek);

////////////////////////////////////////////////////////////////////////////////
// MAIN LOOP:
  time_t start,end; //mid,
  time (&start);

  std::vector< std::vector< std::vector<double> > > skew;
  std::vector< std::vector< std::vector<double> > > kurt;
  std::vector<std::string> clock_list;//hold the clock names "Cs-65" or "H-AMC2"

  // Which clocks to simulate (only used when isim==1)
  std::vector<int> numClks = {4,1,5,1,2,2,1,1};
  // Order for clocks: Rb: F,A,R,II, Cs:F,A,II, white


  std::vector<std::string> type_list = {"Rb II ","Rb IIA","Rb IIR","Rb IIF",
                                        "Cs II ","Cs IIA",         "Cs IIF",
                                        "H R ","Rb R","Cs R"};
  std::vector<std::string> type_name = {"RbII","RbIIA","RbIIR","RbIIF",
                                        "CsII","CsIIA","CsIIF",
                                        "H-R","Rb-R","Cs-R"};


  int minday=0,maxday=7;
  for(int iw=minweek; iw<=maxweek; iw++){ //loops through weeks
    MSC_progressBar("Calculateing",35,iw-minweek,maxweek-minweek+1);
    for(int id=minday; id<maxday; id++){ //loops through each day for given week

      int days = iw*7 + id; //days since JPL week 0. (1980-01-06)

      //Read in the data:
      JplGpsData data;
      int jplok=0;
      if(simulate){
        jplok = data.gpsSimulator(psd_path,psd_label,numClks,0.05,0,0,
          sim_ref_clock);
      }else{
        if(use_binary) jplok = data.binaryJpl30s(path,iw,id);
        else           jplok = data.readJplBiasData(path,iw,id);
      }

      if(jplok==2) std::cout<<"\nXXX FAILURE XXX 126\n";
      if(jplok==1) continue; //this file doesn't exist

      if(data.refprn!=ref_clock && ref_clock!="all") continue;

      if(remove_poly) data.polynomialDetrend(i_poly);
      data.differenceData(i_dif);
      data.calculateStdDev();

      for(int i=0; i<data.num_clocks; i++){
        if(data.clk[i]!="H" && data.clk[i]!="Cs" && data.clk[i]!="Rb"
            && data.clk[i]!="Wh") continue;
        if(data.sdev[i]>0.2) continue; //skip problem clocks
        if(data.sdev[i]<=0) continue; //skip problem clocks //?? DAV

        double k3=0,dk3=0;
        if(do_skew) calc_k3(data, i, k3, dk3);
        else        calc_k4(data, i, k3, dk3);
        std::vector<double> k3_lst = {(double)days,k3,dk3,data.sdev[i]};

        std::string clk_name = data.clk[i]+" ";
        if(data.blk[i]=="AR") clk_name += "R "+data.prn[i];
        else clk_name += data.blk[i]+" "+data.svn[i];

        //Store the skew in the correct list (for given SVN)
        bool new_clock=true;
        int isvn=-1;
        for(size_t k=0; k<clock_list.size(); k++){
          if(clock_list[k]==clk_name){
            new_clock=false;
            isvn = k;
            break;
          }
        }
        if(new_clock){
          clock_list.push_back(clk_name);
          std::vector< std::vector<double> > tmpk3;
          tmpk3.push_back(k3_lst);
          skew.push_back(tmpk3);
        }else{
          skew[isvn].push_back(k3_lst);
        }

      }//END clock


    }//END day
  }//END week
  MSC_progressBar("Calculateing",35);

  //Vectors to hold the averaged (over days) skew
  std::vector<double> av_k3;
  std::vector<double> av_sd;
  std::vector<double> av_dk3; //error, from Sqrt(6/N) formula
  std::vector<double> k3_sd;  //Standard deviation (day-by-day)

  // Find the average skew (for each SVN)
  for(size_t i=0; i<skew.size(); i++){// loop through each clock (SVN)
    double t_k3=0;
    double t_sd=0;
    double t_dk3=0;
    for(size_t j=0; j<skew[i].size(); j++){ //over each day
      t_sd  += skew[i][j][3];
      t_k3  += skew[i][j][1];
      t_dk3 += pow(skew[i][j][2],2); //quadrature? Appropriate?
    }
    int clock_count = skew[i].size();
    av_sd.push_back(t_sd/clock_count);
    av_k3.push_back(t_k3/clock_count);
    av_dk3.push_back(sqrt(t_dk3)/clock_count);
  }

  //calculate standard deviation in skewness (day-by-day)
  for(size_t i=0; i<skew.size(); i++){
    double t_k3_sd=0;
    for(size_t j=0; j<skew[i].size(); j++){
      t_k3_sd += pow(skew[i][j][1]-av_k3[i],2);
    }
    t_k3_sd = sqrt(t_k3_sd/skew[i].size());
    k3_sd.push_back(t_k3_sd);
  }


  //write out daily skew for each clock
  for(size_t n=0; n<type_list.size(); n++){ //each 'block'

    double wav_skew=0, wav_w=0, wav_dk3=0, wav_k3_sd=0, wav_sd=0;
    int wav_clkdays=0;

    std::ofstream ofile,o2file;
    std::string ofname= file_pref+"_day-"+type_name[n]+"-"+ref_clock+"-"+label
      +".out";
    std::string of2name=file_pref+"_avg-"+type_name[n]+"-"+ref_clock+"-"+label
      +".out";
    ofile.open(ofname.c_str());
    o2file.open(of2name.c_str());
    ofile<< "# SVN:    day     k3  dk3  sd(d1) \n";
    ofile<< "#  av:    day     k3  dk3  N-days(only for av)\n";
    //ofile<< "# Average: N_days  k3  dk3  sd(k3)  sd(d1)\n";
    o2file<<"# Clock    N_days  k3  dk3  sd(k3)  sd(d1)\n";

    ofile<<"\""<<type_list[n]<<" av\"\n";
    for(int idate=minweek*7; idate<(maxweek+1)*7; idate++){//loop over each date
      double k3db=0;
      double dk3db=0;
      //double k3db_sd=0;
      int N_db=0;
      for(size_t isvn=0; isvn<skew.size(); isvn++){//loop over svns
        for(size_t j=0; j<skew[isvn].size(); j++){//loop over days
          //pick out correct block:
          if(clock_list[isvn].substr(0,6)!=type_list[n]
            && clock_list[isvn].substr(0,4)!=type_list[n]) continue;
          if(skew[isvn][j][0]>(double)idate) break;
          if(skew[isvn][j][0]!=idate) continue; //pick correct date

          k3db += skew[isvn][j][1];
          dk3db += pow(skew[isvn][j][2],2);
          //k3db_sd += pow(skew[isvn][j][1]-av_k3[isvn],2); //??
          N_db ++;

        }
      }
      //loop over blocks, output
      if(N_db!=0){
        ofile<<idate;
        ofile<<" "<<k3db/N_db<<" "<<sqrt(dk3db)/N_db
                <<" "<<N_db;
        ofile<<"\n";
      }
    }
    ofile<<"\n\n";


    for(size_t i=0; i<skew.size(); i++){
      //pick out correct block:
      if(clock_list[i].substr(0,6)!=type_list[n]
        && clock_list[i].substr(0,4)!=type_list[n]) continue;

      wav_skew += av_k3[i]*skew[i].size();
      wav_w    += skew[i].size();
      wav_dk3  += pow(av_dk3[i]*skew[i].size(),2);
      wav_clkdays += skew[i].size();

      wav_k3_sd += skew[i].size()*pow(k3_sd[i],2);
      wav_sd += skew[i].size()*pow(av_sd[i],2);


      ofile<<"\""<<clock_list[i]<<"\"\n";
      // ofile<<"av: "<<skew[i].size()<<" "<<av_k3[i]<<" "
      //     <<av_dk3[i]<<" "<<k3_sd[i]<<" "<<av_sd[i]<<"\n";
      o2file<<clock_list[i]<<" "<<skew[i].size()<<" "<<av_k3[i]<<" "
              <<av_dk3[i]<<" "<<k3_sd[i]<<" "<<av_sd[i]<<"\n";
      for(size_t j=0; j<skew[i].size(); j++){
        for(size_t k=0; k<skew[i][j].size(); k++){
          ofile<<skew[i][j][k]<<" ";
        }
        ofile<<"\n";
      }
      ofile<<"\n\n";
    }

    wav_skew/=wav_w;
    wav_dk3 = sqrt(wav_dk3)/wav_w;
    wav_k3_sd = sqrt(wav_k3_sd/wav_clkdays);
    wav_sd = sqrt(wav_sd/wav_clkdays);

    o2file<<type_name[n]<<" av "<<wav_clkdays<<" "<<wav_skew<<" "<<wav_dk3<<" "
    <<wav_k3_sd<<" "<<wav_sd<<"\n";

    std::cout<<type_name[n]<<" "<<wav_skew<<" "<<wav_dk3<<"\n";
    ofile.close();
    o2file.close();
  }


////////////////////////////////////////////////////////////////////////////////

  time (&end);
  double dif = difftime (end,start);
  printf ("Elasped time is %.0f seconds.\n", dif );
  return 0;
}
//******************************************************************************
//******************************************************************************
