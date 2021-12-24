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
  std::string odir;     //directory for output files
  int minweek,maxweek;  //first/last week of data to consider
  int weekday;          //which day to consider? 0=sunday. 7=all
  int iunc, ip, iac, iav, ih, iccc; //which functions to calc.
  double dexclsd; //exclude clocks?
  int iexclfer;
  double dexclO; //exclude based on outliers. 0=no, dbl=nanoseconds
  int ipoly,iweight; //remove poly?
  int iswapref;          //swap reference?
  std::string out_label;
  //////
  int iSim;   //use simulated data?
  std::string psd_label; // PSD label

  //read in the input file:
  std::string junk;//holds the "junk" after the input variables.
  std::ifstream fInput;
  fInput.open ("processNoise.dat");
    fInput >> path;                                   getline(fInput,junk);
    fInput >> odir;                                   getline(fInput,junk);
    fInput >> minweek >> maxweek;                     getline(fInput,junk);
    fInput >> weekday;                                getline(fInput,junk);
    fInput >> iunc >> ip >> iac >> iav >> ih >> iccc; getline(fInput,junk);
    fInput >> dexclsd >> iexclfer >> dexclO;          getline(fInput,junk);
    fInput >> ipoly >> iweight;                       getline(fInput,junk);
    fInput >> iswapref;                               getline(fInput,junk);
    fInput >> out_label;                              getline(fInput,junk);
    fInput >> iSim >> psd_label;                      getline(fInput,junk);
  fInput.close();

  //Formats input, & reduces liklihood of errors
  if(minweek<1060||minweek>2000) minweek=1060;
  if(maxweek>3000) maxweek=2000;
  if(maxweek==0) maxweek=minweek;
  if(maxweek<minweek) maxweek=minweek;
  if(weekday<0||weekday>7) weekday=7;

  bool remove_poly=true;
  if(ipoly<=0)remove_poly=false;
  if(iweight!=1)iweight=0;

  bool single_day=false;
  if(weekday!=7)single_day=true;

  bool day_para=true; //parallelise by day? Seems to work fine! :)

  //format output label:
  if(out_label=="na")out_label="";
  else out_label="-"+out_label;

  //XX also: if_swap_references... (iswp_min)
  bool do_psd=true;
  bool do_acf=true;
  bool do_avr=true;
  bool do_hist=true;
  bool sd_only=false;
  bool do_sdUnc=true;
  bool do_clkcor = true;
  if(ip==0)do_psd=false;
  if(iac==0)do_acf=false;
  if(iav==0)do_avr=false;
  if(ih==0)do_hist=false;
  if(iccc==0)do_clkcor=false;
  if(iunc==0)do_sdUnc=false;

  if(!do_psd&&!do_acf&&!do_avr&&!do_hist&&!do_clkcor)sd_only=true;

  //Exlude clocks?
  bool exclsd=false;
  bool exclfer=false;
  bool exclO=false;
  if(dexclsd>0)exclsd=true;
  if(iexclfer>0)exclfer=true;
  if(dexclO>0)exclO=true;     //exlude outliers (based on 1st-order difference

  //Make the output directory (if it doesn't exist already)!
  MSC_execute("mkdir -p "+odir);
  if(!MSC_direxists(odir)){
    std::cout<<"ERROR 59: Directory "<<odir<<" does not exist. "
             <<"Please create it, and then try again.\n";
    return 1;
  }

  // For running for simulated data (to test GPS simulator)
  bool test_simulated_data=false;
  if(iSim==1){
    test_simulated_data=true;
    day_para=false; // bottle-neck is generating data
  }
  std::string PSDpath=path; //use same line!!
  // Number of each type of clock. Would be better as input, but don't want to
  // clutter up the input file.
  std::vector<int> numClks(8);
  numClks[0]=10;  //RbIIF
  numClks[1]=10;  //RbIIA
  numClks[2]=15;  //RbIIR
  numClks[3]=5;   //RbII
  numClks[4]=10;  //CsIIF
  numClks[5]=10;  //CsIIA
  numClks[6]=5;   //CsII
  numClks[7]=0; // white clocks

  std::cout<<"\n          ###########################################";
  std::cout<<"\n                      Process Clock Noise            ";
  std::cout<<"\n          ###########################################\n\n";

  printf("Processing noise for weeks %i - %i\n\n",minweek,maxweek);
  std::cout<<"Calculating ";
  if(do_acf)std::cout<<"ACF, ";
  if(do_psd)std::cout<<"PSD, ";
  if(do_hist)std::cout<<"Histogram, ";
  if(do_avr)std::cout<<"AVAR, ";
  if(!sd_only)std::cout<<"and ";
  std::cout<<"standard deviations";
  if(do_sdUnc)std::cout<<" (with uncertainties)";
  std::cout<<"\n\n";

  if(dexclsd==0 && iexclfer==0 && dexclO == 0){
    printf("Not excluding any days.\n\n");
  }else{
    printf("Exluding days with:\n");
    if(dexclsd!=0)  printf("  sigma(1) > %.2f ns\n",dexclsd);
    if(iexclfer!=0) printf("  %i or more missing data points\n",iexclfer);
    if(dexclO!=0)   printf("  any outliers > %.1f ns\n",dexclO);
    std::cout<<"\n";
  }

  if(test_simulated_data)
    printf("Testing simulated data!\n\n");

  //----------------------------------------------------------------------------

  const int NUMEPOCHS=2880;
  const int NUMSVNS=70;     //allows up to ""+min_svn=83
  const int NUMREFS=25;     //Num reference clocks!* see below
  const int NUMCLOCKS=2;
  const int NUMDIFF=3;

  const int ACFPOINTS=128; //don't need that many!
  const int AVRPOINTS=1440;
  const int PSDPOINTS=1441; //epochs/2 +1 (+1 for DC component)

  const int NUMBINS=2001;//2001; //nb: must be ODD number!

  double max_s1=1.;
  double bin_width=2.*max_s1/(NUMBINS); //width of each bin

  //* first 7 (0->6) are the Cs/Rb clocks!
  // "other" reserved as last one! (just in case)
  std::string refname[NUMREFS]=
      {"RbII", "RbIIA", "RbIIR", "RbIIF", "CsII", "CsIIA", "CsIIF",
       "USN3", "AMC2", "USN7", "IRKT", "USNO", "PTBB", "NRC1", "CRO1", "YEBE",
       "ONSA", "KOKB", "TWTF", "BOR1", "NYA1", "MIZU", "WTZR", "NIST", "OTHER"};
  int iHref=7; //"initial H-maser reference (or # of non-H refs)

  //Swap reference clocks (or not!):
  int iswp_min=0;
  if(iswapref==0)iswp_min=iHref;

  //The lowest SVN that exists in our data
  // (introduced so the array can start at 0, then each svn index is
  // off-set from the actual svn by the value of min_svn)
  int min_svn=13;

  //These arrays are static because they are very large,
  // and need to be stored in the 'heap' instead of the 'stack'
  static float PSD[NUMDIFF][NUMCLOCKS][NUMREFS][NUMSVNS][PSDPOINTS]={0};
  static float ACF[NUMDIFF][NUMCLOCKS][NUMREFS][NUMSVNS][ACFPOINTS]={0};
  static float AVR[NUMCLOCKS][NUMREFS][NUMSVNS][AVRPOINTS]={0};
  static float SD[2][NUMCLOCKS][NUMREFS][NUMSVNS]={0};
  static double HIST[2][NUMCLOCKS][NUMREFS][NUMSVNS][NUMBINS]={0};
  int clockdays[NUMCLOCKS][NUMREFS][NUMSVNS]={0};
  int excldays[NUMCLOCKS][NUMREFS][NUMSVNS]={0};

  static float CLKCOR[2][NUMCLOCKS][NUMREFS][NUMSVNS]={0};

  //Prepare some things for the GSL
  //Prepare trigonometric lookup tables for an FFT of size n real elements:
  //Note: We only need to form this once (cen be re-used)
  gsl_fft_real_wavetable * real_wt_0;
  gsl_fft_real_wavetable * real_wt_1;
  gsl_fft_real_wavetable * real_wt_2;
  if(do_psd)
  {
    real_wt_0 = gsl_fft_real_wavetable_alloc (NUMEPOCHS);
    real_wt_1 = gsl_fft_real_wavetable_alloc (NUMEPOCHS-1);
    real_wt_2 = gsl_fft_real_wavetable_alloc (NUMEPOCHS-2);
  }

  //*************************************
  //Loop through each day of data. Read it in, calculate PSD,ACF etc.
  int init_day=minweek*7;
  int max_day=(maxweek+1)*7;
  #pragma omp parallel for if(day_para)
  for(int days=init_day; days<max_day; days++){
    MSC_progressBar("Processing Noise:",40,days-init_day,max_day-init_day,
                    day_para);

    //convert to week+day
    int iw=(int)floor(days/7);
    int id=days - iw*7;
    if(single_day&&(id!=weekday))continue;

    //create "master" data object, that will be duplicated
    //(duplication means we can difference 0,1,2 etc)
    JplGpsData master_data;
    if(test_simulated_data){
      //Use GPS simulator to generate data to test!
      //nb: don't remove poly from this!
      std::string s_ref="USN3";
      double dref=MFS_randDouble(0,2);
      if(dref>1)s_ref="AMC2";

      int jplok = master_data.gpsSimulator(PSDpath,psd_label,numClks,0.,0,0,
                                           s_ref,NUMEPOCHS);

      if(jplok!=0)continue; //didn't find data:
    }else{
      //use real data!
      int jplok = master_data.readJplBiasData(path,iw,id);
      if(jplok!=0)continue; //didn't find data:
      //remove polynomial?
      if(remove_poly) master_data.polynomialDetrend(ipoly,iweight);
    }


    //Exclude clocks (or not) based on missed data, SD, or
    std::vector<bool> bexcl(master_data.num_clocks);
    {
      JplGpsData test_data;
      test_data.makeACopyOf(master_data);
      test_data.differenceData(1);
      test_data.calculateStdDev();
      for(int i=master_data.num_receivers; i<master_data.num_clocks; i++){
        bexcl[i]=false;
        if(exclsd && test_data.sdev[i]>dexclsd){
          bexcl[i]=true;
          continue; //go to next clock, no point checking further!
        }
        int num_missed=0;
        for(int j=1; j<test_data.num_epochs; j++){
          if(test_data.ferr[i][j]==0)num_missed++;
          if(num_missed>=iexclfer && exclfer) bexcl[i]=true;
          if(test_data.bias[i][j]>dexclO && exclO) bexcl[i]=true;
          if(bexcl[i]) break;
        }
      }
    }


    //Prepare workspace for GSL Fourier transforms
    gsl_fft_real_workspace * work_0=0;
    gsl_fft_real_workspace * work_1=0;
    gsl_fft_real_workspace * work_2=0;
    if(do_psd){
      work_0 = gsl_fft_real_workspace_alloc (NUMEPOCHS);
      work_1 = gsl_fft_real_workspace_alloc (NUMEPOCHS-1);
      work_2 = gsl_fft_real_workspace_alloc (NUMEPOCHS-2);
    }

    //loop over differecing orders:
    for(int dif=0; dif<=2; dif++){

      //loop over "swap refs". iswp=iHref means don't swap
      for(int iswp=iswp_min; iswp<=iHref; iswp++){

        JplGpsData data;
        data.makeACopyOf(master_data);
        data.differenceData(dif); //either 0,1,2
        if(dif>0)data.calculateStdDev();

        //Swap reference clock (or not)
        if(iswp<iHref){
          //bool use_svn=false;
          //if(dif==0)
          bool use_svn=true; //need to use same svn for (0,1,2)!
          std::string s_clk=refname[iswp].substr(0,2);
          std::string s_blk=refname[iswp].substr(2,refname[iswp].length());
          bool swpok=data.swapReference(s_clk,s_blk,use_svn);
          if(!swpok)continue; //didn't find reference, try next
        }

        //identify the reference clock
        int iref=NUMREFS-1; //reference index. default "other"
        {
          std::string ref="X";
          if(data.refblk=="AR")ref=data.refprn;   //H-maser is reference
          else ref=data.refclk+data.refblk;       //sat clock is reference
          for(int ir=0; ir<NUMREFS-1; ir++){
            if(ref==refname[ir]){
              iref=ir;
              break;
            }
          }
        }

        //loop through all satellite clocks:
        for(int i=data.num_receivers; i<data.num_clocks; i++){

          //identify clock type: //0=Rb, 1=Cs?
          int iC=0;
          if(data.clk[i]=="Rb")iC=0;
          else if(data.clk[i]=="Cs")iC=1;
          else continue;

          //If we're swapping the ref clock, don't mix clock types!
          //And, skip the svn of the ref clock!
          if(iswp<iHref){
            if(data.clk[i]!=data.refclk)continue;
            if(data.svn[i]==data.refsvn)continue;
          }

          //identify SVN.
          //NB: index offset by min_svn (allows less wasted array)
          int isvn=0;
          isvn=std::stoi(data.svn[i]);
          if(isvn==0)continue;
          else isvn-=min_svn;



          //Exlucde 'bad' clocks:
          //nb: which clocks were are excluded is worked out above
          if(bexcl[i]){
            if(dif==0) excldays[iC][iref][isvn]++; //only count once
            continue;
          }

          //count the clocks (only dif=0)
          if(dif==0) clockdays[iC][iref][isvn]++;

          //Transfer s.d.
          if(dif>0) SD[dif-1][iC][iref][isvn]+=data.sdev[i];
          if(sd_only)continue;

          //transfer data to array (efficiency)
          double d[NUMEPOCHS];
          for(int j=0; j<NUMEPOCHS; j++) d[j]=data.bias[i][j];

          // *** Calculate the noise profile functions/properties ***

          //Power Sectrum:
          if(do_psd){
            double psd[PSDPOINTS];

            if(dif==0)
              powerSpectrumFFT(d,NUMEPOCHS,psd,PSDPOINTS,dif,work_0,real_wt_0);
            else if(dif==1)
              powerSpectrumFFT(d,NUMEPOCHS,psd,PSDPOINTS,dif,work_1,real_wt_1);
            else
              powerSpectrumFFT(d,NUMEPOCHS,psd,PSDPOINTS,dif,work_2,real_wt_2);

            //transfer calc'd functions to 'global' arrays
            for(int k=0; k<PSDPOINTS; k++)
              PSD[dif][iC][iref][isvn][k]+=psd[k];
          }

          //Autocorrelation:
          if(do_acf){
            double acf[ACFPOINTS]; //nb: can have smaller number!
            autoCorrelation(d,NUMEPOCHS,acf,ACFPOINTS,dif,true);
            //transfer calc'd functions to 'global' arrays
            for(int t=0; t<ACFPOINTS; t++)
              ACF[dif][iC][iref][isvn][t]+=acf[t];
          }

          //Allan Variance:
          if(do_avr&&(dif==0)){
            //AVAR only for un-differenced data!
            double av[AVRPOINTS]; //nb: can have smaller number!
            allanVariance(d,NUMEPOCHS,av,AVRPOINTS,false);
            for(int t=0; t<AVRPOINTS; t++)
              AVR[iC][iref][isvn][t]+=av[t];
          }

          //Historgram:
          if(do_hist&&dif>0){
            double hist[NUMBINS]={0}; //must be zerod!
            histogram(d,NUMEPOCHS,hist,NUMBINS,bin_width,dif);
            //transfer calc'd functions to 'global' arrays
            for(int b=0; b<NUMBINS; b++)
              HIST[dif-1][iC][iref][isvn][b]+=hist[b];
          }

          //Cross-Clock Correlation (at 0 lag only)
          if(do_clkcor&&dif>0){
            int t_num_clk=0;
            double t_clk_corr_tot=0;
            double s1=data.sdev[i];
            for(int i2=data.num_receivers; i2<data.num_clocks; i2++){
              //Only include same clock/block
              if(data.clk[i2]!=data.clk[i]) continue;
              if(data.blk[i2]!=data.blk[i]) continue;
              if(i2==i) continue; //ignore self correlation
              t_num_clk++;
              double s2=data.sdev[i2];
              double t_clk_corr=0;
              int t_num_eps=0;
              for(int j=0; j<data.num_epochs; j++){
                if(data.ferr[i][j]*data.ferr[i2][j]==0)continue;
                double d1=data.bias[i][j];
                double d2=data.bias[i2][j];
                t_num_eps++;
                t_clk_corr += d1*d2;
              }
              if(t_num_eps*s2!=0)
                t_clk_corr_tot+=t_clk_corr/(t_num_eps*s2);
            }
            if(t_num_clk*s1!=0)
              CLKCOR[dif-1][iC][iref][isvn]+=t_clk_corr_tot/(t_num_clk*s1);
          }

        }//END loop over (sat) clocks

      }//END swap ref
    }//END loop over dif

    //Clear the memory ascociated with GSL FFT workspace
    gsl_fft_real_workspace_free (work_0);
    gsl_fft_real_workspace_free (work_1);
    gsl_fft_real_workspace_free (work_2);

  }//END day
  MSC_progressBar("Processing Noise:",40);
  //*************************************

  //Clear the memory ascociated with GSL FFT trig tables
  gsl_fft_real_wavetable_free (real_wt_0);
  gsl_fft_real_wavetable_free (real_wt_1);
  gsl_fft_real_wavetable_free (real_wt_2);

  //Output the results!
  // -> average over SVN's, reference clocks (H only)
  // work out block! (by SVN)

  const int NUMCLKBLK=7;   //XXX may need updateing in future!
    // RbII, IIA, R, F
    // CsII, IIA, IIF
  if(NUMCLKBLK != iHref)std::cout<<"\nFAIL 368!! XXX ??\n\n";

  //Arrays to hold the results that are averaged over SVNs:
  //NOTE: now ClockBlock, not just Clock!
  float PSDsvn[NUMDIFF][NUMCLKBLK][NUMREFS][PSDPOINTS]={0};
  float ACFsvn[NUMDIFF][NUMCLKBLK][NUMREFS][ACFPOINTS]={0};
  float AVRsvn[NUMCLKBLK][NUMREFS][AVRPOINTS]={0};
  float SDsvn[2][NUMCLKBLK][NUMREFS]={0};
  static double HISTsvn[2][NUMCLKBLK][NUMREFS][NUMBINS]={0};
  int clockdayssvn[NUMCLKBLK][NUMREFS]={0};
  int excldayssvn[NUMCLKBLK][NUMREFS]={0};

  float CLKCORsvn[2][NUMCLKBLK][NUMREFS]={0};

  //Arrays to hold the results that are averaged over SVNs AND refs (H-only):
  float PSDsvnH[NUMDIFF][NUMCLKBLK][PSDPOINTS]={0};
  float ACFsvnH[NUMDIFF][NUMCLKBLK][ACFPOINTS]={0};
  float AVRsvnH[NUMCLKBLK][AVRPOINTS]={0};
  float SDsvnH[2][NUMCLKBLK]={0};
  static double HISTsvnH[2][NUMCLKBLK][NUMBINS]={0};
  int clockdayssvnH[NUMCLKBLK]={0};
  int excldayssvnH[NUMCLKBLK]={0};

  float CLKCORsvnH[2][NUMCLKBLK]={0};

  //NOTE: remember that svn is off-set!

  //Sum up all the functions, and average over days/SVNs/references:
  std::cout<<"Done reading the JPL files, now averaging the functions.."
            <<std::flush;
  for(int dif=0; dif<3; dif++){
    for(int nr=0; nr<NUMREFS; nr++){
      for(int iC=0; iC<2; iC++){
        for(int svn=0; svn<NUMSVNS; svn++){
          int clock_days=clockdays[iC][nr][svn];
          int excl_days =excldays[iC][nr][svn];
          if(clock_days==0 && excl_days==0) continue;

          //work out iCB (block+clock)
          int iCB=-1;
          if(iC==0){//Rb
            int in_svn=svn+min_svn; //remember the SVN off-set
            std::string blk=svnToBlock(in_svn);
            if(blk=="II")iCB=0;
            else if(blk=="IIA")iCB=1;
            else if(blk=="IIR")iCB=2;
            else if(blk=="IIF")iCB=3;
            else continue;
            //else std::cout<<"FAIL! 258-1!? "<<in_svn<<"\n";
          }else{//Cs
            int in_svn=svn+min_svn; //remember the SVN off-set
            std::string blk=svnToBlock(in_svn);
            if(blk=="II")iCB=4;
            else if(blk=="IIA")iCB=5;
            else if(blk=="IIF")iCB=6;
            else if(blk=="IIR")continue; //no CsIIR clocks
            else continue;
            //else std::cout<<"FAIL! 263-1!? "<<in_svn<<"\n";
          }
          if(iCB==-1)return 1;

          if(dif==0){
            //sum up the number of clocks (only once! [dif])
            if(nr>=iHref){
              clockdayssvnH[iCB]+=clock_days;
              excldayssvnH[iCB]+=excl_days;
            }
            clockdayssvn[iCB][nr]+=clock_days;
            excldayssvn[iCB][nr]+=excl_days;
          }
          if(clock_days==0) continue;
          //Average the 'main' functions (devide by # clocks) over days
          // (each SVN, each reference)
          // and sum up (but not devide yet) the avd'd over svn (etc.) fns
          //ACF:
          for(int j=0; j<ACFPOINTS; j++){
            if(nr>=iHref)ACFsvnH[dif][iCB][j]+=ACF[dif][iC][nr][svn][j];
            ACFsvn[dif][iCB][nr][j]+=ACF[dif][iC][nr][svn][j];
            ACF[dif][iC][nr][svn][j]/=clock_days;
          }
          //PSD:
          for(int j=0; j<PSDPOINTS; j++){
            if(nr>=iHref)PSDsvnH[dif][iCB][j]+=PSD[dif][iC][nr][svn][j];
            PSDsvn[dif][iCB][nr][j]+=PSD[dif][iC][nr][svn][j];
            PSD[dif][iC][nr][svn][j]/=clock_days;
          }
          if(dif==0){
            //AVAR:
            for(int j=0; j<AVRPOINTS; j++){
              if(nr>=iHref)AVRsvnH[iCB][j]+=AVR[iC][nr][svn][j];
              AVRsvn[iCB][nr][j]+=AVR[iC][nr][svn][j];
              AVR[iC][nr][svn][j]/=clock_days;
            }
          }else{ //dif=1,2
            //Cross-Clock Correlation:
            if(nr>=iHref)CLKCORsvnH[dif-1][iCB]+=CLKCOR[dif-1][iC][nr][svn];
            CLKCORsvn[dif-1][iCB][nr]+=CLKCOR[dif-1][iC][nr][svn];
            CLKCOR[dif-1][iC][nr][svn]/=clock_days;
            //HIST:
            for(int j=0; j<NUMBINS; j++){
              if(nr>=iHref)HISTsvnH[dif-1][iCB][j]+=HIST[dif-1][iC][nr][svn][j];
              HISTsvn[dif-1][iCB][nr][j]+=HIST[dif-1][iC][nr][svn][j];
              HIST[dif-1][iC][nr][svn][j]/=clock_days;
            }
            //SD:
            if(nr>=iHref)SDsvnH[dif-1][iCB]+=SD[dif-1][iC][nr][svn];
            SDsvn[dif-1][iCB][nr]+=SD[dif-1][iC][nr][svn];
            SD[dif-1][iC][nr][svn]/=clock_days;
          }
        }//svn
      }//Cs/Rb
      //Average over SVNs (for each reference):
      for(int iCB=0; iCB<NUMCLKBLK; iCB++){
        int clk_days_svn=clockdayssvn[iCB][nr];
        if(clk_days_svn==0)continue;
        //ACF, PSD:
        for(int j=0; j<ACFPOINTS; j++) ACFsvn[dif][iCB][nr][j]/=clk_days_svn;
        for(int j=0; j<PSDPOINTS; j++) PSDsvn[dif][iCB][nr][j]/=clk_days_svn;
        if(dif==0){
          //AVAR:
          for(int j=0; j<AVRPOINTS; j++) AVRsvn[iCB][nr][j]/=clk_days_svn;
        }else{ //dif=1,2
          //Cross-clock correlation:
          CLKCORsvn[dif-1][iCB][nr]/=clk_days_svn;
          //HIST, SD:
          for(int j=0; j<NUMBINS; j++)HISTsvn[dif-1][iCB][nr][j]/=clk_days_svn;
          SDsvn[dif-1][iCB][nr]/=clk_days_svn;
        }
      }//clk-blk
    }//ref
    //Average over the H-maser reference clocks:
    for(int iCB=0; iCB<NUMCLKBLK; iCB++){
      int clk_days_H=clockdayssvnH[iCB];
      if(clk_days_H==0)continue;
      //ACF, PSD:
      for(int j=0; j<ACFPOINTS; j++) ACFsvnH[dif][iCB][j]/=clk_days_H;
      for(int j=0; j<PSDPOINTS; j++) PSDsvnH[dif][iCB][j]/=clk_days_H;
      if(dif==0){
        //AVAR:
        for(int j=0; j<AVRPOINTS; j++) AVRsvnH[iCB][j]/=clk_days_H;
      }else{ //dif=1,2
        //Cross-clock correlation:
        CLKCORsvnH[dif-1][iCB]/=clk_days_H;
        //HIST, SD:
        for(int j=0; j<NUMBINS; j++) HISTsvnH[dif-1][iCB][j]/=clk_days_H;
        SDsvnH[dif-1][iCB]/=clk_days_H;
      }
    }//clk-blk
  }//dif
  std::cout<<"done.\n";


  //print breif summary to screen
  int total_clock_days=0;
  int total_excl_days=0;
  printf("Total of:\n");
  for(int icb=0; icb<7; icb++){
    total_clock_days+=clockdayssvnH[icb];
    total_excl_days+=excldayssvnH[icb];
    if(clockdayssvnH[icb]==0&&excldayssvnH[icb]==0)continue;
    printf(" %5i %5s (sd=%.4f ns).  Excluded %i\n",
          clockdayssvnH[icb],refname[icb].c_str(),SDsvnH[0][icb],
          excldayssvnH[icb]);
          //nb: use ref-name, because it matches clk-blk, not actually ref!
  }
  printf(" %i overall clock-days included, %i excluded\n",
         total_clock_days,total_excl_days);


////////////////////////////////////////////////////////////////////////////////

  // ***** Write out the functions to file *****

  std::string str_weeks="# weeks: "+std::to_string(minweek)+"-"
                       +std::to_string(maxweek)+"\n";

  // ** For each SVN and reference ------------------------------
  for(int iC=0; iC<2; iC++){
    std::string s_clk="Rb";
    if(iC==1)s_clk="Cs";
    for(int nr=0; nr<NUMREFS; nr++){
      for(int svn=0; svn<NUMSVNS; svn++){
        if(clockdays[iC][nr][svn]==0)continue;
        //clock name, and svn strings ready:
        //remember the SVN off-set
        std::string clk_blk=s_clk+svnToBlock(svn+min_svn);
        std::string str_svn=std::to_string(svn+min_svn);
        std::ofstream ofile;
        std::string ofname_base=odir+clk_blk+"-"+str_svn+"-"+refname[nr]
                                +out_label+".";
        std::string str_days="# "+std::to_string(clockdays[iC][nr][svn])
                            +" clock-days, "
                            +std::to_string(excldays[iC][nr][svn])
                            +" excluded\n";

        //Writes out the ACF (for each SVN)
        if(do_acf){
          std::string header1 ="# Std deviation (dif=1,2), "
                               "Autocorrelation (0,1,2)\n";
          std::string ofname=ofname_base+"acf";
          ofile.open(ofname.c_str());
          ofile<<"# "<<ofname<<"\n";
          ofile<<header1;
          ofile<<str_weeks;
          ofile<<str_days;
          //output the s.d.:
          ofile<<"# sd:\n";
          ofile<<SD[0][iC][nr][svn]<<" "<<SD[1][iC][nr][svn]<<"\n";
          //output the ACF
          ofile<<"# ACF: "<<ACFPOINTS;
          for(int j=0; j<ACFPOINTS; j++){
            ofile<<"\n"<<ACF[0][iC][nr][svn][j]<<" "
                 <<ACF[1][iC][nr][svn][j]<<" "<<ACF[2][iC][nr][svn][j];
          }
          ofile.close();
        }

        //Writes out the PSD (for each SVN)
        if(do_psd){
          std::string header1 ="# Power spectrum* (dif=0,1,2) "
                               "(really just |FT|^2)\n"
                               "# NOTE: only the first half! (symmetric)\n";
          std::string ofname=ofname_base+"psd";
          ofile.open(ofname.c_str());
          ofile<<"# "<<ofname<<"\n";
          ofile<<header1;
          ofile<<str_weeks;
          ofile<<str_days;
          //output the PSD
          ofile<<"# PSD: "<<PSDPOINTS;
          for(int j=0; j<PSDPOINTS; j++){
            ofile<<"\n"<<PSD[0][iC][nr][svn][j]<<" "
                 <<PSD[1][iC][nr][svn][j]<<" "<<PSD[2][iC][nr][svn][j];
          }
          ofile.close();
        }

        //Writes out the AVAR (for each SVN)
        if(do_avr){
          std::string header1 ="# Allan Variance (dif=0)\n"
                               "# Units: [ns^2 / epochs^2].\n"
                               "# --> to convert, multiply by 1.111x10^-21\n";
          std::string ofname=ofname_base+"avr";
          ofile.open(ofname.c_str());
          ofile<<"# "<<ofname<<"\n";
          ofile<<header1;
          ofile<<str_weeks;
          ofile<<str_days;
          //output the AVAR
          ofile<<"# AVAR: "<<AVRPOINTS-1; //'zero' point doesn't count!
          for(int j=1; j<AVRPOINTS; j++){
            ofile<<"\n"<<AVR[iC][nr][svn][j];
          }
          ofile.close();
        }

        //Writes out the HISTogram (for each SVN)
        if(do_hist){
          std::string header1 =
            "# Histogram (dif=1,2)\n# Units: per bin per epoch.\n";
          std::string ofname=ofname_base+"hst";
          ofile.open(ofname.c_str());
          ofile<<"# "<<ofname<<"\n";
          ofile<<header1;
          ofile<<"# bin width: "<<bin_width<<"\n";
          ofile<<str_weeks;
          ofile<<str_days;
          //output the HIST
          ofile<<"# HIST: "<<NUMBINS;
          for(int j=1; j<NUMBINS; j++){
            ofile<<"\n"<<HIST[0][iC][nr][svn][j]<<" "<<HIST[1][iC][nr][svn][j];
          }
          ofile.close();
        }

      }//SVNs
    }//References
  }//End Cs/Rb


  // ** For each reference, averaged over SVNs -----------------------
  for(int iCB=0; iCB<7; iCB++){
    for(int nr=0; nr<NUMREFS; nr++){
      if(clockdayssvn[iCB][nr]==0)continue;
      //clock name, and svn strings ready:
      //remember the SVN off-set
      std::string clk_blk=refname[iCB]; //not reference, just same array!
      std::ofstream ofile;
      std::string ofname_base=odir+clk_blk+"-av-"+refname[nr]
                              +out_label+".";
      std::string str_days="# "+std::to_string(clockdayssvn[iCB][nr])
                          +" clock-days, "
                          +std::to_string(excldayssvn[iCB][nr])
                          +" excluded\n";

      //Writes out the ACF (for each reference, avgd over svns)
      if(do_acf){
        std::string header1 ="# Std deviation (dif=1,2), "
                             "Autocorrelation (0,1,2)\n";
        std::string ofname=ofname_base+"acf";
        ofile.open(ofname.c_str());
        ofile<<"# "<<ofname<<"\n";
        ofile<<header1;
        ofile<<str_weeks;
        ofile<<str_days;
        //output the s.d.:
        ofile<<"# sd:\n";
        ofile<<SDsvn[0][iCB][nr]<<" "<<SDsvn[1][iCB][nr]<<"\n";
        //output the ACF
        ofile<<"# ACF: "<<ACFPOINTS;
        for(int j=0; j<ACFPOINTS; j++){
          ofile<<"\n"<<ACFsvn[0][iCB][nr][j]<<" "
               <<ACFsvn[1][iCB][nr][j]<<" "<<ACFsvn[2][iCB][nr][j];
        }
        ofile.close();
      }

      //Writes out the PSD (for each reference, avgd over svns)
      if(do_psd){
        std::string header1 ="# Power spectrum* (dif=0,1,2) "
                             "(really just |FT|^2)\n"
                             "# NOTE: only the first half! (symmetric)\n";
        std::string ofname=ofname_base+"psd";
        ofile.open(ofname.c_str());
        ofile<<"# "<<ofname<<"\n";
        ofile<<header1;
        ofile<<str_weeks;
        ofile<<str_days;
        //output the PSD
        ofile<<"# PSD: "<<PSDPOINTS;
        for(int j=0; j<PSDPOINTS; j++){
          ofile<<"\n"<<PSDsvn[0][iCB][nr][j]<<" "
               <<PSDsvn[1][iCB][nr][j]<<" "<<PSDsvn[2][iCB][nr][j];
        }
        ofile.close();
      }

      //Writes out the AVAR (for each reference, avgd over svns)
      if(do_avr){
        std::string header1 ="# Allan Variance (dif=0)\n"
                             "# Units: [ns^2 / epochs^2].\n"
                             "# --> to convert, multiply by 1.111x10^-21\n";
        std::string ofname=ofname_base+"avr";
        ofile.open(ofname.c_str());
        ofile<<"# "<<ofname<<"\n";
        ofile<<header1;
        ofile<<str_weeks;
        ofile<<str_days;
        //output the AVAR
        ofile<<"# AVAR: "<<AVRPOINTS-1; //'zero' point doesn't count!
        for(int j=1; j<AVRPOINTS; j++){
          ofile<<"\n"<<AVRsvn[iCB][nr][j];
        }
        ofile.close();
      }

      //Writes out the HISTogram (for each reference, avgd over svns)
      if(do_hist){
        std::string header1 =
          "# Histogram (dif=1,2)\n# Units: per bin per epoch.\n";
        std::string ofname=ofname_base+"hst";
        ofile.open(ofname.c_str());
        ofile<<"# "<<ofname<<"\n";
        ofile<<header1;
        ofile<<"# bin width: "<<bin_width<<"\n";
        ofile<<str_weeks;
        ofile<<str_days;
        //output the HIST
        ofile<<"# HIST: "<<NUMBINS;
        for(int j=1; j<NUMBINS; j++){
          ofile<<"\n"<<HISTsvn[0][iCB][nr][j]<<" "<<HISTsvn[1][iCB][nr][j];
        }
        ofile.close();
      }

    }//References
  }//End Clk-Blk


  // ** Averaged over SVNs and (H-maser) references ----------------
  for(int iCB=0; iCB<7; iCB++){
    if(clockdayssvnH[iCB]==0)continue;
    //clock name, and svn strings ready:
    //remember the SVN off-set
    std::string clk_blk=refname[iCB]; //not reference, just same array!
    std::ofstream ofile;
    std::string ofname_base=odir+clk_blk+"-av-Hmas"
                            +out_label+".";
    std::string str_days="# "+std::to_string(clockdayssvnH[iCB])
                        +" clock-days, "
                        +std::to_string(excldayssvnH[iCB])
                        +" excluded\n";

    //Writes out the ACF (avgd over refs+svns)
    if(do_acf){
      std::string header1 ="# Std deviation (dif=1,2), "
                           "Autocorrelation (0,1,2)\n";
      std::string ofname=ofname_base+"acf";
      ofile.open(ofname.c_str());
      ofile<<"# "<<ofname<<"\n";
      ofile<<header1;
      ofile<<str_weeks;
      ofile<<str_days;
      //output the s.d.:
      ofile<<"# sd:\n";
      ofile<<SDsvnH[0][iCB]<<" "<<SDsvnH[1][iCB]<<"\n";
      //output the ACF
      ofile<<"# ACF: "<<ACFPOINTS;
      for(int j=0; j<ACFPOINTS; j++){
        ofile<<"\n"<<ACFsvnH[0][iCB][j]<<" "
             <<ACFsvnH[1][iCB][j]<<" "<<ACFsvnH[2][iCB][j];
      }
      ofile.close();
    }

    //Writes out the PSD (avgd over refs+svns)
    if(do_psd){
      std::string header1 ="# Power spectrum* (dif=0,1,2) "
                           "(really just |FT|^2)\n"
                           "# NOTE: only the first half! (symmetric)\n";
      std::string ofname=ofname_base+"psd";
      ofile.open(ofname.c_str());
      ofile<<"# "<<ofname<<"\n";
      ofile<<header1;
      ofile<<str_weeks;
      ofile<<str_days;
      //output the PSD
      ofile<<"# PSD: "<<PSDPOINTS;
      for(int j=0; j<PSDPOINTS; j++){
        ofile<<"\n"<<PSDsvnH[0][iCB][j]<<" "
             <<PSDsvnH[1][iCB][j]<<" "<<PSDsvnH[2][iCB][j];
      }
      ofile.close();
    }

    //Writes out the AVAR (avgd over refs+svns)
    if(do_avr){
      std::string header1 ="# Allan Variance (dif=0)\n"
                           "# Units: [ns^2 / epochs^2].\n"
                           "# --> to convert, multiply by 1.111x10^-21\n";
      std::string ofname=ofname_base+"avr";
      ofile.open(ofname.c_str());
      ofile<<"# "<<ofname<<"\n";
      ofile<<header1;
      ofile<<str_weeks;
      ofile<<str_days;
      //output the AVAR
      ofile<<"# AVAR: "<<AVRPOINTS-1; //'zero' point doesn't count!
      for(int j=1; j<AVRPOINTS; j++){
        ofile<<"\n"<<AVRsvnH[iCB][j];
      }
      ofile.close();
    }

    //Writes out the HISTogram (avgd over refs+svns)
    if(do_hist){
      std::string header1 =
        "# Histogram (dif=1,2)\n# Units: per bin per epoch.\n";
      std::string ofname=ofname_base+"hst";
      ofile.open(ofname.c_str());
      ofile<<"# "<<ofname<<"\n";
      ofile<<header1;
      ofile<<"# bin width: "<<bin_width<<"\n";
      ofile<<str_weeks;
      ofile<<str_days;
      //output the HIST
      ofile<<"# HIST: "<<NUMBINS;
      for(int j=1; j<NUMBINS; j++){
        ofile<<"\n"<<HISTsvnH[0][iCB][j]<<" "<<HISTsvnH[1][iCB][j];
      }
      ofile.close();
    }

  }//End clk-blk


  ///////////////////////////////////////////////////////////////////

  static float SDunc[2][NUMCLOCKS][NUMREFS][NUMSVNS]={0};
  //Arrays to hold the results that are averaged over SVNs:
  //NOTE: now ClockBlock, not just Clock!
  float SDuncsvn[2][NUMCLKBLK][NUMREFS]={0};
  //Arrays to hold the results that are averaged over SVNs AND refs (H-only):
  float SDuncsvnH[2][NUMCLKBLK]={0};


  // NOTE: this needs to be indentical to above loop!
  //*************************************
  //Loop through each day of data. Read it in, calculate PSD,ACF etc.
  if(do_sdUnc){

    #pragma omp parallel for if(day_para)
    for(int days=init_day; days<max_day; days++){
      MSC_progressBar("Calculating uncertainty in sd",30,days-init_day,
                      max_day-init_day,day_para);

      //convert to week+day
      int iw=(int)floor(days/7);
      int id=days - iw*7;
      if(single_day&&(id!=weekday))continue;

      //create "master" data object, that will be duplicated
      //(duplication means we can difference 0,1,2 etc)
      JplGpsData master_data;
      int jplok=master_data.readJplBiasData(path,iw,id);
      if(jplok!=0)continue;

      //remove 2-order polynomial (weigtd)
      if(remove_poly) master_data.polynomialDetrend(ipoly,iweight);

      //Exclude clocks (or not) based on missed data, SD, or
      std::vector<bool> bexcl(master_data.num_clocks);
      {
        JplGpsData test_data;
        test_data.makeACopyOf(master_data);
        test_data.differenceData(1);
        test_data.calculateStdDev();
        for(int i=master_data.num_receivers; i<master_data.num_clocks; i++){
          bexcl[i]=false;
          if(exclsd && test_data.sdev[i]>dexclsd){
            bexcl[i]=true;
            continue; //go to next clock, no point checking further!
          }
          int num_missed=0;
          for(int j=1; j<test_data.num_epochs; j++){
            if(test_data.ferr[i][j]==0)num_missed++;
            if(num_missed>=iexclfer && exclfer) bexcl[i]=true;
            if(test_data.bias[i][j]>dexclO && exclO) bexcl[i]=true;
            if(bexcl[i]) break;
          }
        }
      }

      //loop over differecing orders:
      for(int dif=1; dif<=2; dif++){ //only 1 and 2!

        //loop over "swap refs". iswp=iHref means don't swap
        for(int iswp=iswp_min; iswp<=iHref; iswp++){

          JplGpsData data;
          data.makeACopyOf(master_data);
          data.differenceData(dif); //either 0,1,2
          if(dif>0)data.calculateStdDev();

          //Swap reference clock (or not)
          if(iswp<iHref){
            bool use_svn=true; //need to use same svn for (0,1,2)!
            std::string s_clk=refname[iswp].substr(0,2);
            std::string s_blk=refname[iswp].substr(2,refname[iswp].length());
            bool swpok=data.swapReference(s_clk,s_blk,use_svn);
            if(!swpok)continue; //didn't find reference, try next
          }

          //identify the reference clock
          int iref=NUMREFS-1; //reference index. default "other"
          {
            std::string ref="X";
            if(data.refblk=="AR")ref=data.refprn;   //H-maser is reference
            else ref=data.refclk+data.refblk;       //sat clock is reference
            for(int ir=0; ir<NUMREFS-1; ir++){
              if(ref==refname[ir]){
                iref=ir;
                break;
              }
            }
          }

          //loop through all satellite clocks:
          for(int i=data.num_receivers; i<data.num_clocks; i++){

            //identify clock type: //0=Rb, 1=Cs?
            int iC=0;
            if(data.clk[i]=="Rb")iC=0;
            else if(data.clk[i]=="Cs")iC=1;
            else continue;

            //If we're swapping the ref clock, don't mix clock types!
            //And, skip the svn of the ref clock!
            if(iswp<iHref){
              if(data.clk[i]!=data.refclk)continue;
              if(data.svn[i]==data.refsvn)continue;
            }

            //identify SVN.
            //NB: index offset by min_svn (allows less wasted array)
            int isvn=0;
            isvn=std::stoi(data.svn[i]);
            if(isvn==0)continue;
            else isvn-=min_svn;

            //Exlucde 'bad' clocks:
            //nb: which clocks were are excluded is worked out above
            if(bexcl[i]){
              continue;
            }

            double avsd=SD[dif-1][iC][iref][isvn];
            SDunc[dif-1][iC][iref][isvn]+=pow(data.sdev[i]-avsd,2);
            //Note: SDunc still needs Sqrt[SDunc/N_clocks]! (afterwards)

          }//END loop over (sat) clocks
        }//END loop over 'swap references'
      }//END loop over dif
    }//END day
    MSC_progressBar("Calculating uncertainty in sd",30);

    //*************************************

    //Devide  by num_clocks, and variance -> deviation

    //Sum up all the functions, and average over days/SVNs/references:
    std::cout<<"Done calculating sd Uncertainty, now averaging over SVNs etc.."
             <<std::flush;
    for(int dif=1; dif<=2; dif++){
      for(int nr=0; nr<NUMREFS; nr++){
        for(int iC=0; iC<2; iC++){
          for(int svn=0; svn<NUMSVNS; svn++){
            int clock_days=clockdays[iC][nr][svn];
            int excl_days =excldays[iC][nr][svn];
            if(clock_days==0 && excl_days==0) continue;

            //work out iCB (block+clock)
            int iCB=-1;
            if(iC==0){//Rb
              int in_svn=svn+min_svn; //remember the SVN off-set
              std::string blk=svnToBlock(in_svn);
              if(blk=="II")iCB=0;
              else if(blk=="IIA")iCB=1;
              else if(blk=="IIR")iCB=2;
              else if(blk=="IIF")iCB=3;
              else continue;
              //else std::cout<<"FAIL! 258-2!? "<<in_svn<<"\n";
            }else{//Cs
              int in_svn=svn+min_svn; //remember the SVN off-set
              std::string blk=svnToBlock(in_svn);
              if(blk=="II")iCB=4;
              else if(blk=="IIA")iCB=5;
              else if(blk=="IIF")iCB=6;
              else if(blk=="IIR")continue; //no CsIIR clocks
              else continue;
              //else std::cout<<"FAIL! 263-2!? "<<in_svn<<"\n";
            }
            if(iCB==-1)return 1;

            if(clock_days==0) continue;
            //Average the 'main' functions (devide by # clocks) over days
            // (each SVN, each reference)
            // and sum up (but not devide yet) the avd'd over svn (etc.) fns
              //SD:
            SDuncsvn[dif-1][iCB][nr] += (
                SDunc[dif-1][iC][nr][svn] //SDunc has not been sqrt or /N yet!
              + clock_days*pow(SD[dif-1][iC][nr][svn]
                               -SDsvn[dif-1][iCB][nr],2)
            );
            if(nr>=iHref) SDuncsvnH[dif-1][iCB] += (
                SDunc[dif-1][iC][nr][svn] //SDunc has not been sqrt or /N yet!
              + clock_days*pow(SD[dif-1][iC][nr][svn]
                               -SDsvnH[dif-1][iCB],2)
            );
            //Now, take sqrt and devide the "each svn" SDunc
            SDunc[dif-1][iC][nr][svn]=
              sqrt(SDunc[dif-1][iC][nr][svn]/clock_days);
          }//svn
        }//Cs/Rb
        //Average over SVNs (for each reference):
        for(int iCB=0; iCB<NUMCLKBLK; iCB++){
          int clk_days_svn=clockdayssvn[iCB][nr];
          if(clk_days_svn==0)continue;
          SDuncsvn[dif-1][iCB][nr]=sqrt(SDuncsvn[dif-1][iCB][nr]/clk_days_svn);
        }//clk-blk
      }//ref
      //Average over the H-maser reference clocks:
      for(int iCB=0; iCB<NUMCLKBLK; iCB++){
        int clk_days_H=clockdayssvnH[iCB];
        if(clk_days_H==0)continue;
        SDuncsvnH[dif-1][iCB]=sqrt(SDuncsvnH[dif-1][iCB]/clk_days_H);
      }//clk-blk
    }//dif
    std::cout<<"done.\n";

  }// end if do uncertainty


  //Write the Summary to file for standard deviations with uncertianty
  std::ofstream s2File;
  std::string sum2Name;
  if(do_sdUnc)sum2Name=odir+"sdWunc"+out_label+".txt";
  else        sum2Name=odir+"sd"+out_label+".txt";
  s2File.open(sum2Name.c_str());
  if(do_sdUnc)s2File<<"# Summary (with uncertainties): "<<out_label<<"\n";
  else        s2File<<"# Summary (no uncertainties): "<<out_label<<"\n";
  s2File<<"# weeks "<<minweek<<" - "<<maxweek<<" \n";
  s2File<<"# clk blk svn ref days excluded s1 +/- s2 +/-\n";
  s2File.precision(5);
  for(int iC=0; iC<2; iC++){
    for(int isvn=0; isvn<NUMSVNS; isvn++){
      for(int nr=0; nr<NUMREFS; nr++){
        if((clockdays[iC][nr][isvn]==0)&&(excldays[iC][nr][isvn]==0))continue;
        int svn=isvn+min_svn;
        std::string blk=svnToBlock(svn);
        std::string clk="Rb"; if(iC==1)clk="Cs";
        s2File<<clk<<" "<<blk<<" "<<svn<<" "<<refname[nr]<<" "
              <<clockdays[iC][nr][isvn]<<" "<<excldays[iC][nr][isvn]<<"  "
              <<SD[0][iC][nr][isvn]<<" "<<SDunc[0][iC][nr][isvn]<<"  "
              <<SD[1][iC][nr][isvn]<<" "<<SDunc[1][iC][nr][isvn]<<"\n";
      }//ref
    }//svn
  }
  s2File<<"# Averaged over all SVNs\n";
  for(int iCB=0; iCB<7; iCB++){
    if(clockdayssvnH[iCB]==0&&excldayssvnH[iCB]==0)continue;
    std::string clk=refname[iCB].substr(0,2);
    std::string blk=refname[iCB].substr(2,refname[iCB].length());
    for(int nr=0; nr<NUMREFS; nr++){
      if((clockdayssvn[iCB][nr]==0)&&(excldayssvn[iCB][nr]==0))continue;
      s2File<<clk<<" "<<blk<<" av "<<refname[nr]<<" "
            <<clockdayssvn[iCB][nr]<<" "<<excldayssvn[iCB][nr]<<"  "
            <<SDsvn[0][iCB][nr]<<" "<<SDuncsvn[0][iCB][nr]<<"  "
            <<SDsvn[1][iCB][nr]<<" "<<SDuncsvn[1][iCB][nr]<<"\n";
    }
    s2File<<clk<<" "<<blk<<" av Hmas "
          <<clockdayssvnH[iCB]<<" "<<excldayssvnH[iCB]<<"  "
          <<SDsvnH[0][iCB]<<" "<<SDuncsvnH[0][iCB]<<"  "
          <<SDsvnH[1][iCB]<<" "<<SDuncsvnH[1][iCB]<<"\n";
  }
  s2File.close();



  //Write the Summary to file for standard deviations with uncertianty
  std::ofstream cccFile;
  std::string cccName;
  cccName=odir+"crossClkCorr"+out_label+".txt";
  cccFile.open(cccName.c_str());
  cccFile<<"# Cross-Clock Correlation: "<<out_label<<"\n";
  cccFile<<"# weeks "<<minweek<<" - "<<maxweek<<" \n";
  cccFile<<"# clk blk svn ref days excluded s1 s2 \n";
  cccFile.precision(5);
  for(int iC=0; iC<2; iC++){
    for(int isvn=0; isvn<NUMSVNS; isvn++){
      for(int nr=0; nr<NUMREFS; nr++){
        if((clockdays[iC][nr][isvn]==0)&&(excldays[iC][nr][isvn]==0))continue;
        int svn=isvn+min_svn;
        std::string blk=svnToBlock(svn);
        std::string clk="Rb"; if(iC==1)clk="Cs";
        cccFile<<clk<<" "<<blk<<" "<<svn<<" "<<refname[nr]<<" "
              <<clockdays[iC][nr][isvn]<<" "<<excldays[iC][nr][isvn]<<"  "
              <<CLKCOR[0][iC][nr][isvn]<<" "
              <<CLKCOR[1][iC][nr][isvn]<<" "<<"\n";
      }//ref
    }//svn
  }
  cccFile<<"# Averaged over all SVNs\n";
  for(int iCB=0; iCB<7; iCB++){
    if(clockdayssvnH[iCB]==0&&excldayssvnH[iCB]==0)continue;
    std::string clk=refname[iCB].substr(0,2);
    std::string blk=refname[iCB].substr(2,refname[iCB].length());
    for(int nr=0; nr<NUMREFS; nr++){
      if((clockdayssvn[iCB][nr]==0)&&(excldayssvn[iCB][nr]==0))continue;
      cccFile<<clk<<" "<<blk<<" av "<<refname[nr]<<" "
              <<clockdayssvn[iCB][nr]<<" "<<excldayssvn[iCB][nr]<<"  "
              <<CLKCORsvn[0][iCB][nr]<<" "
              <<CLKCORsvn[1][iCB][nr]<<" "<<"\n";
    }
    cccFile<<clk<<" "<<blk<<" av Hmas "
              <<clockdayssvnH[iCB]<<" "
              <<excldayssvnH[iCB]<<"  "
              <<CLKCORsvnH[0][iCB]<<" "
              <<CLKCORsvnH[1][iCB]<<" "<<"\n";
  }
  cccFile.close();











  //change to H:mm:ss
  time (&end);
  double dif = difftime (end,start);
  printf ("Total elasped time is %.0f seconds.\n", dif );

  return 0;
}










//******************************************************************************
std::string svnToBlock(int in_svn){

  //minimum and maximum SVNs for each block:
  int IImin=13,IImax=21;  //17 is max in our data
  int IIAmin=22,IIAmax=40;
  int IIRmin=41,IIRmax=61;
  int IIFmin=62,IIFmax=73;//Cs only has 63.. may cause problems?

  if(in_svn<IImin){
    return "I"; //nb: this should never happen, may cause errors above.
  }else if(in_svn>=IImin && in_svn<=IImax){
    return "II";
  }else if(in_svn>=IIAmin && in_svn<=IIAmax){
    return "IIA";
  }else if(in_svn>=IIRmin && in_svn<=IIRmax){
    return "IIR"; //Note: break into IIR and IIRM ??
  }else if(in_svn>=IIFmin && in_svn<=IIFmax){
    //any 'higher' SVN will be called IIF!
    return "IIF";
  }else if(in_svn>IIFmax){
    return "III"; //this is a future-proofing "catch-all" - needs checking!
  }

  return "??";
}
