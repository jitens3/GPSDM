/*
PROGRAM: searchClass

Still under developement.

nts in Chi^2 sum) is important.
Open question...


=======Change log=====
...

****** To Do ******


*/

#include "SearchClass.h"  

SearchClass::SearchClass(std::string path_to_data , int NMC , 
        int jplweek , int jplday, 
        bool detrend, bool difference, 
        double o , double s , double h) 
        : JplGpsData(path_to_data, jplweek, 
                      jplday, detrend, 
                      difference),
        odds_cut(o),
        SNR_cut(s),
        hmin(h)
  {
      m_snr = new double [MAXEPOCHS];
      m_snr_temp = new double[NMC*MAXEPOCHS];
      tday_snr = new double[NMC*MAXEPOCHS];       // 12/17/20 tday
      m_odds = new double [MAXEPOCHS];
      m_Einv_e = new double [MAXJW*MAXJW*num_clocks*num_clocks];
      m_Einv_p = new double [MAXJW*MAXJW*num_clocks];
      //m_Einv_e = new double [MAXJW*MAXJW*29*29];
      //m_Einv_p = new double [MAXJW*MAXJW*29];
      m_chi = new double [MAXEPOCHS];
      m_chi_temp = new double [NMC*MAXEPOCHS];

      for(int i=0; i<MAXEPOCHS; i++){
        *(m_snr+i)=0;
        *(m_odds+i)=0;
      }

  //    m_investigate = 0;

      m_likelihood = new double [MAXEPOCHS*NMC*2]; //we need to print max size here every day so 
                   //we know that we can store all the necessary vals
      dHd = new double [MAXEPOCHS];
      
      if(m_snr == NULL || m_odds == NULL || m_Einv_e == NULL || m_Einv_p == NULL 
              || m_snr_temp == NULL || dHd == NULL || m_likelihood == NULL){
        std::cout << "\n\n ERROR! Class memory allocation failed \n\n";
      }
  };

  SearchClass::~SearchClass(){
    delete [] m_snr;
    delete [] m_snr_temp;
    delete [] tday_snr;         // 12/17/20 tday
    delete [] m_odds;
    delete [] m_Einv_e;
    delete [] m_Einv_p;
    delete [] m_likelihood;
    delete [] dHd;
    delete [] m_chi;
    delete [] m_chi_temp;
    //free(m_Einv);
  }

int SearchClass::MCSearchStats(
  std::string model,
  SearchClass &obj,
  int icov, int JW,
  int jbeg, int jend,
  int NMC,
  NumericCdfInverse v_prior, 
  NumericCdfInverse psi_prior,
  int exact,
  int week,
  int day,
  int ilimit_cores,
  int sanity
)
{
  std::cout << "\n\n :::::::: SEARCH ::::::::\n\n";
  int check = 0;
  const int hJW = (JW-1)/2; 
  if(jbeg<=hJW) jbeg=hJW;
  if(jend>= num_epochs -hJW) jend=num_epochs-hJW-1;
  jend = 2849; // added by tday 06/24/21

  double weight_norm=0;
  for(int i = num_receivers; i < num_clocks; i++){
    weight_norm += 1 / sdev[i];
  };
 
  double acf_cut=0.025;
  int max_cov=61;

if(check == 0){
  if(icov == 1){
    if(exact == 0){
        std::cout << "hello before genHijl" << "\n";
      genHijl(obj, m_Einv_p, JW, max_cov,acf_cut);
        //======================== test, delete later ==========================
        /**/std::ofstream Hmat("./results2/Hn61matrices.txt");
        /**/    //Hmat<<"Hessian:\n";
        /**/    for(int i=0;i<num_clocks;i++){
            /**/      //Hmat<<"Clock "<<i<<"\n";
            /**/      for(int j=0;j<JW;j++){
                /**/        for(int k=0;k<JW;k++){
                    /**/          Hmat<< *(m_Einv_p + i*JW*JW + j*JW + k) <<" ";
                    /**/        }
                /**/        Hmat<<"\n";
                /**/      }
            /**/      Hmat<<"\n";
            /**/    }
        /**/Hmat.close();
        //======================== test, delete later ==========================
        
    }else{
      calcEinv(m_Einv_e, obj, JW, icov, max_cov, acf_cut, day, week,0.005,1, sanity);
      // reference stdev set to 0.005
      // calc ccc by hand = 1 does not use inputted reference stdev
    }
  }//end icov if
//invert both the exact and perturbative 
}
    std::cout << "hello after genHijl" << "\n";


  double logNorm = likelihoodLogNorm(icov, num_clocks, JW, sdev, m_Einv_p);
  double sb0 = b0[0];

    std::cout << "num_clocks: " << num_clocks << "\n";
    std::cout << "num_satellites: " << num_satellites << "\n";
    
    std::cout << "this is start of template bank: " << "\n";
    int tevent = 1194;
/*    
    templatebanksnr(model, *this, icov, 1, 10*week+day, tevent, JW, NMC,
    hmin, v_prior, psi_prior, exact, m_Einv_e);
*/

  omp_set_num_threads(ilimit_cores);

  printf("Parallelised (OMP) using %i cores\n\n",ilimit_cores);

 // clock_t totalStart=clock();
double start;
double end;
    
    std::string chiPath = "./results3/chi2_" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream chi(chiPath);
    
    std::string chi2Path = "./results2/chihit" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream chi2(chi2Path);
    
    std::string chi3Path = "./results3/t2chids" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream test(chi3Path);
    
    std::string chi4Path = "./results3/chi1k" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream chi3(chi4Path);

    std::string chi5Path = "./results3/chiPEs" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream chi4(chi5Path);
    
    std::string sigPath = "./results2/tnewsignal_1024.out";
    std::ofstream smatt(sigPath);
    
    int increment;
    increment = JW;

start = omp_get_wtime();
//  #pragma omp parallel for shared(jend, num_clocks, obj, m_snr, m_Einv_e, num_satellites, m_Einv_p, sdev, icov, sb0, dHd, logNorm, NMC, v_prior, psi_prior, model, vdiff, keff, hmin, m_odds, check)
  #pragma omp parallel for default(shared) shared(jbeg, jend, obj, icov, sb0, logNorm, NMC, v_prior, psi_prior, model, check)
//  for(int j0 = jbeg; j0 <= jend; j0+=increment){
  for(int j0 = jbeg; j0 <= jend; j0++){
    MSC_progressBar("Epochs",30,j0,jend,1);
      //std::cout << "hello tyler" << "\n";
      //std::cout << "num clocks: " << num_clocks << "\n";

    double dataStream[num_clocks][MAXJW];
    double rsat[num_clocks][3];
    double rref[3];
    int jmin;

    prepareData(hJW,j0,jmin,dataStream,rsat,rref, obj); 

    // this makes a vector of d^T*E^-1 11/26/19 tday
    double vec[num_clocks*MAXJW];
    for(int i=0;i<num_clocks*JW;i++){
          vec[i] = 0;
      }
    calcdE(dataStream, m_Einv_e, vec, JW, num_clocks);      // this function makes the vector
    // end of making new vector on 11/26/19 tday
   
    double chi_dd = 0.;

    if(exact == 1){
        // num_satellites = num Cs clocks and not number of total clocks
      //chi_dd = calcChi_BLAS(m_Einv_e, *dataStream, *dataStream, num_satellites, JW);
      chi_dd = calcChi_BLAS(m_Einv_e, *dataStream, *dataStream, num_clocks, JW);
      //XXX why are you passing *dataStream?? isn't 'dataStram' a ptr to 1st element of array?
    }
    else{
        //std::cout << "num_satellites: " << num_satellites << "\n";
      //chi_dd = calcChi(dataStream, dataStream, m_Einv_p, sdev, JW, num_satellites, icov);
      //chi_dd += calcChiW(dataStream,dataStream,sb0,sdev,JW,num_satellites);
      chi_dd = calcChi(dataStream, dataStream, m_Einv_p, sdev, JW, num_clocks, icov);
        //std::cout << "chi-dd: " << chi_dd << "\n";
      chi_dd += calcChiW(dataStream,dataStream,sb0,sdev,JW,num_clocks);
    }
    #pragma omp critical
    {
       // test << j0 << " " << chi_dd <<std::endl;
     *(dHd+j0) = chi_dd;
    }

    double Pdm = 0;
    double logPno = (logNorm-0.5*(chi_dd));
      
      //std::cout << "hello in front of NMC loop" << "\n";

      
      
      // NMC loop
    for (int ix=0; ix< NMC; ix++){

      double t0,v,d,R,a;
      double n[3]; 

      randomParameters(v_prior,psi_prior,model,j0,t0,v,n,d,R,a);
        
        

        double signal[num_clocks][MAXJW] = {0};
      generateSignal(model,signal,rsat,rref, num_clocks, jmin, JW, vdiff, keff,
            t0, v, n, d, R, a);    //XXX forming signal later as 1's || WHAT DOES THIS MEAN
     
/*
    // injecting my own template signal, only 1 template for the 21 clocks
    
    std::string path_to_sig = "./results/injSignalTemplate13101.txt";
        
        double ss[num_clocks][MAXJW];
        
    std::ifstream inFilet;
    inFilet.open (path_to_sig.c_str());
        if(inFilet.is_open()){
            for(int i = 0; i<num_clocks; i++){
                for(int j = 0; j<JW;j++){
                    inFilet >> ss[i][j];
//                    std::cout << ss[i][j] << " ";
                    signal[i][j] = ss[i][j];
                }
//                std::cout << "\n";
            }
        }
        inFilet.close();
*/     
        
    // added on april 21 to make 1024 templates
  /*  
        if (j0 == jbeg){
            
            for(int a=0; a<num_clocks;a++){
              for(int b =0; b<JW; b++){
                  smatt << signal[a][b] << " ";
//                  std::cout << signal[a][b] << " ";
              }
                  smatt << "\n";
  //              std::cout << "\n";
            }
//            smatt.close();
//            std::cout << "NMC: " << ix << "\n";
            
        }
    */     
        // end of april 21 addition for making 1024 template bank
        

      double chi_ds=0., chi_ss=0.;

      if(exact == 1){
        // this is a new function to perform this dHs multiplication 11/26/19 tday
        //chi_ds = calcdEs(vec, signal, JW, num_clocks);
        chi_ss = calcChi_BLAS(m_Einv_e, *signal, *signal, num_clocks, JW);
        chi_ds = calcChi_BLAS(m_Einv_e, *signal, *dataStream, num_clocks, JW);
      }
      else{
        //chi_ds = calcChi(signal, dataStream, m_Einv_p, sdev, JW, num_satellites, icov);
        //chi_ss = calcChi(signal, signal, m_Einv_p, sdev, JW, num_satellites, icov);
        
        //chi_ds += calcChiW(signal, dataStream, sb0, sdev, JW, num_satellites);
        //chi_ss += calcChiW(signal, signal, sb0, sdev, JW, num_satellites);
   //       std::cout << "Hello" << "\n";
          // added 08/20/20
        chi_ds = calcChi(signal, dataStream, m_Einv_p, sdev, JW, num_clocks, icov);
        chi_ss = calcChi(signal, signal, m_Einv_p, sdev, JW, num_clocks, icov);
   //       std::cout << "chi-ss: " << chi_ss << "\n";
        
        chi_ds += calcChiW(signal, dataStream, sb0, sdev, JW, num_clocks);
        chi_ss += calcChiW(signal, signal, sb0, sdev, JW, num_clocks);
          
  //        std::cout << "chi-ss: " << chi_ss << "\n";
      }
        
        if(chi_ss<0){
          //XXX Stand-in. ALSO: "other" error??? XXX
          std::cout<<"\nERROR 345 in xx: ss=0?\n";
          continue;
        }

        double chi_val;
      #pragma omp critical
      {
          
          // need to make a new chi value  and then save to file
              
              //chi_val = chi_dd + chi_ss - 2*chi_ds;
              //chi << chi_val << "\n";
          test << j0 << " " << ix << "  " << chi_dd << "  " <<  chi_ds << "  " << chi_ss << std::endl;
              
          // end of new additions for the chi value
      double tempPdm = hAnalytic(chi_ds,chi_ss, hmin);
      Pdm += tempPdm;

      *(m_likelihood + 2*NMC*j0 + 2*ix) = chi_ds;
      *(m_likelihood + 2*NMC*j0 + 2*ix+1) = chi_ss;

	    *(m_snr_temp + j0*NMC + ix) = chi_ds / sqrt(chi_ss);
          
        double h;
          h = chi_ds/chi_ss;
          
        *(m_chi_temp + j0*NMC + ix) = chi_dd + (h*h*chi_ss) - (2*h*chi_ds);
   //       chi << *(m_chi_temp + j0*NMC + ix) << "\n";
          
          if(*(m_chi_temp + j0*NMC + ix)<=0){
            std::cout << "This is negative chi: " << *(m_chi_temp + j0*NMC + ix) << "\n";
          }
          
      }
      //std::cout << "hello end of NMC loop" << "\n";
    } // end of MCMC loop
      //std::cout << "hello out of NMC loop" << "\n";

    #pragma omp critical
    {
    Pdm/=NMC;
    *(m_odds + j0)=Pdm;
    //*(m_odds + 2*j0 + 1)= log(Pdm);
    // *(m_odds + 2*j0 + 1)= log(Pdm)+logPno;
    }

  } // end of epoch loop
    std::cout << "hello out of epoch loop" << "\n";
   //chi.close();

  
  MSC_progressBar("Epochs",30);
  end = omp_get_wtime();

for (int i = MAXJW/2+1; i < num_epochs-MAXJW/2; i++){
  for(int j = 0; j < NMC; j++){
   if(fabs(m_snr_temp[i*NMC + j]) > fabs(m_snr[i])){
          m_snr[i]=m_snr_temp[i*NMC + j];
          m_chi[i] = m_chi_temp[i*NMC + j];
  //        chi << i << " " << m_snr[i] << "  " << m_chi[i] << "\n";
          //std::cout << "\n" << i*NMC + j;
	    }
  }
    chi << i << " " << m_snr[i] << "  " << m_chi[i] << "\n";
}
   
// this is just an investigation of snr 12/17/20 tday
  //  for (int i = MAXJW/2+1; i < num_epochs-MAXJW/2; i++){
    for (int i = MAXJW/2+1; i < num_epochs-MAXJW/2; i+=increment){
      for(int j = 0; j < NMC; j++){
       //if(fabs(m_snr_temp[i*NMC + j]) > fabs(m_snr[i])){
              tday_snr[i*NMC+j]=m_snr_temp[i*NMC + j];
              chi3 << i << " " << j << " " << tday_snr[i*NMC+j] << "\n";
              //std::cout << "\n" << i*NMC + j;
       //     }
      }
    }
// end of stuff on 12/17/20 for snr investigation
    std::cout << "Hello at end of program but before setting snr to zero if chi too big" << "\n";
    
// this is for the chi-squared filter 09/16/20
    // making some bogus results just to check filter works!
 //   m_snr[100] = 8.6;
 //   m_chi[100] = 2350;
    
/* need to make a regular function for calculating the correct chi-squared threshold value
 depending on the number of clocks in the network,
 09/21/20 tday
 */
    /*
    int chi_thresh;
    if (num_clocks == 21){
        chi_thresh = 1603;
    }
    if (num_clocks == 22){
        chi_thresh = 1672;
    }
    if (num_clocks == 23){
        chi_thresh = 1740;
    }
    if (num_clocks == 24){
        chi_thresh = 1808;
    }
    if (num_clocks == 25){
        chi_thresh = 1876;
    }
    if (num_clocks == 26){
        chi_thresh = 1944;
    }
    if (num_clocks == 27){
        chi_thresh = 2012;
    }
    if (num_clocks == 28){
        chi_thresh = 2080;
    }
    if (num_clocks == 29){
        chi_thresh = 2147;
    }
    if (num_clocks == 30){
        chi_thresh = 2215;
    }
    if (num_clocks == 31){
        chi_thresh = 2282;
    }
    
    */
    // end of the chi-squared psuedo function to determine thresh value
    
    double alpha = pow(10,-10);
    double beta = 2*(1-alpha)-1;
    int chi_thresh = num_clocks*MAXJW + sqrt(2)*sqrt(2*num_clocks*MAXJW) * MFS_inverseErf(beta);
    
    int fpcount;
    fpcount = 0;
    
    int fpcount2;
    fpcount2 = 0;
    
    std::cout << "chi-thresh: " << chi_thresh << "\n";
    
 //   std::cout << "SNR: " << m_snr[100] << " " << "chi: " << m_chi[100] <<"\n";

   // adding new threshold determination from SNR-max paper (Oct 2021)
    std::string path_to_sig2 = "./results2/network_std_gpsday.out";
    double avg_std;
    std::ifstream inFilet2;
    inFilet2.open (path_to_sig2.c_str());
        if(inFilet2.is_open()){
         //   for(int i = 0; i<num_clocks; i++){
         //       for(int j = 0; j<JW;j++){
                    inFilet2 >> avg_std;
//                    std::cout << ss[i][j] << " ";
               //     signal[i][j] = ss[i][j];
            //    }
//                std::cout << "\n";
           // }
        }
        inFilet2.close();    

    double avg_sig;
    avg_sig = avg_std*avg_std;

    double snr_thresh;
    double xi;
    xi = num_clocks*(b0[0]/avg_sig);
    snr_thresh = 1/(2 + xi*(1-(1/num_clocks))); 
 
    // end of new threshold calculations from SNR-max paper   

 
    for (int i = MAXJW/2+1; i < num_epochs-MAXJW/2; i++){
        
        if(fabs(m_snr[i]) > snr_thresh){
        //if(fabs(m_snr[i]) > 6.0){
    	    chi2 << i << "  " << m_chi[i] << " " <<  num_clocks << " " << chi_thresh  << " " << m_snr[i] <<"\n";
            std::cout << "epoch: " << i << "    " << "snr[i]: " << m_snr[i] << "chi[i]: " << m_chi[i] << "\n";
            fpcount2++;
        }


	if(fabs(m_snr[i]) > snr_thresh && m_chi[i] < chi_thresh){
        //if(fabs(m_snr[i]) > 6.0){
            chi4 << i << "  " << m_chi[i] << " " <<  num_clocks << " " << chi_thresh  << " " << m_snr[i] <<"\n";
            std::cout << "epoch: " << i << "    " << "snr[i]: " << m_snr[i] << "chi[i]: " << m_chi[i] << "\n";
            //fpcount2++;
        }
        
        // full network (ns = 30)
//        if(m_chi[i] > 2215 && fabs(m_snr[i]) > 7.7){
//            m_snr[i] = 0.0;     // let us just set the SNR to zero, no event now.
//        }
        
        // this is for a network of 21 clocks (only Rb)
        if(m_chi[i] > chi_thresh && fabs(m_snr[i]) > 6.8){
  //          m_snr[i] = 0.0;     // let us just set the SNR to zero, no event now.
            fpcount++;
        }
        
    }
    
    std::cout << "this is P.E. count above SNR: " << fpcount2 << "\n";
    std::cout << "this is F.P. count above SNR and Chi: " << fpcount << "\n";
    std::cout << "SNR: " << m_snr[2849] << " " << "chi: " << m_chi[2849] <<"\n";
    
// end of chi-squared filter 09/16/20
    
    chi.close();
    chi2.close();
    test.close();
    chi3.close();
    smatt.close();
    chi4.close();
    
    m_investigate = 0;
if(m_investigate){//param est snr
    std::cout <<"hello in investigate" <<"\n";
    int m_estNMC2 = 11000;
      int event = 1815;
      std::cout << "event: " << event << "\n";
      snrParameterEstimation(model, *this, icov, 1, 0, 10*week+day, event, JW, m_estNMC2,
                          hmin, v_prior, psi_prior, exact, m_Einv_e, m_investigate);
      std::vector<double> best_params;
      likeParameterEstimation( model, *this, icov, 1, 0, 10*week+day, event, JW, m_estNMC2,
          hmin, v_prior, psi_prior, best_params, m_Einv_e, exact);
      
}

printf("Day took %f seconds to complete using %i cores.\n", end-start,ilimit_cores);
    // double totaldiff=(double)(clock() - totalStart)/CLOCKS_PER_SEC;
    // printf("Entire loop over the day took: %.10fs\n", totaldiff);
    
  return 0;

}

//need to put likelihoods shit in here
void SearchClass::serialize(int week, int day, int sanity, int NMC, int jbeg, int jend, int chi_thresh)
{
  
    std::cout << "Hello in serialize" << "\n";
    
  int event  = jbeg+(jend-jbeg)/2;
  std::string wd = std::to_string(week) + std::to_string(day);

  // if(m_investigate == 0){
  //   std::string chiPath = "./results/chi" + std::to_string(week) + std::to_string(day) + ".out";
  //   std::ofstream chi(chiPath);
  //   chi << "Epoch     MCMC it.      dHs     sHs";
  //   for(int i = jbeg; i <= jbeg; i++){
  //     chi << i << dHd[i] << '\n';
  //     for(int j = 0; j < 2*NMC; j++){
  //       chi << i << "   " << j/2 << "   " << m_likelihood[i*NMC + j] 
  //           << "    " << m_likelihood[i*NMC + j + 1] << '\n';
  //           j++;
  //     }
  //   }
  //   chi.close();
  // }


  std::string out_file1 = "./results/"+wd+"/info" + wd + ".out";
  std::ofstream fOutput1(out_file1);

  fOutput1 << "Week: " << week << " Day: " << day;
  
  fOutput1 << "\nReference used: " << refname << "\nSatellites: \n";
  for(int i = 0; i < num_clocks; i++){
    fOutput1 << clk[i] << blk[i] << "-" << prn[i] <<" ";
    //if( (i-num_receivers+1 ) % 10 == 0) fOutput << "\n";
  }



  std::ofstream flag("./results/flag.out",std::ofstream::app);

//  m_investigate = 1;

  std::string out_file;
  if(m_investigate){
    out_file = "./results/"+wd+"/searchData" + wd + "-" + std::to_string(event) + ".out";
  }else{
    //out_file = "./results/"+wd+"/searchData" + wd + ".out";
    out_file = "./results/searchData" + wd + ".out";
  } 
  
	std::ofstream fOutput(out_file);

  //fOutput << "Epoch\tOdds\tLogO\tSNR\n";
  //int count = 0;
  for(int i = jbeg; i <= jend; i++){
    //if( ( fabs(m_snr[i]) > SNR_cut || *(m_odds + i)> odds_cut) ) {
    if( ( fabs(m_snr[i]) > SNR_cut || *(m_odds + i)> odds_cut) && m_chi[i] < chi_thresh) {
      flag << week << day << "\t"<< i << std::endl;
      //count+=1;
    }
    fOutput << i << "\t" << *(m_odds + i) << "\t" << log10(*(m_odds + i))
            << "\t" << m_snr[i] << std::endl;
  }
  flag.close();
	fOutput.close();
  
}
