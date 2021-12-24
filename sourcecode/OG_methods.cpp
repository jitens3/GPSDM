#include "methods.h"
#include <time.h>


/*=============================================================================
========================= Link to LAPACK functions ============================
=============================================================================*/
extern "C" {
  extern void dpotrf_(char *uplo,int *jb,double *A,int *lda,int *info);
  extern void dpotri_(char *uplo,int *jb,double *A,int *lda,int *info);
  extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
  extern void dsymv_(char *uplo, int *n,double *alpha,double *A,int *lda,double *x, int *incx, double *beta, double *Y, int *incy);
  extern void dsymm_(char *side, char *uplo, int *M, int *N, double * alpha, double *A,
                     int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC);
}

/*=============================================================================
============================== MClieklihoods ==================================
=============================================================================*/
int MClikelihoods(
  std::string model,
  JplGpsData data,
  int icov, bool use_ccc, bool sub_w_mean,
  std::vector< std::vector<double> > &odds,
  int jbeg, int jend, int JW,
  int NMC,
  double h0,
  NumericCdfInverse v_prior, 
  NumericCdfInverse psi_prior, 
  double *Einv,
  int exact
)
/*
170519.
Randomised Montr-Carlo approach for the integrations to calculate Odds ratio.
Uses either 'mixed' MC (where the h integral is done normally) or 'Full MC'
(where h integral is done with others in randomised MC grid).
Can also do the h integral analytically.
Incorperates the covariance matrix; for now, assumes no correlation between
clocks.

=======Change log=====
170619- Made fullMC an option
170622- Added "fast" versions of exp/erf/inverseErf functions
      - Made genSignal more efficient (not finished)
170623- Starting to play with analytic h integral
170624- Proper K normalisation (Gaussian likelihood)
170627- Fixed error with v. Too small v means chi_ss becomes zero?? [analyt h]
        *! but now, the v spread is too small?
        --> for now, using flat v distro! (needs testing)
170625- Changed Pno etc to log! (odds not log)
170702- Change s array to static double!
      - calcDeterminant, tempH uses std::vector
170705. I sometimes get errors where chi_ss is zero! (even though s ins't).
        This isn't possible. May be indicative of a very bad error.
        Also, when this happens the s signal doesn't look correct!
        Problem only seems to occur when using openMP !!! [w/ // over its]
170705- Starting to change to "shorter" arrays, for each JW !!
        Removed all static arrays from in here, they don't work!!
        Only 1 static array left (dData), outside, in testLikelihoods
        Seems to work...
        BUT. removed thick wall, monopole signals. Need to re-code!
170710- Incorperating Ian's addition for parameter estimation
        WARNING: this will slow down the program! So, should be removed
        after testing is complete! (in likelihoods.cpp)
170715-  Starting to add thick wall/monopole signals. NOTE: still need to
        include correct "relative" h's etc., and correct rho's (for mono)
170816- Adding the impact parameters for monopole/string signal.
170831- Starting updated to work with class. Major changes.
      - Killed 'calcLikelihoods' (non-MC method)
      - Works. But, lots still to do!
170913- Changed vmin/vmax
171018- Cross-clock-correlation now included (optionally) into likelihoods
171020- Includes proper priors for v and anlges! (inluding 'rotation')
171023- Updated Signal generator. Nicer, and now allows "mixed differencing"
171106- Removed FullMC and MixedMC
171107- Parameter estimation routine.
171115- JplGpsData Class may have non-zero lag Wik terms, but the Wikjl Matrix
        is too large, so it's not used here.
171120- Finishing the "subtract mean" routine for s
180118- Starting to convert from Wik to b0 method. Note: CCC terms same for
        each pair of clocks!
190619- Takes in int exact for option to use numerical 
        inverse of covar matrix instead of approximation.
        Still need to implement routine that marginalizes likelihood
        with the exact inversion XXX
TO-DO
-----
      XXXStill need to implement routine that marginalizes likelihood
      with the exact inversion XXX
*/
{

  //Just to shorten some lines
  int num_epochs = data.num_epochs;
  int num_clocks = data.num_clocks;

  //Copy s.d. and difference arrays into new arrays.
  //This adds a small amount of over-head, and uses a small amount of extra
  //memory, but makes code a little nicer?
  // NB: is it better to store in a c-array? Better cache use?? Test! XXX
  std::vector<int> dif; //differencing level
    dif = data.vdiff;
  std::vector<float> sdev; //standard deviations
    sdev = data.sdev;
  if(data.keff.size()==0) data.formRelativeCouplings(); //for safety
  std::vector<float> keff; //effective (relative) couplings
    keff = data.keff;

  //XXX better to use arrays? or vectors? Some performance testing is needed!
    
  //Work out beginning and end-points of integration
  const int hJW = (JW-1)/2; //nb: JW must be odd! [not huge problem if not]
  //hJW is "half" window - window on either side of j0
  if(jbeg<=hJW) jbeg=hJW+1;//XXX Check if outside data range.
  if(jend>=num_epochs-hJW-1) jend=num_epochs-hJW-2;


  //"Wegithed Mean" normalisation constant
  //(For subtrating the 'mean' from each s)
  double weight_norm=0;
  for(int i=0; i<data.num_clocks; i++){
    weight_norm += 1/data.sdev[i];
  }

  //Generate the Hijl "Hessian" matrix (Inverts the E^i_jl matrix)
  double acf_cut=0.025; //XXX make input!? Needs testing!
  // (any ACF points smaller than acf_cut will be set to zero!)
  int max_cov=5;  //Maximum number of points allowed in ACF [0 after]
  double *Hijl;
  Hijl = new double[num_clocks*JW*JW];

  if(icov==1){
    if(exact == 0){
      genHijl(data,Hijl,JW,max_cov,acf_cut);
      //Gaussian likelihood log-normalisation constant log(K)
    }
  }

  double logNorm = likelihoodLogNorm(icov,num_clocks,JW,sdev,Hijl);

  //Generate the (0-lag) cross-clock correlation part of Hessian
  double b0=0;
  if(use_ccc) b0 = data.b0[0];

  //===================================================
  // Loop over epochs:
  //#pragma omp parallel for
  for(int j0=jbeg; j0<=jend; j0++){

    //Define small data array, just centered around relevant window
    //Means we have to throw less data around, make better use of cache
    double dD[MAXCLOCKS][MAXJW];//careful, trailing junk!
    double rsat[MAXCLOCKS][3];
    double rref[3];
    int jmin;
    prepareData(hJW,j0,jmin,dD,rsat,rref,data);

    double chi_dd=0.;
    if(exact == 1){
      //do exact stuff
      //chi_dd = calcChiExact(dD, dD, Einv, data.num_clocks, JW);
      //BLAS is a faster exact method
      chi_dd = calcChi_BLAS(Einv, *dD, *dD, data.num_clocks, JW);
    }
    else{
      //Calculate the 'd.d' part of Chi^2
      chi_dd = calcChi(dD,dD,Hijl,sdev,JW,num_clocks,icov);
      //Include the cross-clock-correlation term:
      if(use_ccc) chi_dd += calcChiW(dD,dD,b0,sdev,JW,num_clocks);
    }

    //The likelihood functions:
    double Pdm = 0; //holds positive likelihood
    double logPno = (logNorm-0.5*(chi_dd));//negative likelihood.

    // Monte Carlo integration:
    //#pragma omp parallel for reduction(+:Pdm)
    //Need to use reduction(+:Pdm) ??? TEST!
    //NB: no need to parallelise this loop! Just adds over-head!
    for (int ix=0; ix<NMC; ix++){

      //Randomly choose integration params
      double t0,v,d,R,a;  //time, vel, width, impact param, impact angle
      double n[3];        //incident direction
      randomParameters(v_prior,psi_prior,model,j0,t0,v,n,d,R,a);

      //generate signal:
      double s[MAXCLOCKS][MAXJW];//holds signal
      generateSignal(model,s,rsat,rref,num_clocks,jmin,JW,dif,keff,
        t0,v,n,d,R,a);

      //Subtract the weighted mean
      // s -> s - s_bar
      if(sub_w_mean) weightedMeanSignal(s,sdev,weight_norm,num_clocks,JW);

      double chi_ds=0., chi_ss=0.;
      if(exact == 1){
        //do exact
        //BLAS is a faster exact method
        chi_ds = calcChi_BLAS(Einv, *s, *dD, data.num_clocks, JW);
        chi_ss = calcChi_BLAS(Einv, *s, *s, data.num_clocks,  JW);
        //chi_ds = calcChiExact(s, dD, Einv, data.num_clocks, JW);
        //chi_ss = calcChiExact(s, s, Einv, data.num_clocks,  JW);
      }
      else{
        //Calculate chi_ds=d.s, and chi_ss=s.s. [chi_dd done above!]
        // NB: always call with s first! (more efficient)
        chi_ds = calcChi(s,dD,Hijl,sdev,JW,num_clocks,icov);
        chi_ss = calcChi(s,s,Hijl,sdev,JW,num_clocks,icov);
        // Cross-clock Correlations:
        if(use_ccc && !exact){
          chi_ds += calcChiW(s,dD,b0,sdev,JW,num_clocks);
          chi_ss += calcChiW(s,s,b0,sdev,JW,num_clocks);
        }
      }

      if(chi_ss==0){
        //XXX Stand-in. ALSO: "other" error??? XXX
        std::cout<<"\nERROR 345 in xx: ss=0?\n";
        continue;
      }

      //Do the h-integral analytically:
      double tempPdm = hAnalytic(chi_ds,chi_ss,h0);
      Pdm += tempPdm;

    }//END NMC loop [Monte Carlo]

    //Normalise (integration normalisation)
    Pdm/=NMC;

    //Odds ratio:
    odds[0][j0]=Pdm; //Odds ratio
    //Log-likelihood:
    odds[1][j0]=log(Pdm)+logPno; //Pno=MFS_fastExp(logNorm-0.5*chi_dd)

  }//END loop over j0

  //==========================================
  //Output results to file

  //Prints odds value for each epoch
  std::ofstream epochOdds;
  std::string oddsData = "oddsData.txt";
  epochOdds.open(oddsData, std::ios_base::app);
    for(int i=jbeg;i<jend;i++){
      epochOdds<<odds[0][i]<<"\t";
      epochOdds<<odds[1][i]<<"\n";
    }
  epochOdds<<std::endl;
  epochOdds.close();

  //Prints max odds from window
  //std::ofstream epochOdds;
  //std::string oddsData = "oddsData.txt";
  //epochOdds.open(oddsData, std::ios_base::app);
  //    epochOdds<<max_o<<std::endl;
  //epochOdds.close();

  return 0;
}



/*=============================================================================
=========================== likeParamEstimation ===============================
=============================================================================*/
int likeParameterEstimation(
  std::string model,
  JplGpsData data,
  int icov, bool use_ccc, bool sub_w_mean,
  int itn_or_day,
  int j0, int JW,
  int NMC,
  double h0,
  NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  std::vector<double> &best_params, double *Einv, int exact
)
/*
171107.
New routine that does the parameter estimation.
Doesn't do integration, just finds maximum of likelihood.
Runs for single epoch only
Has some repitition with MClikelihoods... ok?

190619- Takes in int exact and double * Einv for option to use numerical 
        inverse of covar matrix instead of approximation.
TO-DO
-----
      XXX Add method if exact inversion option chosen!!!
*/
{

  //hJW is "half" window - window on either side of j0
  const int hJW = (JW-1)/2; //nb: JW must be odd! [not huge problem if not]

  //Generate the Hijl "Hessian" matrix (Inverts the E^i_jl matrix)
  double acf_cut=0.025; //XXX make input!? Needs testing!
  // (any ACF points smaller than acf_cut will be set to zero!)
  int max_cov=61;  //Maximum number of points allowed in ACF [0 after]
  double *Hijl;
  Hijl = new double[data.num_clocks*JW*JW];
  if(icov==1){
    if (exact==0){
      genHijl(data,Hijl,JW,max_cov,acf_cut);
    }
  }

  //Generate the cross-clock correlation part of Hessian [0 lag only]
  double b0=0;
  if(use_ccc) b0 = data.b0[0];

  //"Wegithed Mean" normalisation constant
  //(For subtrating the 'mean' from each s)
  double weight_norm=0;
  for(int i=0; i<data.num_clocks; i++){
    weight_norm += 1/data.sdev[i];
  }

  //Define small data array, just centered around relevant window
  //Means we have to throw less data around, make better use of cache
  double dD[MAXCLOCKS][MAXJW];//careful, trailing junk!
  double rsat[MAXCLOCKS][3];
  double rref[3];
  int jmin;
  prepareData(hJW,j0,jmin,dD,rsat,rref,data);

  //Best fit parameters
  double best_t0=-1,best_v=0,best_th=-1,best_ph=-1,best_d=0,best_R=-1,best_a=-1
  ,best_h=0;
  double best_O=-1;

  // Monte Carlo integration:
  //#pragma omp parallel for reduction(max:best_O) //XXX what about others!????
  //NB: no point parallelising this loop. Just adds overhead.
  for (int ix=0; ix<NMC; ix++){

    //Randomly choose integration params
    double t0,v,d,R,a;  //time, vel, width, impact param, impact angle
    double n[3];        //incident direction
    randomParameters(v_prior,psi_prior,model,j0,t0,v,n,d,R,a);

    //generate signal:
    double s[MAXCLOCKS][MAXJW];//holds signal
    generateSignal(model,s,rsat,rref,data.num_clocks,jmin,JW,
                   data.vdiff,data.keff,t0,v,n,d,R,a);

    //Subtract the weighted mean
    // s -> s - s_bar
    if(sub_w_mean)
      weightedMeanSignal(s,data.sdev,weight_norm,data.num_clocks,JW);

    double chi_ds=0., chi_ss=0.;
    if(exact == 1){
      //BLAS is a faster exact method
      chi_ds = calcChi_BLAS(Einv, *s, *dD, data.num_clocks, JW);
      chi_ss = calcChi_BLAS(Einv, *s, *s, data.num_clocks,  JW);
      //chi_ds = calcChiExact(s, dD, Einv, data.num_clocks, JW);
      //chi_ss = calcChiExact(s, s, Einv, data.num_clocks,  JW);
    }
    else{
      //Calculate chi_ds=d.s, and chi_ss=s.s. [chi_dd done above!]
      // NB: always call with s first! (more efficient)
      chi_ds=calcChi(s,dD,Hijl,data.sdev,JW,data.num_clocks,icov);
      chi_ss=calcChi(s,s,Hijl,data.sdev,JW,data.num_clocks,icov);
      // Cross-clock Correlations:
      if(use_ccc){
        chi_ds += calcChiW(s,dD,b0,data.sdev,JW,data.num_clocks);
        chi_ss += calcChiW(s,s,b0,data.sdev,JW,data.num_clocks);
      }
    }
    //Just for safety??
    if(chi_ss==0) continue;

    //Do the h-integral analytically:
    double tempPdm = hAnalytic(chi_ds,chi_ss,h0);

    //Perform the "maximisation"
    if(tempPdm>best_O){
      //update best-fit parameters
      best_O  = tempPdm;
      best_t0 = t0;
      best_v  = v;
      best_th = acos(n[2]);
      best_ph = atan2(n[1],n[0]);
      best_d  = d;
      best_R  = R;
      best_a  = a;
      best_h  = chi_ds/sqrt(chi_ss);
    }

  }//END NMC loop [Monte Carlo]

  //Add the parameters to output array
  best_params.clear(); //needed? just tp be safe..
  best_params.push_back(best_t0);
  best_params.push_back(best_O);
  best_params.push_back(best_h);
  best_params.push_back(best_v/30.);
  best_params.push_back(best_th/M_PI);
  double tmp_phi=best_ph;
  if(tmp_phi<0) tmp_phi += 2*M_PI;
  best_params.push_back(tmp_phi/M_PI);
  best_params.push_back(best_d);
  best_params.push_back(best_R);
  best_params.push_back(best_a/M_PI);

  //===============================================
  //parameter estimation output
  std::ofstream paramFile;
  std::string eventParams = "./results/"+std::to_string(itn_or_day)+"/eventParams-odds"+std::to_string(itn_or_day)+".out";
  paramFile<<"t0\todds\th\tv\ttheta\tphi\td\tR\ta\n";
  paramFile.open(eventParams, std::ios_base::app);
    for(int i_prm=0; i_prm<9;i_prm++){//fix!need->10 for mono, etc.
      paramFile<<best_params[i_prm]<<"\t";
      //if(i_prm==8) paramFile<<std::endl;    
    }
  paramFile<<"\n";
  paramFile.close();

  return 0;
}




/*=============================================================================
=================================== MCsnr =====================================
=============================================================================*/
int MCsnr(std::string model, JplGpsData data, int icov, 
	bool use_ccc, bool sub_w_mean, double *snr,
	int jbeg, int jend, int jw,	int NMC, double h0,
	NumericCdfInverse v_prior, NumericCdfInverse psi_prior, 
  int exact, double * Einv) 
/*
190619- Takes in int exact and double * Einv for option to use numerical 
        inverse of covar matrix instead of approximation.
TO-DO
-----
*/
{
  //*****************************************
	//Copy s.d. and difference arrays into new arrays.
	//This adds a small amount of over-head, and uses a small amount of extra
	//memory, but makes code a little nicer?
	// NB: is it better to store in a c-array? Better cache use?? Test! XXX
	std::vector<int> dif; //differencing level
	dif = data.vdiff;
	std::vector<float> sdev; //standard deviations
	sdev = data.sdev;
	if (data.keff.size() == 0) data.formRelativeCouplings(); //for safety
	std::vector<float> keff; //effective (relative) couplings
	keff = data.keff;
  
  //=============================================
  //Work out beginning and end-points of integration
	const int hJW = (jw - 1) / 2; //nb: JW must be odd! [not huge problem if not]
	//hJW is "half" window - window on either side of j0
	if (jbeg <= hJW) jbeg = hJW + 1;//XXX Check if outside data range.
	if (jend >= data.num_epochs - hJW - 1) jend = data.num_epochs - hJW - 2;

  //Generate the Hijl "Hessian" matrix (Inverts the E^i_jl matrix)
  double acf_cut=0.025; //XXX make input!? Needs testing!
  // (any ACF points smaller than acf_cut will be set to zero!)
  int max_cov=5;  //Maximum number of points allowed in ACF [0 after]
  double *Hijl;
  Hijl = new double[data.num_clocks*jw*jw];
  if(icov==1){
    if (exact==0){
      genHijl(data,Hijl,jw,max_cov,acf_cut);
    }
  }

	//Generate the (0-lag) cross-clock correlation part of Hessian
	double b0 = 0;
	if (use_ccc) b0 = data.b0[0];

  //"Wegithed Mean" normalisation constant
  //(For subtrating the 'mean' from each s)
  double weight_norm=0;
  for(int i=0; i<data.num_clocks; i++){
    weight_norm += 1/data.sdev[i];
  }

	double chi_val[data.num_epochs] = {0};

	double max_snr = 0.; //set orig max snr to 0

  //XXX used for printing snr values w/o injected signal
  //Print SNR for each template
  //make this an option? nah.
  //std::ofstream templateSnr;
  //templateSnr.open("./results/snrData.txt", std::ios_base::app);

	//================================================
	// Loop over epochs:
	//#pragma omp parallel for
	for (int j0 = jbeg; j0 <= jend; j0++) {

		//Define small data array, just centered around relevant window
		//Means we have to throw less data around, make better use of cache
		double dD[MAXCLOCKS][MAXJW];//careful, trailing junk!
		double rsat[MAXCLOCKS][3];
		double rref[3];
		int jmin;
		prepareData(hJW, j0, jmin, dD, rsat, rref, data);

		snr[j0]=0.; //make sure 1st comparison of calculated snr is to 0

		double chi_dd=0.;
        	chi_dd = calcChi(dD,dD,Hijl,sdev,jw,data.num_clocks,icov);
        //Include the cross-clock-correlation term:
        	if(use_ccc) chi_dd += calcChiW(dD,dD,b0,sdev,jw,data.num_clocks);

		//Monte Carlo Parameter Generation:
		for (int ix = 0; ix < NMC; ix++) {

			//Randomly choose params
			double t0, v, d, R, a;  //time, vel, width, impact param, impact angle
			double n[3];        //incident direction
			randomParameters(v_prior, psi_prior, model, j0, t0, v, n, d, R, a);

			//generate signal:
			double s[MAXCLOCKS][MAXJW];//holds signal
			generateSignal(model, s, rsat, rref, data.num_clocks, jmin, jw, data.vdiff, data.keff,
				t0, v, n, d, R, a);

			//Subtract the weighted mean
			// s -> s - s_bar
			if (sub_w_mean) weightedMeanSignal(s, sdev, weight_norm, data.num_clocks, jw);

      double dHs=0., sHs=0., snrVal=0.;
      if(exact == 1){
        //
        //dHs=calcChiExact(s, dD, Einv, data.num_clocks, jw);
        //sHs=calcChiExact(s, s, Einv, data.num_clocks, jw);
        
        dHs= calcChi_BLAS(Einv, *s, *dD, data.num_clocks, jw);
        sHs= calcChi_BLAS(Einv, *s, *s, data.num_clocks, jw);
      }
      else{
//          std::cout << "this is data.num_clocks[i]" << data.num_clocks << "\n";
			 // NB: always call with s first! (more efficient)
			 //Calculate dHs and sHs
			 dHs = calcChi(s, dD, Hijl, data.sdev, jw, data.num_clocks, icov);
			 sHs = calcChi(s, s, Hijl, data.sdev, jw, data.num_clocks, icov);

			 //Include the cross-clock-correlation term:
			 if (use_ccc) {
			 	dHs += calcChiW(s, dD, b0, data.sdev, jw, data.num_clocks);
			 	sHs += calcChiW(s, dD, b0, data.sdev, jw, data.num_clocks);
			 }
      }
			if (sHs == 0) {
				//XXX Stand-in. ALSO: "other" error??? XXX
				std::cout << "\nERROR 345 in xx: ss=0?\n";
				continue;
			}

			// new 06/23/20
            double h;
            double chiVal;
            h = dHs/sHs;
            chiVal  = chi_dd -2*h*dHs + (h*h)*sHs;
            // end of new 06/23/20

			//Calculate SNR
			snrVal = dHs / sqrt(sHs);     
      //snrVal1 = dHs1 / sqrt(sHs1); 
      
      //XXX print to std::ofstream object above for loop over j0
      //templateSnr << snrVal << std::endl;//print template snr

			if(fabs(snrVal) > fabs(snr[j0])){
				snr[j0] = snrVal;
				chi_val[j0] = chiVal;
			}
			if(fabs(snr[j0])>fabs(max_snr)){
				max_snr = snr[j0];
			}
			//Parameter estimation and best snr check
		//MSC_progressBar("SNR iterations", 30, ix, NMC,1);
		}//END NMC loop [Monte Carlo]

	
	}//END loop over j0

  //XXXclose std::ostream object above for loop over j0
  //templateSnr.close();//end template snr printing
	

  //===============================================
  //Output results to file

	//Prints SNR value for each epoch
  std::ofstream epochSnr;
  std::string snrData = "./results2/snr4Data.txt";
  epochSnr.open(snrData, std::ios_base::app);
  	for(int i_snr=jbeg;i_snr<jend;i_snr++){
  		epochSnr<<snr[i_snr]<<"\n";
  	}
  epochSnr<<std::endl;
  epochSnr.close();

    std::ofstream epochSnrt;
    std::string snrDatat = "./results/snrtData.txt";
    epochSnrt.open(snrDatat, std::ios_base::app);
      for(int i_snr=jbeg;i_snr<jend;i_snr++){
          epochSnrt<<snr[i_snr]<<"  " << chi_val[i_snr]<<"\n";
      }
    epochSnrt<<std::endl;
    epochSnrt.close();
	//Prints max snr from window
	//std::ofstream epochSnr;
	//std::string snrData = "./results/snrData.txt";
	//epochSnr.open(snrData, std::ios_base::app);
	//		epochSnr<<max_snr<<std::endl;
	//epochSnr.close();

  //Prints snr at injected epoch
  //std::ofstream epochSnr;
  //std::string snrData = "./results/snrData.txt";
  //epochSnr.open(snrData, std::ios_base::app);
  //    epochSnr<<snr[126]<<std::endl;
  //epochSnr.close();

	return 0;
}


/*=============================================================================
========================== snrParameterEstimation =============================
=============================================================================*/
int snrParameterEstimation(
  std::string model,
  JplGpsData data,
  int icov, bool use_ccc, bool sub_w_mean,
  int itn_or_day,
  int j0, int JW,
  int NMC,
  double h0,
  NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
  int exact, double * Einv, int invest)
{

  double dsArray[NMC]={0};
  double ssArray[NMC]={0};
  double bestsig[data.num_clocks][JW]={0};
  //hJW is "half" window - window on either side of j0
  const int hJW = (JW-1)/2; //nb: JW must be odd! [not huge problem if not]

  //Generate the Hijl "Hessian" matrix (Inverts the E^i_jl matrix)
  double acf_cut=0.025; //XXX make input!? Needs testing!
  // (any ACF points smaller than acf_cut will be set to zero!)
  int max_cov=61;  //Maximum number of points allowed in ACF [0 after]
  double *Hijl;
  Hijl = new double[data.num_clocks*JW*JW];
  if(icov==1){
    if (exact==0){
      genHijl(data,Hijl,JW,max_cov,acf_cut);
    }
  }

  //Generate the cross-clock correlation part of Hessian [0 lag only]
  double b0=0;
  if(use_ccc) b0 = data.b0[0];

  //"Wegithed Mean" normalisation constant
  //(For subtrating the 'mean' from each s)
  double weight_norm=0;
  for(int i=0; i<data.num_clocks; i++){
    weight_norm += 1/data.sdev[i];
  }

  //Define small data array, just centered around relevant window
  //Means we have to throw less data around, make better use of cache
  double dD[MAXCLOCKS][MAXJW];//careful, trailing junk!
  double rsat[MAXCLOCKS][3];
  double rref[3];
  int jmin;
  prepareData(hJW,j0,jmin,dD,rsat,rref,data);
    
    
    // this makes a vector of d^T*E^-1 11/1/20 tday
    double vec[data.num_clocks*MAXJW];
    for(int i=0;i<data.num_clocks*JW;i++){
          vec[i] = 0;
      }
    calcdE(dD, Einv, vec, JW, data.num_clocks);      // this function makes the vector
    // end of making new vector on 11/17/20 tday

  // Monte Carlo integration:
  //#pragma omp parallel for reduction(max:best_O) //XXX what about others!????
  //NB: no point parallelising this loop. Just adds overhead.
  double bestParams[11]={0};
    double ng[3]; //predicted DM impact vector
    double phi_lower; //angle between predicted DM direction "ng[]" and direction from data "n[]"
    double ng_dot_n;
    double mag_ng_n;
	//{1,	2,	3,	4,	5,	6,	7,	8,	9}
	//{t0,snr,h,	v,	th,	phi,	n[1],	n[2],	n[3]}
//    std::cout << "Hello Tyler this is num clocks: " << data.num_clocks << "\n";
//    std::cout << "Hello Tyler this is maxclocks: " << MAXCLOCKS << "\n";
    std::cout << "Hello Tyler this is NMC in param est: " << NMC << "\n";
    
  for (int ix=0; ix<NMC; ix++){

    //Randomly choose integration params
    double t0,v,d,R,a;  //time, vel, width, impact param, impact angle
    double n[3];        //incident direction
    randomParameters(v_prior,psi_prior,model,j0,t0,v,n,d,R,a);

    //generate signal:
    double s[MAXCLOCKS][MAXJW];//holds signal
    generateSignal(model,s,rsat,rref,data.num_clocks,jmin,JW,
                   data.vdiff,data.keff,t0,v,n,d,R,a);

    //Subtract the weighted mean
    // s -> s - s_bar
    if(sub_w_mean)
      weightedMeanSignal(s,data.sdev,weight_norm,data.num_clocks,JW);

    double dHs=0., sHs=0.;
    if(exact == 1){
      //dHs = calcChiExact(s, dD, Einv, data.num_clocks, JW);
      dHs = calcdEs(vec, s, JW, data.num_clocks);
      //sHs = calcChiExact(s, s, Einv, data.num_clocks, JW);
      sHs = calcChi_BLAS(Einv, *s, *s, data.num_clocks, JW); // added 11/17/20
    }
    else{
      //Calculate chi_ds=d.s, and chi_ss=s.s. [chi_dd done above!]
      // NB: always call with s first! (more efficient)
      dHs = calcChi(s,dD,Hijl,data.sdev,JW,data.num_clocks,icov);
      sHs = calcChi(s,s,Hijl,data.sdev,JW,data.num_clocks,icov);
      // Cross-clock Correlations:
      if(use_ccc){
        dHs += calcChiW(s,dD,b0,data.sdev,JW,data.num_clocks);
        sHs += calcChiW(s,s,b0,data.sdev,JW,data.num_clocks);
      }
    }

    //Just for safety??
    if(sHs==0) continue;

    dsArray[ix]=dHs;
    ssArray[ix]=sHs;
    double tempSnr=dHs/sqrt(sHs);
    double tempsnrsd=1/sqrt(sHs);   // added 12/8/20 to store the h-hat sd.
    //Perform the "maximisation"
    //update best-fit parameters
    if(fabs(tempSnr) > fabs(bestParams[1])){
	    bestParams[0]=t0;
      bestParams[1]=tempSnr;
	    bestParams[2]=dHs/sHs;
	    bestParams[3]=v/30.;
	    bestParams[4]=acos(n[2])/M_PI;
      double tmp_phi=atan2(n[1],n[0]);
      if(tmp_phi<0) tmp_phi += 2*M_PI;
      bestParams[5]=tmp_phi/M_PI;
	    bestParams[6]=n[0];
	    bestParams[7]=n[1];
	    bestParams[8]=n[2];
        
        //checking components of n vectors against predicted DM direction ng
        ng[0] = 0.46;
        ng[1] = -0.49;
        ng[2] = 0.74;

        ng_dot_n = ng[0]*n[0] + ng[1]*n[1] + ng[2]*n[2];
        mag_ng_n = sqrt(pow(ng[0],2)+pow(ng[1],2)+pow(ng[2],2)) * sqrt(pow(n[0],2)+pow(n[1],2)+pow(n[2],2));

        phi_lower = acos(ng_dot_n/mag_ng_n);
        
        bestParams[9] = phi_lower;
        // end of directionality angle updated on 11/20/20
        bestParams[10] = tempsnrsd;
	  
      for(int a=0; a<data.num_clocks;a++){
		    for(int b=0; b<JW;b++){
			    bestsig[a][b] = s[a][b];
		    }
	    }
    }
    
  }//END NMC loop [Monte Carlo]

  //===============================================
  //if using search.cpp
  //Ouput parameters to file
  //std::ofstream paramFile;
    if(itn_or_day <= 12000){ //tests if in testMethods program or search program
      std::string eventParams = "./results/eventParams-snr.out";
      std::ofstream paramFile;
      paramFile.open(eventParams, std::ios_base::app);
      if(itn_or_day == 0){
        paramFile<<"itn\tt0\tSNR\th\tv\ttheta\tphi\tvector\n";
      }
        paramFile<<itn_or_day<<"\t";
        for(int i_prm=0; i_prm<11;i_prm++){//fix!need->10 for mono, etc.
          paramFile<<bestParams[i_prm]<<"\t";
          if(i_prm==10) paramFile<<std::endl;
        }
        paramFile.close();

    }else{
      std::string eventParams = "./results/"+std::to_string(itn_or_day)+"/eventParams-snr-"+std::to_string(itn_or_day)+".out";
      std::ofstream paramFile;
      paramFile.open(eventParams, std::ios_base::app);
    //  if(invest){
        paramFile<<"t0\tSNR\th\tv\ttheta\tphi\tvector\n";
        for(int i_prm=0; i_prm<11;i_prm++){//fix!need->10 for mono, etc.
          paramFile<<bestParams[i_prm]<<"\t";
          if(i_prm==10) paramFile<<std::endl;
        }
        paramFile<<"\n";
  
        std::ofstream sigFile("./results/"+std::to_string(itn_or_day)+"/signalTemplate"+std::to_string(itn_or_day)+"-"+std::to_string(j0)+".out");
        if(invest == 1){
          for(int a=0;a<data.num_clocks;a++){
   //         if(data.clk[a]=="Rb"){
              for(int b=0;b<JW;b++){
              sigFile<<bestsig[a][b]<<" ";
              }
              sigFile<<std::endl;
   //         }
          }
          sigFile.close();
        }
        else{
          for(int a=0;a<data.num_clocks;a++){
              for(int b=0;b<JW;b++){
              sigFile<<bestsig[a][b]<<" ";
              }
              sigFile<<std::endl;
          }
          sigFile.close();
        }

        std::ofstream chiFile;
        std::string chiVals = "./results/"+std::to_string(itn_or_day)+"/chi-vals-snr-"+std::to_string(itn_or_day)+"-"+std::to_string(j0)+".out";
        chiFile.open(chiVals, std::ios_base::app);
        for(int i=0; i<NMC;i++){//fix!need->10 for mono, etc.
            chiFile<<dsArray[i]<<"\t"<<ssArray[i]<<"\n";
            //if(i_prm==8) paramFile<<std::endl;    
        }
        chiFile.close();
    //  }
      paramFile.close();
    }

  return 0;
}




//******************************************************************************
double likelihoodLogNorm(int icov, int Nclk, int JW, std::vector<float> sdev,
                         double * Hijl)
/*
171106.
Gaussian likelihood log-normalisation constant log(K)
See "Normalisation and numerical stability" in Bayes memo.
nb: doesn't include cross-clock correlation part..
*/
{
  double mylogDet=0;
  if(icov==1){
    //For covariance case only: [uses det(H)]
    for(int i=0; i<Nclk; i++){
      std::vector< std::vector<double> > tempH(JW,std::vector<double>(JW));
      for(int l=0; l<JW; l++){
        for(int m=0; m<JW; m++) tempH[l][m]=*(Hijl +i*JW*JW + l*JW + m); //[i][l][m];
      }
      double tempDet=MFS_calcDeterminant(tempH);
      if(tempDet!=0) mylogDet+=JW*log(sdev[i])-0.5*log(tempDet);
      if(tempDet==0) mylogDet+=JW*log(sdev[i]); //??
    }
  }else{
    for(int i=0; i<Nclk; i++) mylogDet+=JW*log(sdev[i]);
  }
  return -0.5*Nclk*JW*log(2.*PI)-mylogDet;
}


//******************************************************************************
void prepareData(int hJW, int j0, int &jmin,
                 double dD[MAXCLOCKS][MAXJW],
                 double rsat[MAXCLOCKS][3], double rref[3],
                 JplGpsData data)
/*
171106.
Define small data array, just centered around relevant window.
This helps with performance...maybe? XXX
window: [j0-hJW, j0+hJW]
*/
{
    //std::cout << "Hello from prepareData" << "\n";
    
  jmin=j0-hJW;
  int JW=2*hJW+1; //temp

  //Define small data array, just centered around relevant window
  //Means we have to throw less data around, make better use of cache
    //std::cout << "num clocks: " << data.num_clocks << "\n";       // checks out
  for(int i=0; i<data.num_clocks; i++){
    for(int j=0; j<JW; j++){
      dD[i][j]=data.bias[i][j+jmin];
    }
  }
  //Form  rsat and rref arrays:
  //Here, I just take the "slice" of positions that I need to calculate
  // the signals.
  for(int i=0; i<data.num_clocks; i++){
    for(int ix=0; ix<3; ix++)
      rsat[i][ix]=data.pos[i][j0][ix];
  }
  for(int ix=0; ix<3; ix++)
    rref[ix]=data.refpos[j0][ix];

}

//******************************************************************************
void prepareDataSearch(int hJW, int j0, int &jmin,
                 double dD[MAXCLOCKS][MAXJW],
                 double rsat[MAXCLOCKS][3], double rref[3],
                 JplGpsData data)
/*
190806.
Changes from prepareData to work on real Jpl data. Only care about
the satellite clocks!

Define small data array, just centered around relevant window.
This helps with performance...maybe? XXX
window: [j0-hJW, j0+hJW]
*/
{
  jmin=j0-hJW;
  int JW=2*hJW+1; //temp

  //Define small data array, just centered around relevant window
  //Means we have to throw less data around, make better use of cache
  for(int i=data.num_receivers; i<data.num_clocks; i++){
    for(int j=0; j<JW; j++){
      dD[i-data.num_receivers][j]=data.bias[i][j+jmin];
    }
  }
  //Form  rsat and rref arrays:
  //Here, I just take the "slice" of positions that I need to calculate
  // the signals.
  for(int i=data.num_receivers; i<data.num_clocks; i++){
    for(int ix=0; ix<3; ix++)
      rsat[i-data.num_receivers][ix]=data.pos[i][j0][ix];
  }
  for(int ix=0; ix<3; ix++)
    rref[ix]=data.refpos[j0][ix];

}


//******************************************************************************
void randomParameters(NumericCdfInverse v_prior, NumericCdfInverse psi_prior,
                        std::string model, int j0,
                        double &t0, double &v, double n[3],
                        double &d, double &R, double &a)
/*
171106.
Randomly chooses the integration params.
Uses given numerical priors for v and angles (n).
Logarithmic prior for d.
Uses geometric priors for R, a
Flat t0 = (j0-1, j0)
*/
{

  // choose t0 (collission time)
  t0 = MFS_randDouble(double(j0-1),double(j0));

  //Choose the incident |v| using the numeric priors
  {
    double uv=MFS_randDouble(0.001,1.);
    v = v_prior.inverseCdf(uv)*30.; //double check x30! tau_0
  }

  // Angles (direction) eDM = n
  //Note: we choose the angles in the frame that has z along n_gal
  // this makes the priors very easy. Then, we just rotate back to
  // the ECI frame.
  {
    //Generate angles in DM (psi/omega) frame:
    double u_psi = MFS_randDouble(0.,1.);
    double psi   = psi_prior.inverseCdf(u_psi);
    double omega = MFS_randDouble(0.0001,2.*PI);//omega has flat prior
    //Rotate back into ECI (theta/phi) frame:
    double cp   = cos(psi);
    double cosp = cos(omega)*sin(psi);
    double spso = sin(psi)*sin(omega);
    n[0] = -0.46332*cp + 0.179476*cosp + 0.867827*spso; //x
    n[1] =  0.49003*cp + 0.867827*cosp + 0.082144*spso; //y
    n[2] = -0.73838*cp + 0.463320*cosp - 0.490030*spso; //z
  }

  //For thick walls, monopoles, strings: d, R, and alpha
  // d: object width
  // R: "Overall" impact parameter, R
  // a: "Overall" incident angle (in plane per. to n), alpha
  if(model!="twall"){
    //Choose d on a log grid. (use Inverse transform sampling)
    double dmin = 1.e2;    // ?? CHECK! XXX input!?
    double dmax = 1.e6; // ?? CHECK!    XXX input!?
    double ud=MFS_randDouble(0.,1.);
    d=dmin*pow((dmax/dmin),ud);
    if(model!="dwall"){//string + monopole
      //choose alpha on a ragular grid
      //Choose R on a "power" grid P(R)~R. [0->Rgps]
      a=MFS_randDouble(0.,2.*PI);
      double ur=MFS_randDouble(0.,1.);
      R=RGPS*sqrt(ur);
    }
  }

}


//******************************************************************************
double hAnalytic(double ds, double ss, double h0)
/*
Do the h-integral analytically:
Note: I take h_max -> infinity in the integral, but take it to be some
finite constant to work out normalisation 9arbitrarily chosen).
Also: I include a h_min in the integral, which is h0. Can be zero.
*/
{
  double hmax_norm=0.1; //Only for normalisation!!
  double kh=1/(2.*hmax_norm);
  double tempA=MFS_fastExp(pow(ds,2)/(2*ss));
  double tempss=sqrt(2*ss);
  double tempB=kh*1.77245/tempss; //NB: Sqrt(Pi) = 1.77245
  double tempC = 2;
  // "full" formula for 'H_min' part:
  tempC += MFS_fastErf( (ds - h0*ss)/tempss )
          -MFS_fastErf( (ds + h0*ss)/tempss );
  //take hmax->infty in int (but not norm)!
  return tempA*tempB*tempC;
}


//******************************************************************************
double calcChi_Exact_fast(double *E_inverse, double *X, 
  double *Y, int num_clocks, int JW)
{

  double chi = 0;
    for(int iclock = 0; iclock < num_clocks; iclock++){
     for(int iwindow = 0; iwindow < JW; iwindow++){
      double temp2 = 0;
      for(int jclock = 0; jclock < num_clocks; jclock++){
        for(int jwindow = 0; jwindow < JW; jwindow++){
          if (*(X + jclock*num_clocks + jwindow)==0) continue; //just don't bother with these calculations
          else{
            temp2 = 0;
            double hijl= *(E_inverse + iclock*JW*JW*num_clocks + iwindow*num_clocks*JW + jclock*JW + jwindow);
          if(hijl==0)break; //just don't bother with these calculations
          else{
            temp2+=hijl * *(X + jclock*num_clocks + jwindow);
            // std::cout << temp2 << " clock row & column: " << iclock << "    " << jclock 
            //           << " window row and column: " << iwindow << "   " << jwindow << "\n";
            }
          }
         }//end of jwindow
        }//end of jclock
        if(*(Y + iclock*num_clocks + iwindow) == 0){
          continue; //just don't bother with these calculations
        }else{
          chi += temp2* *(Y + iclock*num_clocks + iwindow);
        }
       }//end of iwindow
      }//end of iclock




  //   for(int i=0;i<num_clocks;i++){
  //     double temp=0;
  //     for(int j=0;j<JW;j++){
  //       if(*(X + i*num_clocks + j)==0){
  //           continue;
  //       }
  //       double temp2=0;
  //       temp2 += *(E_inverse + i*JW*JW + j*JW + j)*d2[i][j];
  //       for(int l=j-1; l>=0; l--){
  //         double hijl=*(hessian + i*JW*JW + j*JW + l);
  //         if(hijl==0)break;
  //         temp2+=hijl*d2[i][l];
  //       }
  //       for(int l=j+1; l<JW; l++){
  //         double hijl=*(hessian + i*JW*JW + j*JW + l);
  //         if(hijl==0)break;
  //         temp2+=hijl*d2[i][l];
  //       }//END loop over epochs l
  //       temp+=d1[i][j]*temp2;
  //     }//END loop over epochs j
  //     chi_xx+=temp/pow(sig,2);
  //   }//END loop over clocks.
  // }//END IF covariance
  return chi;
}


//******************************************************************************
double calcChi(
  double d1[MAXCLOCKS][MAXJW], //input "first" vector (d or s)
  double d2[MAXCLOCKS][MAXJW], //input "second" vector (d or s)
  double *hessian, //Input Hessian matrix
  //double sdev[MAXCLOCKS], //input clock standard deviation
  std::vector<float> sdev, //input clock standard deviation
  //int jmin, int jmax,
  int JW,
  int iClocks, int icov)
/*
170511.
Calculates single term of the argument of the likelihood function:
-Does NOT include the negative, or the (1/2) factor!
-Does include the 'sigma^2' factor though.
Note: putting the "if" check for covariance outside the loops makes the code a
little more chunky, but makes it a little but faster (haven't tested though).

NOTE: Function will be faster when called with (s,d), instead of (d,s)!!!
This is due to the fact that if |d1|<'cut', then: skip!
Only works because of "acf_cut" in genHijl function.

==== Change Log ====
170705- Change to "shorter version" of array
171012- Made more efficient. After H->0, can stop! [but +/-]
        Only works because of "acf_cut" in genHijl function.
*/
{
  //'Cut-off' for d; below this, will be =~0. Speeds up greatly!
  //double d_cut = 0.01; //in units of 'sigma'
  //XXX XXX XXX XXX NOTE: this likely NOT valid for monopoles/thick walls!!

    //std::cout << "Hello from calcChi 1" << "\n";
    
    //std::cout << "number of clocks: " << iClocks << "\n";
    
  double chi_xx=0;
  if(icov==0){
    // No covariance:
    for(int i=0;i<iClocks;i++){
      double sig=sdev[i];
      double temp=0;
      for(int j=0;j<JW;j++){
        temp+=d1[i][j]*d2[i][j];
      }
      chi_xx+=temp/pow(sig,2);
    }
  }else{
    //Multivariate case (with covariance):
    for(int i=0;i<iClocks;i++){
        //std::cout << "Hello from calcChi 1.5" << "\n";
      double sig=sdev[i];
      double temp=0;
      for(int j=0;j<JW;j++){
        //if(fabs(d1[i][j])<d_cut*sig)continue;// esp if d=s, mostly 0!
        // XXX HERE! XXX
        if(d1[i][j]==0)continue; //faster for thin walls
        double temp2=0;
        //First, calculate the zero-lag term (l=j)
        temp2+=*(hessian + i*JW*JW + j*JW + j)*d2[i][j];
        //Calculate the step "backwards" terms (l<j)
        for(int l=j-1; l>=0; l--){
          double hijl=*(hessian + i*JW*JW + j*JW + l);
          if(hijl==0)break; //reach effective "end" of H matrix
          temp2+=hijl*d2[i][l];
        }//END "backwards" loop over epochs l
        //Calculate the step "forwards" terms (l>j)
        for(int l=j+1; l<JW; l++){
          double hijl=*(hessian + i*JW*JW + j*JW + l);
          if(hijl==0)break; //reach effective "end" of H matrix
          temp2+=hijl*d2[i][l];
        }//END "forwards" loop over epochs l
        temp+=d1[i][j]*temp2;
      }//END loop over epochs j
      chi_xx+=temp/pow(sig,2);
    }//END loop over clocks.
  }//END IF covariance
    //std::cout << "Hello from calcChi 2" << "\n";
  return chi_xx;
}

//******************************************************************************
double calcChiW(
  double d1[MAXCLOCKS][MAXJW], //input "first" vector (d or s)
  double d2[MAXCLOCKS][MAXJW], //input "second" vector (d or s)
  double b0, //Input cross-clock correlation term
  //double sdev[MAXCLOCKS], //input clock standard deviation
  std::vector<float> sdev, //input clock standard deviation
  int JW,
  int num_clocks)
/*
171017.
Function that calculates the d1*W*d2 [cross-clock corr.] part of chi^2
-Does NOT include the negative from Gaussian, or the (1/2) factor!
--(Does include the ccc negative)
-Does include the 'sigmai^2 sigmak^2' factor though.
Assumes no cross-clock correlations for lag>0.
XXX In the case where d1 = d2 [sWs term], can be 2x faster!!(?)
  ===== Change Log  =====
180118- Updated to use b0 method, instead of Wik. This is because CCC is the
        same for each pair of clocks
*/
{

  if(b0==0) return 0;

  double chi_W_xx=0;
  // Calculate the d1*W*d2 term:
  for(int i=0; i<num_clocks; i++){
    for(int k=0; k<num_clocks; k++){
      if(k==i) continue; // k != i in this term! //is there a better way?
      double t_d1d2=0;
      for(int j=0; j<JW; j++){
        t_d1d2 += d1[i][j]*d2[k][j];
      }//j
      chi_W_xx += t_d1d2/pow(sdev[i]*sdev[k],2);
    }//k
  }//i

  return -1*b0*chi_W_xx;
}

//******************************************************************************
double calcdE(
              double d1[MAXCLOCKS][MAXJW],       // this is the data vector
              double *m_Einv_e,     // this is the hessian matrix
              double dEij[MAXCLOCKS*MAXJW],         // this is the saved vector variable
              int JW,
              int num_clocks)
/*
 this function is to compute the d^T * (E^-1) multiplication prior to the MC loop
 this is done in order to save some computation time when using the exact inverse method
 this function returns the variable dEij, or maybe it should return zero and the saved vector
 is saved in a variable passed to the function?
 11/25/19 tday
 */
{
    double* dtran;
    int size = num_clocks*JW;
    //std::cout << "size: " << size << "\n";
    //dEij = new double [size];
    dtran = new double [size];
    double dD[size] = {0};
    
    // initilaize the transpose of the data
    for(int i=0; i<size; i++){
        *( dtran + i ) = 0;
    }
    
    int k = 0;
    for(int i=0; i<num_clocks; i++){
        for(int j=0; j<JW; j++){
            dD[k] = d1[i][j];
            k++;
        }
    }
    
    
    // lets make the transpose of data vector
    for(int i=0; i<size; i++){
        *(dtran + (i) ) = dD[i];
    }
    
    
    // now lets try to multiply the matrices together dtran and m_Einv_e
    for(int i=0; i<size; i++){
        for(int j=0;j<size;j++){
            //*(dEij + i) += (*(dtran + j )) * (*(m_Einv_e + j*size + i));                  // added 11/25/19
            dEij[i] += (*(dtran + j )) * (*(m_Einv_e + j*size + i));                  // added 11/25/19
        }
    }
    
    
    //delete [] dEij;
    delete [] dtran;
    
    return 0;
    
}

//******************************************************************************
double calcdEs(
               double dEij[MAXCLOCKS*MAXJW],       // this is the saved vector
               double d2[MAXCLOCKS][MAXJW],         // this is the saved vector variable
               int JW,
               int num_clocks)
/*
 this function is to multiply the vector (d^T*E^-1) with the signal vector
 this is for the exact inverse method and takes the place of calcChiBLAS function for that
 particular case of dHs.
 11/26/19 by tday
 */
{
    double chi_xx = 0;
    
    int k = 0;
    // now lets multiply the two vectors together
    for(int i=0; i<num_clocks; i++){
        for(int j=0; j<JW; j++){
            chi_xx += dEij[k]*d2[i][j];
            k++;
        }
    }
    
    
    return chi_xx;
}

//******************************************************************************
int genHijl(
  JplGpsData data,
  double *hessian, //output Hessian matrix
  int JW,
  int max_cov,
  double acf_cut //ACF/Hessian cut-off threshold
  )
/*
170512.
Takes in the ACFs, and generates the Hessian matrices.
Has a cutoff, 'acf_cut'. Any value lower than this will be set to zero. Applies
both the the input ACF/Eij matrix, and the output Hij matrix. Can set this to
zero to have no cut-off.
===== Change Log =====
170702- Uses std::vector for Ejl matrix
171013- Updated to use ACF from JplGpsData class.
*/
{
    std::cout << "hello inside genHijl 1" << "\n";
    std::cout << "this is max_cov: " << max_cov << "\n";
    
  // added 09/22/20 tday to check if getting correct inverse
  double Hijl[MAXCLOCKS][MAXJW][MAXJW]={0};     // added on 11/22/19 by tday
  std::string atACFpath = "./results2/new61_abtACF13101.out";
  std::ofstream atACF_file(atACFpath);
    
  std::string aitACFpath = "./results2/new61_aiACF13101.out";
  std::ofstream aitACF_file(aitACFpath);
  // end of new additions on 09/22/20 tday
    
// check on the acf values
    for(int i=0;i<data.num_clocks;i++){
        for(int j=0;j<max_cov;j++){
            std::cout << "this is acf[i][j]: " << data.acf[i][j] << "\n";
        }
        std::cout << "\n";
    }
    
    

  int H_test = 0;
  if(JW>MAXJW) JW=MAXJW;// Just for safety?

  if(H_test == 0){
      
      // added 08/20/20
      int n;
      n = 0;
      //

  for(int i=0; i<data.num_clocks; i++){
 //std::cout << "line 1159 \n";
    //temp E_ij cov. matrix. To be inverted.
    std::vector< std::vector<double> > Ejl(JW,std::vector<double>(JW)); //this might be slow
    // std::cout << "line 1162 \n";
    // Form the "modified" Covariance. [modified since s.d. factored out]
    for(int j=0; j<JW; j++){
      for(int l=0; l<JW; l++){
        double eval=0;    //temp. hold value for this Eij value.
        int tau=abs(j-l); //"lag"
        if(tau<max_cov) {
            
            
          eval=data.acf[i][tau];
          //eval=data.acf[n][tau];        // added on 08/20/20
          //std::cout << "ACF to be put into H: " << eval << std::endl;
          } ///XXX check for over-flow!
        if(fabs(eval)<acf_cut) eval=0;
        Ejl[j][l]=eval;
          Hijl[i][j][l] = Ejl[j][l];          // added on 09/22/20
          atACF_file << (data.sdev[i]*data.sdev[i])*Hijl[i][j][l] << "  ";    // added on 09/22/20
      }
        atACF_file << "\n";    // added on 09/22/20
    }
 //     std::cout << "hello inside genHijl 2" << "\n";
   //Invert E_ij
   int iOK=MFS_invertMatrix(Ejl);
   //std::cout << "iok\n";
   if(iOK!=0)std::cout<<"FAILURE 477 in genHijl: Eij matrix ("<<i
                      <<") could not be inverted!\n";
   for(int j=0;j<JW;j++){
      for(int l=0;l<JW;l++){
        *(hessian + i*JW*JW + j*JW + l)=Ejl[j][l];
        if(*(hessian + i*JW*JW + j*JW + l)<acf_cut){
          *(hessian + i*JW*JW + j*JW + l)=0; //OK?
            Ejl[j][l] = 0.;                                // added on 09/22/20
        }
        //std::cout << j << "  " << l << std::endl;
          Hijl[i][j][l] = Ejl[j][l];          // added on 09/22/20
          aitACF_file << Hijl[i][j][l] << "  ";    // added on 09/22/20
      }
       aitACF_file << "\n";    // added on 09/22/20
   }
      n++;      // added on 08/20/20
 //std::cout << "line 1184 \n" << "loop ctr: " << i << std::endl;
 }//END loop over clocks
  }else{
  //Alternate for testing.
  for (int i = 0; i < data.num_clocks; i++){
    for (int j = 0; j < JW; j++){
      for(int l = 0; l < JW; l++){
        *(hessian + i*JW*JW + j*JW + l) = 0;
      }
      *(hessian + i*JW*JW + j*JW + j) = 1;
    }
  }
  }
 //   std::cout << "hello inside genHijl 3" << "\n";

  return 0;
}



//******************************************************************************
int calcEinv( double * Einv, JplGpsData & data, int jw, int icov,
      double max_cov, double acf_cut, int day, int week, double ref_sig, 
      int by_hand_ccc, int sanity)//for now adding day and week for file.
/* FUNCTION: calcEinv
*   This function fill the upper triangular part of the covariance matrix
* (autocorr and cross-corr) and passes it to the LAPACK cholesky decomposition
* and cholesky inversion functions.
*INPUTS:  double * Einv - pointer to memory where space is allocated for E/Einv
*     JplGpsData & data - class object by reference
*     int jw - time window 
*     double max_cov - max lag for covar
*     double acf_cut - cutoff for autocorr before set to zero
*     double ref_sig - reference clock standard dev
*
*OUTPUTS: double *  - pointer to location of inverse of E
*
*NOTE: Does not compute ccc from simulated data, instead fills ccc block matrices
*       with diagonal elements equal to the ref clock variance.
*
*CHANGELOG:
*190619 - NEW FUNCTION!
*/
{
//  std::cout << "\n =============== computing E matrix! =============\n\n"
//        << "Given " << data.num_clocks <<" clocks"
//        << " and " << jw <<" epochs." << std::endl;

  int size = data.num_clocks*jw;
  //int size = (data.num_clocks - 1)*jw;       // added on 03/23/20
    
// added 03/22/20 for pos. def. cov matrix test
  std::vector< std::vector<double> > oEij(size,std::vector<double>(size));
    

    //Initialize E elements to zero
    for(int i=0; i<size; i++){
      for(int j=0;j<size;j++){
        *( Einv + i*size + j ) = 0;
        oEij[i][j] = 0;         // added 03/22/20
      }
    }
    //=====================================================================
    //=================== Form auto-corr part of E  =======================
    //=====================================================================
    if(icov==1){
      for(int i=0; i<data.num_clocks; i++){
         // if(i==2){
         //     continue;
         // }
        double vari = data.sdev[i]*data.sdev[i];
        for(int j=0; j<jw; j++){
            for(int l=0; l<jw; l++){
              float eval=0;    //temp. hold value for this Eij value.
              int tau = abs(j-l); //lag
              if(tau<max_cov) eval=data.acf[i][tau]; ///XXX check for over-flow!
              if(fabs(eval)<acf_cut) eval=0;
                if(l>=j){
                  *( Einv + (i*jw+j)*size + (i*jw+l) ) = eval*vari;
                  *( Einv + (i*jw+l)*size + (i*jw+j) ) = eval*vari;
                  //*(Einv+j*size+l)=arrayE[j][l]; //conversion
                }else{
                  //*( Einv + (i*jw+j)*size + (i*jw+l) ) = 0;
                }
              }
          }
      }
    }
    else{
      for(int i=0; i<data.num_clocks; i++){
          double vari = data.sdev[i]*data.sdev[i];
          for(int j=0; j<jw; j++){
              for(int l=0; l<jw; l++){
                if(l==j){
                    *( Einv + (i*jw+j)*size + (i*jw+l) ) = vari;
                    //*(Einv+j*size+l)=arrayE[j][l]; //conversion
                }else{
                    //arrayE[(i*jw+j)][i*jw+l] = 0;
                    *( Einv + (i*jw+j)*size + (i*jw+l) ) = 0;
                }
              }
          }
      }
    }

    
  //======================================================================
  //=================== Form cross-corr part of E   ======================
  //======================================================================

    if(by_hand_ccc){
    //------------------------------------------------------
    //-------------- "by-hand" ccc calculation -------------
    //------------------------------------------------------
      
    // testing to make C.C.C. lag 01/21/21 tdaykin
        for(int a=0; a<data.num_clocks; a++){
        //    if(a==2){
        //        continue;
        //    }
          for(int b=a+1; b<data.num_clocks; b++){
              for(int j=0; j<jw; j++){
                  
                  for(int l=0; l<jw; l++){
                      float eval2=0;    //temp. hold value for this Eij value.
                      int tau2 = abs(j-l); //lag
                      
                      if(tau2 < max_cov) eval2 = data.b0[tau2];
                      if(l>=j){
                       
                          *( Einv + (a*jw+l)*size + (b*jw+j) ) = eval2;
                          *( Einv + (a*jw+j)*size + (b*jw+l) ) = eval2;
                          *( Einv + (b*jw+l)*size + (a*jw+j) ) = eval2;
                          *( Einv + (b*jw+j)*size + (a*jw+l) ) = eval2;
                          
                      }
                      
                  }
                  
                  //*(Einv+j*size+l)=arrayE[j][l]; conversion
                }
            }
        }
    // end of testing to make C.C.C. lag 01/21/21 tdaykin
        
      //XXX NEEDS TESTING! XXX
//      for(int a=0; a<data.num_clocks; a++){
      //    if(a==2){
      //        continue;
      //    }
//        for(int b=a+1; b<data.num_clocks; b++){
//            for(int j=0; j<jw; j++){
//                *( Einv + (a*jw+j)*size + (b*jw+j) ) += data.b0[0];
//                *( Einv + (b*jw+j)*size + (a*jw+j) ) += data.b0[0];
                //*(Einv+j*size+l)=arrayE[j][l]; conversion
//              }
//          }
//      }
    }
    else{
    //------------------------------------------------------
    //------------------- simulation ccc -------------------
    //------------------------------------------------------
      //for white noise simulations, just fill ccc with ref clock variance
      for(int a=0; a<data.num_clocks; a++){
        for(int b=a+1; b<data.num_clocks; b++){
            for(int j=0; j<jw; j++){
                *( Einv + (a*jw+j)*size + (b*jw+j) ) = ref_sig*ref_sig;
                *( Einv + (b*jw+j)*size + (a*jw+j) ) = ref_sig*ref_sig;
                //*(Einv+j*size+l)=arrayE[j][l]; conversion
            }
        }
      }
  }

  std::cout<<"b0 = "<<data.b0[0]<<" sx^2 = "<<ref_sig*ref_sig<<std::endl;

// added on 03/22/20 tday for check of pos. def
// cov matrix
    
    for(int a=0; a<size; a++){
        for(int b=0; b<size; b++){
            oEij[a][b] = *(Einv + a*size +b);
            //if(abs(oEij[a][b])<1e-5){
            //    oEij[a][b] = 0.0;
            //}
        }
    }
    
    std::string ocACFpath = "./results/oc3Eij" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream ocACF_file(ocACFpath);
    //
    
    
    for(int a=0; a<size; a++){
        for(int b=0; b<size; b++){
            
            ocACF_file << oEij[a][b] << "  ";
        }
        ocACF_file << "\n";
    }
    
    ocACF_file.close();
// end of new addition on 03/22/20 tday
    

double * Einv_copy;
Einv_copy = new double[size*size];

if(sanity == 1){
std::string Epath = "../results/E" + std::to_string(week) + std::to_string(day) + ".out";
std::ofstream Eout(Epath);  

//making a copy
    for(int a=0; a<size; a++){
      for(int b=0; b<size; b++){
          *( Einv_copy + a*size + b ) = *( Einv + a*size + b );
          Eout << *( Einv + a*size + b ) << ' ';
         }
    }

Eout.close();

}else{
  delete [] Einv_copy;
}


  //========================================================================
  //==================  Matrix inversion below  ==========================
  //======================================================================== 
           
    //Invert E
    int info;
    clock_t tStart=clock();

    //Computes the Cholesky factorization
    //takes in upper sizeXsize upper triangular part of E
    //returns upper cholesky factor U
    char UPLO = 'U';
    //char LPLO = 'L';
    int LDA = size;
    
    // this is to do eithr LU or Cholesky 12/9/19 tday
    int tcheck = 0;         // this is for LU=1, or chol = 0
    if(tcheck ==1){
        MFS_invertMatrix(oEij);
    }else{
    dpotrf_(&UPLO,&size,Einv,&LDA,&info);
        
    // added on 01/04/21 tday
        if(info!=0){
          std::cout<<"FAILURE 477 in genHijl: Eij matrix could not be inverted!\n";
          if(info<0){
            std::cout<<"i = " << info << " illegal value!" << std::endl;
          }
          else{
            std::cout<<"(" << info <<"," <<info << ") element is zero!" << std::endl;
          }
            // need a hard stop here
            std::exit(EXIT_FAILURE);
        }
        // end of added on 01/04/21 tday
    //Computes the inverse
    //takes in cholesky factor U
    //returns lower triangular part of Einv
    dpotri_(&UPLO,&size,Einv,&size,&info);
    }
    // end of new additions


//    dpotrf_(&UPLO,&size,Einv,&LDA,&info);
    //Computes the inverse
    //takes in cholesky factor U
    //returns lower triangular part of Einv
//    dpotri_(&UPLO,&size,Einv,&size,&info);
    


    double tdiff=(double)(clock() - tStart)/CLOCKS_PER_SEC;
    printf("Matrix inversion took: %.10fs\n", tdiff);

    if(info!=0){
      std::cout<<"FAILURE 477 in genHijl: Eij matrix could not be inverted!\n";
      if(info<0){
        std::cout<<"i = " << info << " illegal value!" << std::endl;
      }
      else{
        std::cout<<"(" << info <<"," <<info << ") element is zero!" << std::endl;
      }
        // need a hard stop here
        std::exit(EXIT_FAILURE);
    }

//    //dpotri makes lower triangle elements of E inv
//    //need to set upper!
//    for(int a=0; a<size; a++){
//      for(int b=a+1; b<size; b++){
//          *( Einv + a*size + b ) = *( Einv + b*size + a );
//         }
//    }
    
    
    //dpotri makes lower triangle elements of E inv
         //need to set upper!
         if(tcheck==0){
             for(int a=0; a<size; a++){
               for(int b=a+1; b<size; b++){
                   *( Einv + a*size + b ) = *( Einv + b*size + a );
                  }
             }
         }
      //   for(int a=0; a<size; a++){
      //     for(int b=a+1; b<size; b++){
      //         *( Einv + a*size + b ) = *( Einv + b*size + a );
      //        }
      //   }
         
         // added on 12/9/19 by tday to do LU decomp
         if(tcheck==1){
             for(int a=0; a<size; a++){
                 for(int b=0; b<size; b++){
                     if(abs(oEij[a][b])<acf_cut){
                         oEij[a][b]=0;
                     }
                     *(Einv + a*size + b) = oEij[a][b];
                 }
             }
         }else{
             for(int a=0; a<size; a++){
                 for(int b=0; b<size; b++){
                     oEij[a][b] = *(Einv + a*size + b);
                 }
             }
         }
         // end of new additions on 12/9/19 tday

    std::string Einv_path = "./results/E3inv" + std::to_string(week) + std::to_string(day) + ".out";
    std::ofstream Einvout(Einv_path);
        for(int a=0; a<size; a++){
          for(int b=0; b<size; b++){
              //Einvout << *( Einv + a*size + b ) << ' ';
              Einvout << oEij[a][b] << ' ';
             }
            Einvout << "\n";
        }
    Einvout.close();

if (sanity == 1){
  std::string Einv_path = "../results/Einv" + std::to_string(week) + std::to_string(day) + ".out";
  std::ofstream Einvout(Einv_path);  
      for(int a=0; a<size; a++){
        for(int b=0; b<size; b++){
            Einvout << *( Einv + a*size + b ) << ' ';
           }
      }
  Einvout.close();
      double *identity;
      identity = new double[size*size];
      for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
          *(identity + i*size + j) = 0;
        }
      }
      char side = 'l';
      double alpha = 1;
      double beta = 0;
      int lda = size;

      dsymm_(&side, &UPLO, &size, &size, &alpha, Einv, &lda, Einv_copy, &lda, &beta, identity, &lda);

  int countdiag = 0;
  int countoffdiag = 0;
  int excess = 0;
      for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
          if(*(identity + i*size + j) <= 0.000001 && i != j){
            countoffdiag++;
          }
          if(*(identity + i*size + j) <= 1.000001 && *(identity + i*size + j) >= 0.999999 && i == j){
            countdiag++;
          }else if(i == j && excess <= 5){
            std::cout << "uh - oh " << *(identity + i*size + j) << '\n';
            excess++;
          }else if(excess == 6){
            std::cout << "... supressing all future 'uh-oh's ...\n";
            excess = 7;
          }
        }
      }
      if(countdiag == size && countoffdiag == size*size - size){std::cout << "\nIdentity achieved!\n\n";}
      else {std::cout << "\nIdentity not achieved\n\n";}

      delete [] Einv_copy;
      delete [] identity;
}

  return 0;
}


//******************************************************************************
int calcEinvWhiteNoise( double * Einv, JplGpsData & data, int jw, double sigma, double sigma_x){

  sigma = 0.1;
  //sigma_x = 0.006;
  int size = data.num_clocks*jw;
  double squiggle = data.num_clocks*(sigma_x*sigma_x)/(sigma*sigma);
  double cross = (squiggle/(1+squiggle)); 

    //Initialize E elements to zero
    for(int i=0; i<size; i++){
      for(int j=0;j<size;j++){
       
         *( Einv + i*size + j ) = 0;
       
        
      }
    }
    for(int i = 0; i < data.num_clocks; i++){
      for(int j = 0; j < data.num_clocks; j++){
        for(int l = 0; l < jw; l++){
          for (int k = 0; k < jw; k++){
            if (i == j && l == k){
              *(Einv + i*jw*size + j*jw + l*size + k) = (1/(sigma*sigma))*(1-cross/data.num_clocks);
            }else if(l == k && i != j){
              *(Einv + i*jw*size + j*jw + l*size + k) = -1*(1/(sigma*sigma))*cross/data.num_clocks;
            }
          }
        }
      }
    }
    
  return 0;
}

//******************************************************************************
double calcChiExact(double x[MAXCLOCKS][MAXJW], 
  double y[MAXCLOCKS][MAXJW], 
  double * Einv, int numClks, int jw)
/*
*NEW FUCNTION
* Computes (x^T)(E^-1)(y) exactly.
*
*Change Log -
*190619 - 
*
*To Do
*-----
*   XXX Optimize multiplications!
*/
{

  int size = numClks*jw;

    double chi_xy=0.;

    //set ptr to first element of x/y[i][j]
    double * yprime = y[0];
    double * xprime = x[0];
std::cout << "line 1424 in calcchiexact\n";
    //calc Ey
    double temp[size]; //hold (E^-1)*(y)
    double * tempPtr = temp;
    for(int i=0;i<size;i++){
      *(tempPtr+i)= 0;
      for(int j=0;j<size;j++){
        if(*(xprime+j) == 0)continue; //faster for thin walls
        *(tempPtr+i) += *(Einv+i*size+j) * (*(xprime+j));
        //temp[i] += Einv[i][j] * yprime[j];

      }//END loop over epochs j
    }//END loop over clocks.
std::cout << "line 1437 in calcchiexact\n";
    //calc x(Ey)
    for(int i=0;i<size;i++){
      
      //for(int j=0;j<data.num_clocks*jw;j++){
        if(*(yprime+i)==0)continue;//faster for thin walls

        chi_xy += (*(yprime+i))*( *(tempPtr+i) );
      //}//END loop over epochs j
    }//END loop over clocks.
  std::cout << "line 1447 in calcchiexact\n";
 return chi_xy;

}

double calcChi_BLAS(double *E_inverse, double *X, 
  double *Y, int num_clocks, int jw)
/*
*NEW FUCNTION
* Computes (x^T)(E^-1)(y) exactly.
*
*Change Log -
*190806 - Implimented into methods code
*
*To Do
*-----
*   Maybe use 4tran functions instead?
*/
{
    //CBLAS_ORDER layout = CblasColMajor;
    char Uplo = 'u';
    double alpha = 1.0;
    double beta = 0.0;
    double chi;
    int size = (num_clocks)*(jw);
    int incx = 1;
    int incy = 1;
    double result[num_clocks*jw] = {0};

    //===========================  Multiplication  =================================
    //auto start = high_resolution_clock::now();
  dsymv_(&Uplo, &size, &alpha, E_inverse, &size, Y, &incx, &beta, result, &incy);

  chi = ddot_(&size, X, &incx, result, &incy);

  return chi;
}
  


//******************************************************************************
int weightedMeanSignal(double s[MAXCLOCKS][MAXJW], std::vector<float> sdev,
                       double W, int num_clocks, int JW)
/*
171017.
Forms a regular c-array for the cross-clock correlation Hessian matrix, Wik.
Only the zero-lag part

XXX CARFEUL XXX THIS DOESN"T WORK FOR MIXED DIFFERENCE!!!!! XXX XXX XXX XXX

XXX Should use variance (sigma^2) instead of sigma for the weights??
XXX NB: make sure EXACTLY the same as

*/
{

  for(int j=0; j<JW; j++){
    //a) find weighted mean of all clocks.
    //For each clock, subtract wegithed mean (of all others!)
    double sbar=0;
    for(int i=0; i<num_clocks; i++){
      sbar += s[i][j]/sdev[i];
    }
    sbar/=W;
    for(int i=0; i<num_clocks; i++){
      //s[i][j] -= sbar;
      //Subtract the weighted mean of all OTHER clocks:
      s[i][j] -= sbar + (sbar - s[i][j])/(W*sdev[i]-1);
      //XXX DOUBLE CHECK! XXX
    }
  }
  return 0;
}

//******************************************************************************
//                           Signal Generator
//******************************************************************************

//******************************************************************************
// int generateSignal(std::string model, double s[MAXCLOCKS][MAXJW],
//            double rsat[MAXCLOCKS][3], double rref[3],
//            int num_clocks, int jmin, int JW, int dif[MAXCLOCKS],
//            double t0, double v, double eDM[3], double d, double R, double a)
int generateSignal(std::string model, double s[MAXCLOCKS][MAXJW],
           double rsat[MAXCLOCKS][3], double rref[3],
           int num_clocks, int jmin, int JW,
           std::vector<int> dif, std::vector<float> keff,
           double t0, double v, double eDM[3], double d, double R, double a)
/*
171023. New function.
Generates the signal for wall/monopole/string.
Uses "mixed" differencing (i.e., Rb can be 1-order, Cs 2-order etc.)
Uses effective/relative coupling strengths (keff).
Uses "small" s array, goes from 0 to JW, so need to shift ti and tr by the
value of jmin!
Makes use of "MFS_fastErf" function, that is accurate enough for our purposes
and much faster than the cmath erf function.

INPUTS:
  model     ::  TD type: twall, dwall, monopole, string
  rsat      ::  Projection of sat's along incident vector n (xi in notes)
  rref      ::  Projection of ref along incident vector n (xi in notes)
  num_clocks::  Number of clocks
  jmin      ::  Offset: start position of signal window in epochs
  JW        ::  Length of signal window
  diff      ::  Array that holds the differenc level (=0,1,2) for each clock*1
  keff      ::  Array that holds the effective/relative couplings
  t0        ::  Time event crosses ECI0
  v         ::  Scalar DM velocity. nb: v=|v|>0
  eDM       ::  Incindent DM direction [=\hat n]
  d         ::  TD width (root-mean-square width). Optional*2
  R         ::  "Overall impact parameter" [perp. distance from ECI0]. Optional
  a         ::  Impact angle (in plane perp. to n). Optional
INPUTS/OUTPUTS:
  s       ::  2D Array that holds the signal for each clock

*1. Can be called with single int here instead (overloaded)
*2. All optional parameters are 0 by default.

*/
{
  //Projection of clocks onto eDM:
  double xi_ref=projectClock(rref[0],rref[1],rref[2],eDM);
  double xi_i=0;  //work out inside loop

  //RELATIVE clock jump-height parameters
  //these depend on which coupling we consider, and which clock type
  double rel_hr=keff[num_clocks];
  double rel_hi=1; //worked out inside loop


  double tr=t0-jmin-(xi_ref/v);  //time object centre hits the reference clock
  double ti=0; //time object centre hits the ith clock (done in loop)
  //Note: subtract jmin to shift offset because of short array

  //Thin wall signal much faster, so do this first:
  if(model=="twall"){
    thinWallSignal(s,rsat,num_clocks,jmin,JW,dif,keff,t0,v,eDM,tr,rel_hr);
    // std::cout << "signal matrix" << std::endl;
    // for (int i = 0; i < num_clocks; i++){
    //   for (int j = 0; j < JW; j++){
    //     std::cout << s[i][j] << " ";
    //   }
    //   std::cout << "\n";
    // }
    return 0;
  }

  // If "thick wall" model, but d is small, take a short-cut:
  double min_d=100.;
  if((model=="dwall")&&(d<min_d)){
    thinWallSignal(s,rsat,num_clocks,jmin,JW,dif,keff,t0,v,eDM,tr,rel_hr);
    // std::cout << "signal matrix" << std::endl;
    // for (int i = 0; i < num_clocks; i++){
    //   for (int j = 0; j < JW; j++){
    //     std::cout << s[i][j] << " ";
    //   }
    //   std::cout << "\n";
    // }
    return 0;
  }

  //Use "full" signal.

  // The "impact factor" (actually Exp[-(q/d)^2], where q is impact factor):
  // Just returns 1 in case of walls.
  double expqd2r=expqd2(model,rref[0],rref[1],rref[2],xi_ref,eDM,d,a,R);
  double expqd2i=1.;
  //form this ratio often:
  double vd = v/d;

  for(int i=0; i<num_clocks; i++){
    xi_i=projectClock(rsat[i][0],rsat[i][1],rsat[i][2],eDM);
    ti=t0-jmin-xi_i/v; //time object centre hits the clock
    expqd2i=expqd2(model,rsat[i][0],rsat[i][1],rsat[i][2],xi_i,eDM,d,a,R);
    rel_hi=keff[i];
    if(dif[i]==1){
      for(int j=0;j<JW;j++){
        s[i][j]=firstOrderSignal(rel_hi,rel_hr,j,ti,tr,vd,expqd2i,expqd2r);
      }
    }else if(dif[i]==2){
      for(int j=0;j<JW;j++){
        s[i][j]=secondOrderSignal(rel_hi,rel_hr,j,ti,tr,vd,expqd2i,expqd2r);
      }
    }else{
      for(int j=0;j<JW;j++){
        s[i][j]=zerothOrderSignal(rel_hi,rel_hr,j,ti,tr,vd,expqd2i,expqd2r);
      }
    }
  }//End loop over clocks
  // std::cout << "signal matrix" << std::endl;
  // for (int i = 0; i < num_clocks; i++){
  //   for (int j = 0; j < JW; j++){
  //     std::cout << s[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  return 0;
}

//------------------------------------------------------------------------------
int generateSignal(std::string model, double s[MAXCLOCKS][MAXJW],
           double rsat[MAXCLOCKS][3], double rref[3],
           int num_clocks, int jmin, int JW, int id, std::vector<float> keff,
           double t0, double v, double eDM[3], double d, double R, double a)
/*
171023.
Overloaded function.
So it can be called with a constant 'difference' level (id) instead of array.
*/
{
 // int dif[MAXCLOCKS];
 std::vector<int> dif(num_clocks);
 for(int i=0; i<num_clocks; i++){
   dif[i]=id;
 }
 generateSignal(model,s,rsat,rref,num_clocks,jmin,JW,dif,keff,t0,v,eDM,d,R,a);
 return 0;
}


//------------------------------------------------------------------------------
int thinWallSignal(double s[MAXCLOCKS][MAXJW], double rsat[MAXCLOCKS][3],
                   int num_clocks, int jmin, int JW, std::vector<int> dif,
                   std::vector<float> keff,
                   double t0, double v, double eDM[3],
                   double tr, double rel_hr)
/*
171023.
Generates the signal in the case of a thin wall (Delta function profile)
Uses "relative"/effective couplings, keff.

=======Change log=====
*/
{
    
    std::string bchiPath = "./results2/signal.out";
    std::ofstream smat(bchiPath);

  //Define terms; worked out inside loops:
  double xi_i=0;  //Projection of clocks onto eDM:
  double rel_hi=1; //RELATIVE clock jump-height parameter
  double ti=0; //time object centre hits the ith clock (done in loop)
    
    /* initialize random seed 11/25/19 */
    unsigned int global_seed = 1;
    
    //std::cout << "num_clocks: " << num_clocks << "\n";

  for(int i=0; i<num_clocks; i++){
      
  //    std::cout << "vdiff[i]: " << dif[i] << "\n";
    //Project down clock (onto n), find t_i given t_0 and v
    xi_i=projectClock(rsat[i][0],rsat[i][1],rsat[i][2],eDM);
    ti=t0-jmin-xi_i/v; //time object centre hits the clock
    rel_hi=keff[i];
      //std::cout << "rel_hi: " << rel_hi << "\n";
    if(dif[i]==1){
      for(int j=0;j<JW;j++){
        s[i][j]=0.0;//array isn't zero'd outside. More efficient this way.
        if((j-1<ti)&&(ti<=j))s[i][j]+=rel_hi;        // commented out 09/24/20 for tests
        if((j-1<tr)&&(tr<=j))s[i][j]-=rel_hr;        // commented out 09/24/20 for tests
        /* this is where i make semi-random signal 10/11/19 */
//        s[i][j] = rand_r(&global_seed) % 2;               // check 10/11/19 (tday)
       //   if(s[i][i] > 2){
       //       s[i][j] = 0.0;
       //   }
  //        std::cout << "signal: " << s[i][j] << "\n";
      }
    }else if(dif[i]==2){
      for(int j=0;j<JW;j++){
        s[i][j]=0;
        if((j-1<ti)&&(ti<=j))        {s[i][j]+=rel_hi;}
        else if((j-2<ti)&&(ti<=j-1)) {s[i][j]-=rel_hi;}
        if((j-1<tr)&&(tr<=j))        {s[i][j]-=rel_hr;}
        else if((j-2<tr)&&(tr<=j-1)) {s[i][j]+=rel_hr;}
      }
    }else{
      //0-order differencing
      for(int j=0;j<JW;j++){
        if((j<ti)&&(j<tr)){s[i][j]=0.;}
        else if((j>=ti)&&(j<tr))  {s[i][j]=rel_hi;}
        else if((j<ti)&&(j>=tr))  {s[i][j]=-1*rel_hr;}
        else if((j>=ti)&&(j>=tr)) {s[i][j]=rel_hi-rel_hr;}
      }
    }//End which difference?
  }//End loop over clocks
    
 /*
    for(int i=0; i<num_clocks; i++){
          for(int j=0; j<JW; j++){
              smat << s[i][j] << " ";
           }
              smat << "\n";
    }
    
    smat.close();
   */

  return 0;
}

//------------------------------------------------------------------------------
double firstOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r)
/*
171021.
Little function that generates the general first-order differenced signal
Make inline??
*/
{
  return 0.5*(
    hi*(
      MFS_fastErf((j-ti)*vd) - MFS_fastErf((j-ti-1)*vd)
    )*expqd2i
    +
    hr*(
      MFS_fastErf((j-tr-1)*vd) - MFS_fastErf((j-tr)*vd)
    )*expqd2r
  );
}
//------------------------------------------------------------------------------
double secondOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r)
/*
171021.
Little function that generates the general second-order differenced signal
Make inline??
*/
{
  return 0.5*(
     hi*(1 + MFS_fastErf((j-ti)*vd))*expqd2i
    -hr*(1 + MFS_fastErf((j-tr)*vd))*expqd2r
  );
}
//------------------------------------------------------------------------------
double zerothOrderSignal(double hi, double hr, int j, double ti, double tr,
                         double vd, double expqd2i, double expqd2r)
/*
171021.
Little function that generates the general non-differenced signal
Make inline??
*/
{
  return 0.5*(
    hi*(
        MFS_fastErf((j-ti-2)*vd)
      - 2*MFS_fastErf((j-ti-1)*vd)
      + MFS_fastErf((j-ti)*vd)
    )*expqd2i
    -
    hr*(
        MFS_fastErf((j-tr-2)*vd)
      - 2*MFS_fastErf((j-tr-1)*vd)
      + MFS_fastErf((j-tr)*vd)
    )*expqd2r
  );
}

//------------------------------------------------------------------------------
double expqd2(std::string model, double rx, double ry, double rz, double xi,
              double eDM[3], double d, double a, double R)
/*
171023.
Little function that calculates:
  Exp[-(q/d)^2],
where q is the impact parameter.
Note: impact parameter is 0 for walls.
And it is different for strings/monopoles.
*/
{
  if(model=="dwall"||model=="twall") return 1.;

  double cosgamma=cosGamma(rx,ry,rz,xi,eDM,a);
  double rp=perpendicularClock(rx,ry,rz,xi);
  double qd2=0;
  if(model=="monopole")
    qd2 = (pow(rp,2) + pow(R,2) - 2.*rp*R*cosgamma)/(pow(d,2));
  else if(model=="string")
    qd2 = pow((R - rp*cosgamma)/d,2);

  return MFS_fastExp(-qd2);
}

//******************************************************************************
double cosGamma(double rx, double ry, double rz, double xi, double eDM[3],
                double a)
/*
170816.
Calculates the cosine of the angle gamma, which is needed to calculate the
impact parameter.

INPUT:
  rx  ::  x (and y and z) position of the satellite (or reference clock)
  xi  ::  projection of sat (or ref) onto eDM
  eDM ::  The incident direction of DM object (unit vector)
  a   ::  Incindent "angle" of DM object (impact parameter angle)
  R   ::  "Earth" impact parameter of DM object
OUTPUT (on return):
  Cosine of the angle gamma.
*/
{
  double tempx = rz - eDM[2]*xi;
  double tempy;
  if(fabs(tempx)>0.05){  // double check number??
    tempy = eDM[0]*ry - eDM[1]*rx;
  }else{
    tempx = rx - eDM[0]*xi;
    tempy = eDM[1]*rz - eDM[2]*ry;
  }
  double beta=atan2(tempy,tempx);
  return cos(a-beta);
}

//******************************************************************************
double projectClock(double rx, double ry, double rz, double eDM[3])
/*
170816.
Tiny function that projects clock positions onto eDM (n)
*/
{
    //std::cout << "eDM[0]: " << eDM[0] << "\n";
    //std::cout << "eDM[1]: " << eDM[1] << "\n";
    //std::cout << "eDM[2]: " << eDM[2] << "\n";
  return rx*eDM[0]+ry*eDM[1]+rz*eDM[2];
}



//******************************************************************************
double perpendicularClock(double rx, double ry, double rz, double xi)
/*
170816.
Calculates the perpendicular distance between sat and n
*/
{
  return sqrt(rx*rx+ry*ry+rz*rz - xi);
}



//******************************************************************************
int refOnlySignal(double t0, double v, double rref[3], double eDM[3],
                  int iClocks, int jmin, int JW, double s[MAXCLOCKS][MAXJW])
/*
170816.
Overloaded function.
Designed so that refOnlySignal may be called using t0 OR tr!
.....
INPUTS
   t0       ::  Time event crosses ECI0
   v        ::  Vel of DM object (just to work out tr from t0)
   rref     ::  Ref clock position
   iClocks  ::  Number of clocks
   jmin     ::  Offset: start position of signal window in epochs
   JW       ::  Length of signal window
INPUTS/OUTPUTS
   s        ::  2D Array that holds the signal for each clock

=======Change log=====

*/
{


  //Project reference clock onto DM direction (to find tr)
  double xref=0;
  for(int k=0;k<3;k++)xref+=eDM[k]*rref[k];

  //time "signal" hits the reference clock
  double tr=t0-(xref/v);

  //call overloaded function using tr:
  refOnlySignal(tr,iClocks,jmin,JW,s);

  return 0;
}
//------------------------------------------------------------------------------
int refOnlySignal(double tr, int iClocks, int jmin, int JW,
                  double s[MAXCLOCKS][MAXJW])
/*
170816.
170717. Code from Kelly.
Injects a "reference only" signal.
For now, first-order differenced only.
This will inject a jump at the epoch jr.
NOTE: this function takes jmin (short array) into account.
It should be passed the ACTUAL tr (NOT shofted tr)
Overloaded function.
Designed so that refOnlySignal may be called using t0 OR tr!
.....
INPUTS
   tr       ::  Time event crosses reference clock
   iClocks  ::  Number of clocks
   jmin     ::  Offset: start position of signal window in epochs
   JW       ::  Length of signal window
INPUTS/OUTPUTS
   s        ::  2D Array that holds the signal for each clock

=======Change log=====

*/
{

  //making a jump at the reference clock
  for(int i=0;i<iClocks;i++){
    //Note: subtract jmin to shift offset because of short array
    int jr=(int)ceil(tr)-jmin;
    for(int j=0;j<JW;j++){
      //NB: have to loop through all epochs, since s may not be "zerod"
    	if(j==jr){
        s[i][j] = 1;
      }else{
        s[i][j] = 0;
      }
    }
  }//END loop over clocks

  return 0;
}




//******************************************************************************
int randJumpSignal(double t0, int iClocks, int jmin, int JW,
                   double s[MAXCLOCKS][MAXJW])
/*
170717. Code from Kelly.
Injects one jump per clock at a random epoch centered around t0.
For now, first-order differenced only.
.....
INPUTS:
   t0      ::  Time event crosses ECI0
   iClocks::  Number of clocks
   jmin    ::  Offset: start position of signal window in epochs
   JW      ::  Length of signal window
INPUTS/OUTPUTS:
   s       ::  2D Array that holds the signal for each clock

=======Change log=====

*/
{

  double ti=0; //time signal hits ith clock (random, +/- 4 epochs from t0)
  int ji=0; //the following epoch (discritised sample point)
  //Note: subtract jmin to shift offset because of short array

  //making jumps for random clocks
  for(int i=0;i<iClocks;i++){
    ti=t0-jmin+MFS_randDouble(-4,4);
    ji=(int)ceil(ti); //round up to epoch
    for(int j=0;j<JW;j++){
      //NB: have to loop through all epochs, since s may not be "zerod"
      if(j==ji){
        s[i][j]=1;
      }else{
         s[i][j]=0;
      }
    }
  }

  return 0;
}
