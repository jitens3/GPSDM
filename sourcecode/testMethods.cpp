/*
PROGRAM: testMethods

Still under developement.

This program is a combination of two older programs testLikelihoods and testSNR.
This program is designed to test the "likelihoods" odds ratio program as well
as the matched-filter technique "snr" program.
   * It generates random time-series for a user-specified number of clocks.
The noise can be either white (with given standard-deviation), or it can use the
known PSDs from GPS clocks to simulate realistic time series.
It can do this using either the overall averaged PSD, or by randomly assigning
each sattelite an SVN, and using the correct PSD for that SVN.
   * It places each sattelite in a realistic (though static) position in space
   * Either inject ideal signals into the data (based on user input) to test
     how well program can picj up events, or
   * not inject events, to test for false-positives
   * Then calls 'methods' program, which calculates the odds ratio or SNR.

Note: you may limit the number of cores you use for parallelisation.
If you are running tests on one of the servers this is a good idea, so you
don't just use up all the servers resources. It is also a good idea if you are
running the code on your own PC to keep at least 1 core free, so you can still
use the PC for other tasks without it being too slow.

The program optionally makes use of power spectrums (for simulating clocks) and
the inverse covariance functions (for covariance matrix). These have been
calculated already, and can be downloaded from github (see main README for
more detailed git instructions)
   * https://github.com/benroberts999/PowerSpectrums

e.g.:
   * $ wget http://github.com/benroberts999/PowerSpectrums/archive/master.tar.gz
   * $ tar xf master.tar.gz -C ./ --strip-components=1

NOTE: At the moment, vmin and vmax are hard-coded in.
J_W (the data "window") is worked out using vmin.
Jw (number of points in Chi^2 sum) is important.
Open question...


=======Change log=====
...
170719- Adding "bad signal" options.
170831- Starting updated to work with class. Major changes.
      - Killed 'calcLikelihoods' (non-MC method)
      - Works. But, lots still to do!
171004- Simulates cross-clock-correlations by swapping to white ref clock
171018- Cross-clock-correlation can now be included in the likelihood
171024- Mixed difference, assume relative couplings, input monopole/strings etc
171107- Parameter estimation
190521- Combined testLikelihoods and testSNR to this with method option in .dat
190619- Added exact inversion of cavar matrix option
****** To Do ******
XXX subtract weighted mean of clocks from same type?
XXX Best way to read in data from either side? -- JplGpsData

*/
#include "methods.h"

//***********************************************************************
//***********************************************************************
void printParameters(std::vector<int> numClks, std::string sref, int icov,
    bool use_ccc, bool sub_w_mean, int in_receivers, double sta_sd,
    bool simulate_ccc, double ref_sig,
    int irandSVN, double sig, int idiff, double jh,
    bool bad_ref, int iRand_bad, double bad_amp,
    int NMC, double h0,
    double oddscut, double SNRcut, std::string label, int method);

//***********************************************************************
//***********************************************************************

int main (void){

	//for timing:
	time_t start,end; //total program time
	time (&start);
	long int margTime=0;  //time just the integrations

	//Input Parameters:
	std::string path_to_psd;//location of input PSD and ACF files
	std::string path_to_cdf;//location of input numeric prior CDF files
	int iRbIIF,iRbIIA,iRbIIR,iRbII,iCsIIF,iCsIIA,iCsII; //number of each clock
	std::string sref; //which ref clock to simulate
	std::string acf_label; //which PSD/ACF files (which label) to use
	int iw;     //number of white clocks to simulate
	double sig;//s.d. for white noise
	//Use ACF for specific SVN? Use sd from file?
	int irandSVN; // Use randomised SVNs?
	int icov, read_acf, iccc; //Covariance? Read ACF from file? Use cross-clock
	int isubWmean;
	int in_exact;
	int biasprint;
	int in_receivers; //number of receiver clocks to simulate (white noise)
	double sta_sd;    //s.d. of receiever clocks
	int isimccc; //simulate cross-clock-correlations?
	double ref_sig; //(optional) s.d. for ref clock (when simulate_ccc is true)
	int i_iterations;  //numer of itterations of entire thing (for tests!)
	double jh; //input "jump" height (magnitude)
	//input event:
	double int0,inV,inT,inP,inLd; //t0, vel. (km/s), theta/phi (/pi), log(d)
	double indt0,indV,indTP; //ERROR in vel., theta/phi(/pi)
	std::string model;  //model (wall/string/monopole)
	double inRdm, inadm; //impact parameter terms
	int idiff; //difference (1 or 2)
	std::string coupling; //assume specific coupling?
	int iprior; //use full priors?
	double h0; //
	int jbeg,jend; //which epochs to include.
	std::string label; //label for output files
	int ilimit_cores; //limit the number of cores for openMP XXX THREADS!
	int NMC; //Number grid points for mixed Monte-Carlo.
	int iwhichpara; //paralise over itterations, or inside likelihoods
	int ibad_ref=0;   //Insert a "bad" ref jump? (yes/no)
	int iRand_bad=0;  //Insert "bad" random sat. jumps. (number)
	double bad_amp=0; //Magnitude of "bad" signals (sign randomised)
	int method = 0; //method to use: 0-odds, 1-snr
	double oddscut; //Odds ratio cut-off for "found"
	double SNRcut; //snr cut-off value for "found"
	int i_param;  //parameter estimation

	//Read in the input .dat file
	std::ifstream fInput;
	std::string junk;
	fInput.open ("testMethods.dat");
		fInput >> path_to_psd;                            getline(fInput,junk);
		fInput >> path_to_cdf;                            getline(fInput,junk);
		fInput >> iRbIIF>>iRbIIA>>iRbIIR>>iRbII;          getline(fInput,junk);
		fInput >> iCsIIF>>iCsIIA>>iCsII;                  getline(fInput,junk);
		fInput >> sref >> irandSVN >> acf_label;          getline(fInput,junk);
		fInput >> iw >> sig;                              getline(fInput,junk);
		fInput >> in_receivers >> sta_sd;                 getline(fInput,junk);
		fInput >> isimccc >> ref_sig;                     getline(fInput,junk);
		fInput >> idiff >> icov >> read_acf >> iccc
		       >> isubWmean >> in_exact;               getline(fInput,junk);
		fInput >> i_iterations;                           getline(fInput,junk);
		fInput >> jbeg >> jend;                           getline(fInput,junk);
		fInput >> jh >> inLd;                             getline(fInput,junk);
		fInput >> int0 >> inV >> inT >> inP;              getline(fInput,junk);
		fInput >> indt0 >> indV >> indTP;                 getline(fInput,junk);
		fInput >> model >> inRdm >> inadm;                getline(fInput,junk);
		fInput >> ibad_ref >> iRand_bad >> bad_amp;       getline(fInput,junk);
		fInput >> coupling >> iprior;                     getline(fInput,junk);
		fInput >> NMC >> h0;                              getline(fInput,junk);
		fInput >> method >> oddscut >> SNRcut 
			   >> i_param >> biasprint;					  getline(fInput,junk);
		fInput >> ilimit_cores >> iwhichpara;             getline(fInput,junk);
		fInput >> label;                                  getline(fInput,junk);
	fInput.close();
//-------------------------------------------------------------------
	std::cout<<"\n\n          ###########################################"
             << "\n              Testing search simulation program        "
             <<" \n          ###########################################\n\n";

	//helps remove errors if wrong inputs given
	if(sig<0)sig=fabs(sig);
	if(sig==0)sig=0.05;
	if(i_iterations>MAXITS)i_iterations=MAXITS;
	if(idiff<0||idiff>3) idiff=1;

	//XXX Temporary fix
  	if(idiff!=1&&isubWmean==1){
    	std::cout<<"\nSorry, for now, I cannot use 2nd order differencing and sub "
            	 <<" mean :( \nProblem is in generateSignal.\n";
    	return 1;
  	}

  	//XXX Temporary fix!
  	if(i_param==1 && iwhichpara==1 && jh==0){
    	std::cout<<"Bad input: changeing iwhichpara to 0\n";
    	iwhichpara=0;
  	}

  	if(jend==0) jend=2880;

  	//if invalid inputs, set inversion to approx and method to odds-ratio
  	if (in_exact != 1) in_exact = 0;
  	if(method != 1) method = 0;

  	//Make sure input for models is ok. Will default back to thin-walls!
  	if(model=="wall") model="dwall";
  	if(model=="mono") model="monopole";
  	if(model!="dwall"&&model!="string"&&model!="monopole") model="twall";

	//Makes sure input coupling is ok.
	//If not alpha, m_e, or m_p, assumes all clocks the same.
	if(coupling=="a") coupling="alpha";
	if(coupling=="ep"||coupling=="mep"||coupling=="e") coupling="me";
	if(coupling=="q") coupling="mq";
	if(coupling!="alpha"&&coupling!="me"&&coupling!="mq") coupling="";

	//Use real, or "flat" priors for velocity/angles?
	bool use_flat_priors=false;
	if(iprior==0) use_flat_priors=true;
	
	//Simulate cross-clock-correlations, by generateing a (white-noise) reference
	//clock, [low sd=0.0075 ns], and swapping to it.
	bool simulate_ccc=false;
	if(isimccc==1) simulate_ccc=true;
	
	//Include the cross-clock-correlation into covariance.
	bool use_ccc = false;
	if(iccc==1) use_ccc=true;
	bool sub_w_mean = false;
	if(isubWmean==1) sub_w_mean=true;
	
	//"Bad" signals
	bool bad_ref=false;
	if(ibad_ref>0)bad_ref=true; //insert a "bad" ref jump
	
	//At least for now, cannot both insert bad reference and subtract mean!
	if(bad_ref && sub_w_mean){
	  std::cout<<"\n\n FAIL 222 in testLikelihoods:\nBecause reasons, cannot "
	           <<"insert a Bad_Ref signal, and subtract mean!\n\n";
	  return 1;
	}

	//Parallelise over iterations?
	bool b_it_para=false;
	if(iwhichpara==1)b_it_para=true;

  	//Lets the user limit the number of threads in use.
  	//Use get_.._threads instead?? (threads instead of cores)
  	int i_threads=1;
  	#if defined(_OPENMP)
    	int i_cores=omp_get_num_procs();
    	if(ilimit_cores>0)omp_set_num_threads(ilimit_cores);
    	#pragma omp parallel
    	{
      		i_threads=omp_get_num_threads();
    	}
    	printf("Parallelised (OMP) using %i/%i cores\n\n",i_threads,i_cores);
 	#endif

    if(method == 1){
    	std::cout << "Using SNR-MFT Method ";
    	if (in_exact == 1) {std::cout << "with exact E inversion\n\n";}
    	else {std::cout << "with appriximate E inversion\n\n";}
    }
    else{
    	std::cout << "Using Odds-Ratio Method ";
    	if (in_exact == 1) {std::cout << "with exact E inversion\n\n";}
    	else {std::cout << "with appriximate E inversion\n\n";}
    }

	// //some testing/formatting options: print/plot the data?
	// bool printPlots=false; //plot data for first iteration
	// if(b_it_para)printPlots=false;

	//Label attached to the input PSD/ACF files.
  	//If input is "na", then there is no label.
  	if(acf_label=="na")acf_label="";

  	//simulate random SVN, or use the average?
  	bool useRandSVN=false;
  	if(irandSVN==1)useRandSVN=true;

  	//Store 'avergae' value for s.d. of the satellite clocks (not stations!)
  	double avsig; // used for outputting results only

  	//store number of each sat/clock to simulate
  	//note: in_receivers separate!
  	std::vector<int> numClks(8);
    	numClks[0]=iRbIIF;
    	numClks[1]=iRbIIA;
    	numClks[2]=iRbIIR;
    	numClks[3]=iRbII;
    	numClks[4]=iCsIIF;
    	numClks[5]=iCsIIA;
    	numClks[6]=iCsII;
    	numClks[7]=iw;    //white clocks: must be last entry!

  	//count the clocks (just for progress bar + tests):
  	int temp_clocks=in_receivers;
  	for(size_t iType=0; iType<numClks.size(); iType++)
    	temp_clocks+=numClks[iType];

  	//number of epochs (for GPS simulator + odds array etc.)
  	int in_epochs=2880; //Default. Usually less:
  	if(jend+150<2880) in_epochs = jend+150; //JW is ~85
  	if(iccc /*&& !in_exact*/) in_epochs=2880; //need large J to work out ccc accurately!

  	//if _only_ have white clocks, can simplify some things.
  	bool white=false; //white noise? Or from PSD
  	if(iw+in_receivers==temp_clocks) white=true;
  	if(white) icov=0; //no need to run covariate


  	//Output input params etc. to the screen
  	printParameters(numClks,sref,icov,use_ccc,sub_w_mean,in_receivers,sta_sd,
                  simulate_ccc,
                  ref_sig,irandSVN,sig,idiff,jh,bad_ref,iRand_bad,bad_amp,NMC,
                  h0,oddscut, SNRcut, label, method);
  	if(temp_clocks==0) return 1;

  	//array to hold all the SNR values (for all itterations)
  	static double allSNR[MAXITS][2880] = {0};
  	//array to hold all the odds rations (for all itterations)
  	static double allodds[MAXITS][2][2880] = {0};
  	//stores where and how big the injected event was
  	static double injevent[MAXITS][5];

  	//open the output file for the parameter esimation:
  	bool estimate_params=false;
  	if(i_param==1) estimate_params=true;
  	std::vector<std::vector<double> > all_best_params;


/*===========================================================================
=========================== LOOPING ITERATIONS ==============================
===========================================================================*/
  	//loop through the 'iterations'
  	#pragma omp parallel for //if(b_it_para)
  	for(int itn=0; itn<i_iterations; itn++){

    	MSC_progressBar("Iterations",30,itn,i_iterations,b_it_para);


		// *************************************************************
		// ******* Below this is for generateing simulated data! *******
		// *************************************************************

    	//create blank GPS data object
   		JplGpsData data;

    	//Run GPS simulator:
    	int jplok=data.gpsSimulator(path_to_psd,acf_label,numClks,sig,
                                in_receivers, sta_sd,
                                sref,in_epochs,
                                useRandSVN,simulate_ccc,ref_sig,!b_it_para);
    	if(jplok==2){
      		std::cout<<"FAIL 414 in testLikelihoods: Something wrong! Check PSD?\n\n";
      		i_iterations=0; //try to exit the OMP'd for loop
      		continue;
    	}


    	//Work out the relative coupling sensitivities
    	data.formRelativeCouplings(coupling);


    	//Generate fake event parameters:
    	//Randomise sign of input event
    	int pm=1;
	    //{ /* only uses positive for now */
	    //  double temp=MFS_randDouble(-1,1);
	    //  if(temp<0)pm=-1;
	    //}

    	//event size, in ns:
    	double amp = pm*jh;
    	//time when DM crosses Ref:
    	double t0 = MFS_randDouble(int0-indt0,int0+indt0);
    	//DM speed: km/epoch
    	double vDM = fabs(MFS_randDouble(inV*30.-indV*30.,inV*30.+indV*30.));
    	//theta angle for DM:
    	double thetaDM = MFS_randDouble(inT*PI-indTP*PI,inT*PI+indTP*PI);
    	if(thetaDM<0) thetaDM = -thetaDM;
    	while(thetaDM>PI) thetaDM = 2*PI - thetaDM;
    	// phi angle for DM:
    	double phiDM = MFS_randDouble(inP*PI-indTP*PI,inP*PI+indTP*PI);

      // //time when DM crosses Ref:
      // t0 = 830.5;
      // //DM speed: km/epoch
      // vDM = 300*30;
      // //theta angle for DM:
      // thetaDM = 0.236;
      // if(thetaDM<0) thetaDM = -thetaDM;
      // while(thetaDM>PI) thetaDM = 2*PI - thetaDM;
      // // phi angle for DM:
      // phiDM = 1.741;

    	//DM object width, d:
    	double d=0, Rdm=0, adm=0;
    	if(model!="twall"){
    		//double log_d = MFS_randDouble(inLd-indLd, inLd-indLd);
    		d = pow(10,inLd);
    		//impact parameter: only for monopoles/strings:
    		// no "range" for these
    		Rdm = inRdm*RGPS;
    		adm = inRdm*PI;
    	}

    	//store the generated signal params (for testing/output etc.)
    	//add d, R, a to these!? Nah.. it's OK (done in param est)
    	injevent[itn][0]=t0;
    	injevent[itn][1]=amp;
    	injevent[itn][2]=vDM/30.;     //in km/s
    	injevent[itn][3]=thetaDM/PI;  //in units of 'pi'
    	injevent[itn][4]=phiDM/PI;
    	if(jh!=0&&!b_it_para){
    	  	printf("DM event: h=%.4f ns, t0=%.2f, (v,theta,phi)=%.1fkm/s, "
    	         "%.2fpi,%.2fpi\n",amp,t0,vDM/30.,thetaDM/PI,phiDM/PI);
    	}

    	//========Print input params to file=========
    	if(i_param == 1){
    	  std::ofstream injEvent;
    	  std::string inj = "./results/injE.txt";
    	  injEvent.open(inj, std::ios_base::app);
    	     injEvent << itn <<"\t"<<injevent[itn][0] << "\t"<< injevent[itn][1] << "\t"<< injevent[itn][2] << "\t"<< injevent[itn][3] << "\t"<< injevent[itn][4] << "\n";   
    	  injEvent.close();
    	}
    	// Unit vector, in direction of DM velocity
    	double eDM[3];
    	eDM[0]=sin(thetaDM)*cos(phiDM);  //x
    	eDM[1]=sin(thetaDM)*sin(phiDM);  //y
    	eDM[2]=cos(thetaDM);             //z

    	//Inject Actual event into the clock data:
    	int jmin=int(t0-30);
    	int jmax=int(t0+30);
    	int JW_gen=jmax-jmin+1; //new!
    	int iDGen=0;
    	//Store the signal in 2D array:
    	double s[MAXCLOCKS][MAXJW]; //zeroing not needed!

    	//Form  rsat and rref arrays:
    	//Here, I just take the "slice" of positions that I need to calculate
    	// the signals.
    	// Note: this is the same as is done in the likelihoods program
    	double rsat[MAXCLOCKS][3];
    	double rref[3];
    	int j0=(int)ceil(t0); //index to take
    	int iClocks=data.num_clocks;
    	for(int i=0; i<iClocks; i++){
    	  	for(int ix=0; ix<3; ix++)
    	    	rsat[i][ix]=data.pos[i][j0][ix];
    	}
    	for(int ix=0; ix<3; ix++)
    	  	rref[ix]=data.refpos[j0][ix];

    	//generate the relatvant signal:
    	generateSignal(model,s,rsat,rref,iClocks,jmin,JW_gen,iDGen,data.keff,t0,
                   vDM,eDM,d,Rdm,adm);

   

      //read in fixed signal
      //std::ifstream sigRead("./signal-d0.txt");
      //for(int i=0;i<iClocks;i++){
      //  for(int j=jmin;j<=jmax;j++){
      //      sigRead >> s[i][j-jmin];
      //  }
      //}
      //sigRead.close();

      int Ns = 29;

      //for some reason, needs to be <= for Ns>1
      for(int i=0;i<Ns;i++){
        for(int j=jmin;j<=jmax;j++){
          data.bias[i][j]+=amp*s[i][j-jmin];
        }
      }
      

       std::ofstream templ;
       std::string fileX = "./results2/s.txt";
       templ.open(fileX , std::ios_base::app);
       for(int a=0;a<data.num_clocks;a++){
         for(int b=1;b<MAXJW;b++){
           //templ<<s[a][b]-s[a][b-1]<<" ";
           templ<<s[a][b]<<" ";
         }
         templ<<"\n";
       }
       templ<<std::endl;
       templ.close();


    	// *************************************************************
    	// ******* Above this is for generateing simulated data! *******
    	// *************************************************************

    	// Either subtract a weighted mean and difference, or just difference.
    	// NB: differeing happens inside subtractWeightedMean
    	if(sub_w_mean) data.subtractWeightedMean(idiff);
    	else data.differenceData(idiff);


    	
    	/////////////////////////////////////////////////////////////////
    	//Inject "Bad" signals (must be done after differenceing for now)
    	//"sign" of amplitude randomised:
    	double temp_bad_sign=MFS_randDouble(-1,1);
    	if(temp_bad_sign<0)bad_amp*=-1;
    	if(bad_ref){
    	  	// XXX NOTE: This should NOT be done if 'subtracting mean'
    	  	//inject a "bad" ref-only jump
    	  	//note: sign is opposite for ref jump!
    	  	refOnlySignal(t0,vDM,rref,eDM,iClocks,jmin,JW_gen,s);
    	  	//Once generated the 'bad' signal, inject it:
    	  	for(int i=0;i<iClocks;i++){
    	    	for(int j=jmin;j<=jmax;j++){
    	      		data.bias[i][j]+=-1*bad_amp*s[i][j-jmin];
    	    	}
    	  	}
    	}
    	int nbi=0;
    	while(nbi<iRand_bad){
    	  	//loop, so we can inject a number of "bad" random clock jumps.
    	  	//Each iteration will inject exactly 1 jump per clock, at a random epoch
    	  	randJumpSignal(t0,iClocks,jmin,JW_gen,s);
    	  	//Once generated the 'bad' signal, inject it:
    	  	for(int i=0;i<iClocks;i++){
    	    	for(int j=jmin;j<=jmax;j++){
    	      		data.bias[i][j]+=bad_amp*s[i][j-jmin];
    	    	}
    	  	}
    	  	nbi++;
    	}
    	// End: inject "bad" signals
    	/////////////////////////////////////////////////////////////////

    	//calculate s.d. of data: (must be done /after/ all signals injected):
    	double t_avsig = data.calculateStdDev();

      // std::ofstream sdevout;
      // std::string fileX2 = "./results/sdev.txt";
      // sdevout.open(fileX2 , std::ios_base::app);
      // for(int a=0;a<data.num_clocks;a++){
      //     sdevout<<data.sdev[a]<<" ";
      //   }
      // sdevout<<std::endl;
      // sdevout.close();
    	//nb: This does not include the simulate base station clocks! ???
    	// Doesn't distinguish between s1 or s2
    	// if(!b_it_para){
    	//   std::cout<<"Average sigma="<<t_avsig<<"ns\n";
    	// }
    	//Average s.d. over all itterations:
    	avsig += t_avsig/i_iterations;

    	//Read in the ACF files, form ACF array:
    	if(icov==1){
    	  	if(read_acf==1){
    	    	data.formAutoCorrelation(path_to_psd, acf_label); //from ACF from file
    	  	}
    	  	else{
    	    	data.formAutoCorrelation(); //calulate "on-the-fly"
    	  	}
    	}


    	//Calculate the cross-clock correlations:
    	if(use_ccc){
        data.calculateCrossClockCorrelation();
        //std::cout<<"b0 = "<<data.b0[0]<<" sigma_x^2 = "<<ref_sig*ref_sig<<std::endl;
      }

      if(biasprint){
          std::ofstream biasData("./results/biasData.txt");
        for(int i=0;i<iClocks;i++){
          for(int j=jbeg;j<=jend;j++){
              biasData<<data.bias[i][j+1];
              if(j!=jend){biasData<<" ";}
          }
          biasData<<std::endl;
        }
        biasData.close();
      }

    	//Read in the numeric CDFs (for the priors)
    	std::string path_v_cdf,path_psi_cdf;
    	if(model=="twall"||model=="dwall"||model=="string"){
    	  	//Domain walls and strings
    	  	path_v_cdf  =path_to_cdf+"/vperp.cdf";
    	  	path_psi_cdf=path_to_cdf+"/psiperp.cdf";
    	}
    	else{
    		//Monopoles
    		path_v_cdf  =path_to_cdf+"/v.cdf";
    		path_psi_cdf=path_to_cdf+"/psi.cdf";
    	}
    	NumericCdfInverse v_prior(path_v_cdf);
    	NumericCdfInverse psi_prior(path_psi_cdf);
    	//option to use 'flat' priors for v,angles:
    	if(use_flat_priors){
    		//over-write the above priors with "flat" priors.
    		v_prior = NumericCdfInverse("flat",25,750);
    		psi_prior = NumericCdfInverse("SolidAngle");
    	}

    	//define and initialise SNR (for this itteration)
    	double snr[data.num_epochs]= {0};

    	//define and initialise odds (for this itteration)
    	std::vector< std::vector<double> >
    	  odds(2, std::vector<double>(data.num_epochs));

    	// Work out the window J_W used for the analysis.
    	// vmin_window should be an input!
    	//double vmin_window = 25.*30; //slowest velocity we can 'see' in theory
    	//int jw = 2*((int)ceil(1.1*RGPS/vmin_window))+1; //'window J_W'
    	//if(jw>MAXJW) jw = MAXJW;  //safety-check! should giv3e warning?
      int jw = MAXJW;

      //std::cout<<"jw = "<<jw<<std::endl;

    //Calculate exact inverse fof covariance matrix if exact method chosen
    int size = data.num_clocks*jw;
    double *Einv = (double *)malloc(size * size * sizeof(double));
    if(in_exact){
      double max_cov = 5;
      double acf_cut = 0.025;
      calcEinv(Einv, data, jw, icov, max_cov, acf_cut, 0,0, ref_sig, 0, 0);
      //calcEinvWhiteNoise( Einv, data, jw, sig, ref_sig);
    }
    else{
      delete [] Einv;
    }

		/* =================================================================
		==================	IF ODDS RATIO METHOD CHOSEN ====================
		==================================================================*/

	    //*********************Calc statistics*************************
	    if(method == 0){
	      //Calculate the likelihoods and odds ratio
	    	auto begin = std::chrono::high_resolution_clock::now();//timing
	    	MClikelihoods(model,data,icov,use_ccc,sub_w_mean,odds,jbeg,jend,jw,NMC,h0,
	                    v_prior,psi_prior, Einv, in_exact);
	      	auto end = std::chrono::high_resolution_clock::now();
	      	margTime+=std::chrono::duration_cast
	                <std::chrono::nanoseconds>(end-begin).count();
	
	      	// Transfer the arrays etc.
	      	//for(int j=0; j<data.num_epochs; j++){
	      	for(int j=jbeg; j<=jend; j++){
	        	//Do parameter estimation
	        	if(odds[0][j]>oddscut && estimate_params){
	          		if(jh!=0){
	            	//doing true-positive test! (only do 1 epoch!)
	            	int j0 = (int)ceil(t0);
	            	double temp_best_odds=0;
	            	int temp_best_j=-1;
	            		for(int jj=j0-1; jj<=j0+1; jj++){
	              		
	              			if(odds[0][jj]>temp_best_odds){
	                			temp_best_odds=odds[0][jj];
	                			temp_best_j=jj;
	              			}
	            		}
	            	if(j!=temp_best_j) continue;
	          		}
	          		//Perform the parameter estimation (with 4x grid)
	          		int NMC2=4*NMC;
	          		std::vector<double> best_params;
	          		likeParameterEstimation(model,data,icov,use_ccc,sub_w_mean,itn,j,jw,NMC2,h0,
	                        v_prior,psi_prior,best_params, Einv, in_exact);
	          		best_params.push_back(odds[0][j]);  //add odds ratio (integrated)
	          		best_params.push_back(odds[1][j]);  //add log likelihood
	          		all_best_params.push_back(best_params);
	        	}
	        	//Transfer odds array for each itteration.
	        	for(int l=0; l<2; l++){
	          		allodds[itn][l][j]=odds[l][j];
	       		}
	    	}
		}
		/* =================================================================
		==================	IF SNR MFT METHOD CHOSEN ====================
		==================================================================*/
	    if(method == 1) {
	      	//calculate snr values
	      	auto begin = std::chrono::high_resolution_clock::now();//timing
	      	MCsnr(model, data, icov, use_ccc,
	        		sub_w_mean, snr, jbeg, jend, jw,
	        		NMC, h0, v_prior, psi_prior, in_exact, Einv);
	      	auto end = std::chrono::high_resolution_clock::now();
	      	margTime += std::chrono::duration_cast
	        	<std::chrono::nanoseconds>(end - begin).count();
	
	
	      	//Set allSNR matrix values, for each iteration and epoch
	      	for(int iset=jbeg; iset<jend; iset++) allSNR[itn][iset]=snr[iset];
	    
	      	//open the output file for the parameter esimation:
	      	if (i_param == 1){
	      		//Perform the parameter estimation if event is found
	        	int NMC2=4*NMC;
	            
	            //find largest snr near j0 (over SNRcut)
	            double tempBestSnr = snr[j0-1];
	        	int tempBestEp = j0-1;
	        	if(fabs(snr[j0-1])>SNRcut||fabs(snr[j0])>SNRcut||fabs(snr[j0+1])>SNRcut){
	        		if(fabs(snr[j0])>fabs(tempBestSnr)){
	        			tempBestSnr = snr[j0];
	        			tempBestEp = j0;
	        		}
	        		if(fabs(snr[j0+1])>fabs(tempBestSnr)){
	        			tempBestSnr = snr[j0+1];
	        			tempBestEp = j0+1;
	        		}
	        	}
	        	//perform parameter estimation (increased iterations should better estimates)     
	        	snrParameterEstimation(model,data,icov,use_ccc,sub_w_mean, itn, tempBestEp,jw,NMC2,h0,
	                          v_prior,psi_prior, in_exact, Einv, 0);
	      	}
	    }
	    //*************************************************************
	
	}//END loop through iterations.
	MSC_progressBar("Iterations",30);

  //for(int i = jbeg; i <=jend; i++){
  //  std::cout << allSNR[0][i] << std::endl;
  //}
	////////////////////////////////////////////////////////////////////////////////

  	//re-count the clocks (some may not have been found) - just for safety
  	if(!white){
    	temp_clocks=in_receivers;
    	for(int iType=0;iType<7;iType++)temp_clocks+=numClks[iType];
    	if(temp_clocks==0){
      		std::cout<<"Couldn't generate any clocks! Check PSD?"<<std::endl;
      		return 1;
    	}
  	}

/*===================================================================================
================================= RESULTS ===========================================
===================================================================================*/

	/* =================================================================
	==================	IF ODDS RATIO METHOD CHOSEN ====================
	==================================================================*/
	if(method==0){
  	//Following prints results +related info to the screen after program is done

  		//Plot results (odds):
  		std::string plotname="odds-";
  		if(white)plotname=plotname+MSC_padIntString(temp_clocks)+"WN"
                    +MSC_padIntString(int(sig*100))+"-";
  		if(!white)plotname=plotname+MSC_padIntString(temp_clocks)+sref+"-";
  		plotname=plotname+MSC_padIntString(int(jh*100))+"-"+label;

		//------------------------------
  		//Output results to screen (false positives)
  		if(jh==0&&iRand_bad==0&&!bad_ref){
    		std::cout<<"\n How many false positives?\n";
    		for(int ii=1; ii<=20; ii++){
    			int icut;
    			if(ii<=10) icut=ii;
    			else if(ii<=15)icut=(ii-10)*20;
    			else if(ii==16)icut=500;
    			else if(ii==17)icut=1000;
    			else if(ii==18)icut=10000;
    			else if(ii==19)icut=100000;
    			else if(ii==20)icut=1000000;
    			double cut=double(icut);
    			int numOver=0;
    			for(int it=0; it<i_iterations; it++){
        			for(int j=jbeg; j<=jend; j++){
          				if(allodds[it][0][j]>=cut) numOver++;
        			}
      			}
      			double avgno=double(numOver)/double(i_iterations);//f.p. in window
      			double nopd=2880*(avgno/double(jend-jbeg+1)); //f.p. per day
      			double nopy=nopd*365.25; //f.p. per year
      			printf("cut=%7i: false pos: %8.3f /day = %8.1f /year\n",
             		icut,nopd,nopy);
    		}
  		}
  		std::cout<<"\n";

		//------------------------------
  		//Output results to screen (true positives, & "bad" jumps)
  		int numFound=0;
  		if(jh!=0){
    		double avgodds=0;
    		double avglogodds=0;
    		//int numFound=0;
    		std::cout<<std::endl;
    		for(int k=0; k<i_iterations; k++){
    			double bestodds=0;
      			int besti=0;
      			
      			for(int i=jbeg; i<=jend; i++){
        			if(allodds[k][0][i]>bestodds){
          				bestodds=allodds[k][0][i];
          				besti=i;
        			}
      			}
      			int j0=(int)ceil(injevent[k][0]);
      			printf("%i: %.2fns event, at t0=%.2f. ",k+1,injevent[k][1],
            			injevent[k][0]);
      			printf("v,th,ph=%.1f, %.1f,%.1fpi. ",injevent[k][2],injevent[k][3],
            			injevent[k][4]);
	      		printf("Odds=%.1g or %.1g or %.1g\n",allodds[k][0][j0-1],allodds[k][0][j0],allodds[k][0][j0+1]);
	      		if(allodds[k][0][j0]>oddscut||allodds[k][0][j0-1]>oddscut
	         			||allodds[k][0][j0+1]>oddscut)
	      		{
	        		std::cout<<"--> Found! at: "<<besti<<", odds="<<bestodds<<"   :)\n";
	        		numFound++;
	      		}
	      		else{
	        		std::cout.precision(2);
	        		std::cout<<"--> NOT found. "<<allodds[k][0][j0-1]<<", "
	                		 <<allodds[k][0][j0]
	                		 <<", "<<allodds[k][0][j0+1]<<".   *******  :("<<std::endl;
	      		}
	      		//Sum up average odds (only sum biggest of the 3) ?
	      		// This is probably reasonable..
	      		double tempodds=0,templogodds=0;
	      		for(int jj=j0-1;jj<=j0+1;jj++){
	        		if(allodds[k][0][jj]>tempodds){
	          			tempodds=allodds[k][0][jj];
	          			templogodds=log10(allodds[k][0][jj]);
	        		}
	      		}
	      		avgodds+=tempodds;
	      		avglogodds+=templogodds;
	    	}//End loop over iterations
   			avgodds/=i_iterations;
    		avglogodds/=i_iterations;
    		std::cout<<std::endl;
    		printf("Found %i/%i = %.1f%% of the ~%.3fsig events.\n",numFound,
           		i_iterations,double(numFound*100)/double(i_iterations),jh/avsig);
    		if(iRand_bad!=0||bad_ref){
      			std::cout<<"-> Injected BAD ~"<<fabs(bad_amp)/avsig<<" sigma signals: ";
      			if(bad_ref)std::cout<<"a ref jump";
      			if(bad_ref&&iRand_bad>0)std::cout<<", and ";
      			if(iRand_bad>0)std::cout<<iRand_bad<<" random jumps (per clock)";
      			std::cout<<std::endl;
    		}
    		printf("Average odds:     %.2g\n",avgodds);
    		printf("Average log odds: %.2g (=> %.2g)\n",avglogodds,pow(10,avglogodds));
    		std::cout<<std::endl;
  		}

    	//Just "bad" jumps (only when no signal injected!):
  		if(jh==0&&(iRand_bad!=0||bad_ref)){
		    //Print out jump parameters!
		    std::cout<<"-> Injected BAD ~"<<fabs(bad_amp)/avsig<<" sigma signals: ";
		    if(bad_ref)std::cout<<"a ref jump";
		    if(bad_ref&&iRand_bad>0)std::cout<<", and ";
		    if(iRand_bad>0)std::cout<<iRand_bad<<" random jumps (per clock)";
		    std::cout<<"\nFalse positive rate: \n";
		    const int nlogbins=10; // 1, 10, 100, 1000...
		    const int nlinbins=9; // 1,2,3,4,5,...
		    int badOdds[nlogbins]={0};
		    int badLinOdds[nlinbins]={0}; //2 .. 9
		    for(int k=0;k<i_iterations;k++){
		      	int j0=(int)ceil(injevent[k][0]);
		      	double bestO=0.1;
		      	//Go over entire "window", count spikes!
		      	for (int j=j0-10; j<=j0+10; j++){
		        	if(allodds[k][0][j]>bestO)
		          		bestO=allodds[k][0][j];
		      	}
		      	double lO = log10(bestO);
		      	if(lO <= 0) continue;
		      	int ilO = (int)lO;
		      	if(ilO >= nlogbins) ilO=nlogbins-1;
		      	for (int l=0; l<=ilO; l++) badOdds[l]++; //sum cumulatively!
		      	//"linear" odds (1-9)
		      	for (int ll=0; ll<nlinbins; ll++){
		        	if( bestO > ll+1 ) badLinOdds[ll]++;
		      	}
		    }
		    //Linear odds:
		    for(int ll=0; ll<nlinbins; ll++){
		    	printf(" O > %5i: %.2e\n",ll+1,(double)badLinOdds[ll]/i_iterations);
		    }
		    //Log odds:
		    for(int l=1; l<nlogbins; l++){
		      	if(badOdds[l]==0)break;
		      	printf(" O > %5.5g: %.2e\n",pow(10,l),(double)badOdds[l]/i_iterations);
		      	//XX can add optional computer-readable file write-out here! (with append!
		    }
		    std::cout<<"\n";
  		}
  		//------------------------------

		//Output the parameter estimations
		if(estimate_params){
		  	std::ofstream ofile;
		  	std::string oname="outParameters-";
		  	if(jh!=0) oname+="tp-";
		  	else      oname+="fp-";
		  	oname+=label+".out";
		  	ofile.open(oname.c_str());
		  	ofile<<"# Omax t0 |h| v(km/s) th/PI ph/PI d(km) R a/PI O log(L)\n";
		  	if(jh!=0)
		    ofile<<"# True positive test:\n"
		         <<" "<<int0<<" "<<jh<<" "<<inV<<" "<<inT<<" "<<inP<<" "
		         <<pow(10,inLd)<<" "<<inRdm<<" "<<inadm<<"\n#\n";
		  	for(size_t i=0; i<all_best_params.size(); i++){
		    	for(size_t j=0; j<all_best_params[i].size(); j++){
		      		ofile<<fabs(all_best_params[i][j])<<" ";
		    	}
		    	ofile<<"\n";
		  	}
		  	ofile.close();
		}
	} //end fp for likelihoods

	/* =================================================================
	==================	IF SNR MFT METHOD CHOSEN ======================
	==================================================================*/
	if(method == 1){
  		//Output results to screen & file (false positives)			
  		if (jh == 0 && iRand_bad == 0 && !bad_ref) {
  			std::ofstream fpData("./results/fpData.txt", std::ios_base::app);
  			fpData << "cut\tfp/d\tfp/y\n";
    		std::cout << "\n How many false positives?\n";
    		double cut = 0.5;
    		while (cut < 12.5){
    			cut += 0.5;
    			int numOver = 0;
	      		for (int it = 0; it < i_iterations; it++) {
	        		for (int j = jbeg; j <= jend; j++) {
	          			if (abs(allSNR[it][j]) >= cut) numOver++;
	        		}
	      		}
	    		double avgno = double(numOver) / double(i_iterations);//f.p. in window
	    		double nopd = 2880 * (avgno / double(jend - jbeg + 1)); //f.p. per day
	    		double nopy = nopd * 365.25; //f.p. per year
	    		printf("cut=%2.1f: false pos: %8.3f /day = %8.3f /year\n",
	        			cut, nopd, nopy);
	    		fpData << cut << "\t" << nopd << "\t" << nopy << "\n";
    		}
    		fpData.close();
  		}
  		std::cout << "\n";

		//------------------------------
		//Output results to screen (true positives, & "bad" jumps)
  		int numFound = 0;
		if (jh != 0) {
		    double avgSNR = 0;
		    //int numFound=0;
		    std::cout << std::endl;
		    for (int k = 0; k < i_iterations; k++) {
		    	double bestSNR = 0;
		      	int besti = 0;
		      	for (int i = jbeg; i <= jend; i++) {
		        	if (fabs(allSNR[k][i]) > fabs(bestSNR) ) {
		        		bestSNR = allSNR[k][i];
		          		besti = i;
		        	}
		      	}
		      	int j0 = (int)ceil(injevent[k][0]);
		      	printf("%i: %.2fns event, at t0=%.2f. ", k + 1, injevent[k][1],
		       			injevent[k][0]);
		      	printf("v,th,ph=%.1f, %.1f,%.1fpi. ", injevent[k][2], injevent[k][3],
		      			injevent[k][4]);
		      	if (fabs(allSNR[k][j0]) > SNRcut || fabs(allSNR[k][j0 - 1]) > SNRcut
		        		|| fabs(allSNR[k][j0 + 1]) > SNRcut)
		      	{
		        	std::cout << "--> Found! at: " << j0 << ", SNR=" << allSNR[k][j0] << "   :)\n";       std::cout 
		        			  << "Largest SNR in time window was "<<bestSNR<<" at epoch"<<besti<<std::endl;
		        	numFound++;
		      	}
		      	else {
		        	std::cout.precision(2);
		        	std::cout << "--> NOT found. *******  :(" 
              << "Largest SNR in time window was "<<bestSNR
              <<" at epoch"<<besti<<std::endl;
		      	}
		      	//Find average of the maximum snr values near j0 from the iterations
		       	//if(abs(besti-j0) <=1){
		    	avgSNR += bestSNR;
		      	//}
		      	//double tempSNR = 0;
		      	//for (int jj = j0 - 1; jj <= j0 + 1; jj++) {
		        //	if (allSNR[k][jj] > tempSNR) {
		        //  	tempSNR = allSNR[k][jj];
		        //	}
		      	//}
		      	//avgSNR += tempSNR;
		    }//End loop over iterations
		    
		    avgSNR /= i_iterations;
		    std::cout << std::endl;
		    printf("Found %i/%i = %.1f%% of the ~%.3fsig events.\n", numFound,
		      		i_iterations, double(numFound * 100) / double(i_iterations), jh / avsig);
		    if (iRand_bad != 0 || bad_ref) {
		      	std::cout << "-> Injected BAD ~" << fabs(bad_amp) / avsig << " sigma signals: ";
		      	if (bad_ref)std::cout << "a ref jump";
		      	if (bad_ref&&iRand_bad > 0)std::cout << ", and ";
		      	if (iRand_bad > 0)std::cout << iRand_bad << " random jumps (per clock)";
		      	std::cout << std::endl;
		    }
		    printf("Average SNR:     %.4f\n", avgSNR);
		    std::cout << std::endl;
		}


    	//Just "bad" jumps (only when no signal injected!):
		if (jh == 0 && (iRand_bad != 0 || bad_ref)) {
		    //Print out jump parameters!
		    std::cout << "-> Injected BAD ~" << fabs(bad_amp) / avsig << " sigma signals: ";
		    if (bad_ref)std::cout << "a ref jump";
		    if (bad_ref&&iRand_bad > 0)std::cout << ", and ";
		    if (iRand_bad > 0)std::cout << iRand_bad << " random jumps (per clock)";
		    std::cout << "\nFalse positive rate: \n";
		    const int nlogbins = 10; // 1, 10, 100, 1000...
		    const int nlinbins = 9; // 1,2,3,4,5,...
		    int badSNR[nlogbins] = { 0 };
		    int badLinSNR[nlinbins] = { 0 }; //2 .. 9
		    for (int k = 0; k < i_iterations; k++) {
		      	int j0 = (int)ceil(injevent[k][0]);
		      	double bestO = 0.1;
		      	//Go over entire "window", count spikes!
		      	for (int j = j0 - 10; j <= j0 + 10; j++) {
		        	if (allSNR[k][j] > bestO)
		          		bestO = allSNR[k][j];
		      	}
		      	double lO = log10(bestO);
		      	if (lO <= 0) continue;
		      	int ilO = (int)lO;
		      	if (ilO >= nlogbins) ilO = nlogbins - 1;
		      	for (int l = 0; l <= ilO; l++) badSNR[l]++; //sum cumulatively!
		      	//"linear" SNR (1-9)
		      	for (int ll = 0; ll < nlinbins; ll++) {
		        	if (bestO > ll + 1) badLinSNR[ll]++;
		      	}
		    }
		    //Linear SNR:
		    for (int ll = 0; ll < nlinbins; ll++) {
		      	printf(" O > %5i: %.2e\n", ll + 1, (double)badLinSNR[ll] / i_iterations);
		    }
		    std::cout << "\n";
		}
	}


/*===================================================================================
================================= PRINT OUTPUT ======================================
===================================================================================*/
  	printParameters(numClks,sref,icov,use_ccc,sub_w_mean,in_receivers,sta_sd,
                  simulate_ccc,
                  ref_sig,irandSVN,sig,idiff,jh,bad_ref,iRand_bad,bad_amp,NMC,
                  h0,oddscut, SNRcut, label, method);
  

  	//Output the timing results:
  	time (&end);
  	double dif = difftime (end,start);
  	double difMarg=(double)margTime*1e-6;  //*1000;
  	difMarg/=i_threads;//XX?? not quite correct
  	if(b_it_para)difMarg/=i_threads;//XXX again???
  	std::cout<<"numThreads="<<i_threads<<std::endl;
  	printf ("Total elasped time: %.0f s.\n", dif );
  	printf ("--Per iteration: %.1f s.\n", dif/i_iterations );
  	int numEpochs=jend-jbeg+1;
  	printf ("--Per epoch (total): %2.3f ms.\n",
          	dif*1000./i_iterations/numEpochs);
 	 printf("--Per epoch (Integ): %2.3f ms.\n",difMarg/i_iterations/numEpochs);


return 0;

}








//******************************************************************************
//int plotOdds(double odds[MAXITS][3][2880], int jbegin, int jend,
//             std::string filename, int i_iterations, int idiff, double h0,
//             double dv,
//             int nth, int nt0, double jh, double sig, std::string ClkBLK[7],
//             std::string sref, int numClks[8], bool white)
///*

//XXX 170727 - won't work (changed ClkBlk....

//*/
//{

////  if(jend<=0||jend>iEpochs){//Plot ALL epochs!
////    jend=iEpochs;
////  }
////  if(jbegin<0||jbegin>=iEpochs){//Plot ALL epochs!
////    jbegin=0;
////  }
////
////
////  //std::ofstream gnuPlot;
////

////  execute("mkdir -p output-TestLikelihoods");
////  std::string plotfile="output-TestLikelihoods/"+filename+".out";
////  std::ofstream pfile;
////  pfile.open (plotfile.c_str());
////
////  pfile<<"# "<<filename<<". "<<i_iterations<<" iterations."<<std::endl;

////  pfile<<"# Used "<<idiff<<"-order differencing. "<<std::endl;
////  pfile<<"# Injected a "<<jh<<"ns event"<<std::endl;
////  pfile<<"# h0="<<h0<<"ns, "
////       <<"dv="<<dv<<" km/s, "
////       <<"nth="<<nth<<", "
////       <<"nt0="<<nt0<<"."<<std::endl;
////  if(white){
////    pfile<<"# Used white data: sigma="<<sig<<" ns."<<std::endl;
////  }else{
////    pfile<<"# Ref="<<sref<<".";
////    for(int iType=0;iType<8;iType++){
////      if(numClks[iType]!=0)
////        pfile<<" "<<numClks[iType]<<" "<<ClkBLK[iType];
////    }
////    pfile<<"."<<std::endl;
////  }
////

////  //Writes the data to be plotted
////  for(int it=0;it<i_iterations;it++){
////    for(int j=jbegin;j<=jend;j++){
////      pfile<<j;
////      for(int k=0;k<3;k++){
////        pfile<<" "<<odds[it][k][j];
////      }
////      pfile<<std::endl;
////    }
////    pfile<<" "<<std::endl;
////  }

////  pfile.close();
////
////  std::ofstream gnufile;
////  std::string gnuname="output-TestLikelihoods/temp.gnu";
////  gnufile.open (gnuname.c_str());
////    gnufile<<"set logscale y"<<std::endl;
////    gnufile<<"set format y '%g'"<<std::endl;
////    gnufile<<"set xlabel 'Epochs'"<<std::endl;
////    gnufile<<"set ylabel 'Odds'  offset 2.5"<<std::endl;
////    gnufile<<"set key top right"<<std::endl;
////    gnufile<<"set fit quiet"<<std::endl;
////    gnufile<<"set fit logfile '/dev/null'"<<std::endl;
////    gnufile<<"set border linewidth 2.5"<<std::endl;
////    gnufile<<"set xrange ["<<jbegin-1<<":"<<jend+1<<"]"<<std::endl;
////    gnufile<<"file1='./output-TestLikelihoods/temp.out'"<<std::endl;
////    gnufile<<"plot \\"<<std::endl;
////    gnufile<<"for [i=0 : 25 : 1] file1 using 1:4 every:::i::i "
////           <<"notitle w lp ls i lw 3,\\"<<std::endl;
////    gnufile<<"1,10,100,1000"<<std::endl;
////  gnufile.close();

////  execute("cp output-TestLikelihoods/"+filename
////          +".out output-TestLikelihoods/temp.out");
////  execute("gnuplot output-TestLikelihoods/temp.gnu  -persist > /dev/null 2>&1");



//  return 0;
//}




//******************************************************************************
void printParameters(std::vector<int> numClks, std::string sref, int icov,
    bool use_ccc, bool sub_w_mean, int in_receivers, double sta_sd,
    bool simulate_ccc, double ref_sig,
    int irandSVN, double sig, int idiff, double jh,
    bool bad_ref, int iRand_bad, double bad_amp,
    int NMC, double h0,
    double oddscut, double SNRcut, std::string label, int method)
{


  int temp_clocks=in_receivers;//count the total number of clocks
  for(int iType=0;iType<8;iType++)temp_clocks+=numClks[iType];

  bool white=false; //white noise? Or from PSD
    if(numClks[7]==temp_clocks)white=true;

  //names of each clock/block (only used to format output)
  std::string ClkBlk[7];
    ClkBlk[0]="RbIIF";
    ClkBlk[1]="RbIIA";
    ClkBlk[2]="RbIIR";
    ClkBlk[3]="RbII";
    ClkBlk[4]="CsIIF";
    ClkBlk[5]="CsIIA";
    ClkBlk[6]="CsII";

  std::cout<<"Using "<<temp_clocks<<" clocks:";
  for(int iType=0;iType<7;iType++){ //non-white:
    if(numClks[iType]!=0)
      std::cout<<" "<<numClks[iType]<<" "<<ClkBlk[iType];
  }
  if(numClks[7]!=0){
    printf(" %i White (sig=%.2f)",numClks[7],sig);
  }
  if(in_receivers!=0)
    printf(" %i W-STAs (sig=%.2f)",in_receivers,sta_sd);
  std::cout<<" [with "<<sref<<"].\n";
  if(simulate_ccc){
    std::cout<<"Simulating cross-clock-correlations, with ref_sig = "
             <<ref_sig<<"ns \n";
  }

  if(sub_w_mean)
    std::cout<<"-> Subtracting the weighted mean.\n";

  if(!white){
    if(irandSVN==1){
      std::cout<<"Using random SVNs to generate data."<<std::endl;
    }else{
      std::cout<<"Using averaged PSD (over SVNs) to generate data."<<std::endl;
    }
    if(icov==0)std::cout<<"Not using covariance."<<std::endl;
    if(icov>0)std::cout<<"Multivariate."<<std::endl;
  }
  if(use_ccc)std::cout<<"Including cross-clock-correlations into covariance.\n";

  if(idiff==3){
    std::cout<<"Using mixed differencing (2 for CsII,IIA, 1 for all others).\n";
  }else{
    std::cout<<"Using "<<idiff<<"-order differencing.\n";
  }
  if(jh!=0)std::cout<<"Injecting a "<<jh<<"ns event\n";
  if(jh==0)std::cout<<"Not injecting events (false positive test)\n";

  //"Bad" jumps - robustness test
  if(bad_ref||(iRand_bad>0)){
    std::cout<<"Robust Test:\n --Injecting";
    if(bad_ref) std::cout<<" a BAD reference jump";
    if(iRand_bad>0){
      if(bad_ref) std::cout<<" and";
      std::cout<<" "<<iRand_bad<<" random jumps (per sat)";
    }
    std::cout<<": "<<bad_amp<<"ns\n";
  }

  if(method == 0){
    std::cout<<"Analytic h, hmin="<<h0<<" ns ";
    std::cout<<"with N="<<NMC<<"."<<std::endl;
  
    std::cout<<"Ocut="<<oddscut<<std::endl;
    std::cout<<"Output label="<<label<<std::endl;
    std::cout<<std::endl;
  }
  if(method == 1){
    std::cout << "N=" << NMC << "." << std::endl;

    std::cout << "SNRcut=" << SNRcut << std::endl;
    std::cout << "Output label=" << label << std::endl;
    std::cout << std::endl;
  }

}
