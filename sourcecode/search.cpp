/*
 PROGRAM: search

 Still under developement.

 This program will search through archival gps data file in compute 
 likelihoods and odds ratio values and snr values.

 See search.dat for options and features.


 =======Change log=====
 ...

 ****** To Do ******

 */


#include "SearchClass.h"
#include <string>
#include <time.h>
#include <fstream>
#include <iostream>
#include <vector>


struct options 
{
	std::string path_to_data;
	std::string path_to_psd; //location of input PSD files
	std::string path_to_acf; //location of input ACF files (often the same as PSD)
	std::string path_to_cdf; //location of input numeric prior CDF files
	std::string acf_label; //which PSD/ACF files (which label) to use
	int sanity;
  	int binary;
	int startday, startweek, endday, endweek; //The range of week+day combination
	int invest, event, hJW;
	std::string model;  //model (wall/string/monopole)
	double iRgps;
	int ialpha;
	std::string coupling; //assume specific coupling?
	int iexact; //difference (1 or 2) and exact vs approx. (0 or 1)
	int iccc=1;
	int icov=1;
	int NMC;
  	double hmin;
	double oddscut, snrcut;
	int iestimate;
	int iestNMC;
	int ilimit_cores, iwhichpara;
	bool detrend, idiff;
	std::string label;
	bool ibias;
} input ;


void readIn(options & obj);

void errorCheck(options & input);

int main(){
    /****************************************************************************
    ************************** ESTABLISHING options ***************************
    ****************************************************************************/

   	//for timing:
	
  //for timing:
  time_t start, end; //total program time
  time (&start);
  long int margTime=0;  //time just the integrations
 
  //options v;
  readIn(input);
  errorCheck(input);

  std::string path_to_avgccc = "./all_network_average_ccc.txt";


    /****************************************************************************
    ************************** LOOPING ITERATIONS *******************************
    ****************************************************************************/

	  for(int iw = input.startweek; iw <= input.endweek; iw++){ // loop through weeks
		for(int id = 0; id < 7; id++){ // loop through days in week
			if(iw == input.startweek && id == 0) id = input.startday; // if start day, don't start at id=0

				std::string wd = std::to_string(iw) + std::to_string(id);
  				const std::string path = "mkdir ./results/"+wd;
				const char * newpath = path.c_str();

  				int open = system(newpath);
                
            std::cout << "Here in search 1" << "\n";
  				//establish SearchClass object
  				SearchClass statistic(input.path_to_data, input.NMC,
			 	 					    iw, id, input.detrend, 
									    input.idiff, input.oddscut, 
										input.snrcut, input.hmin);
            
            std::cout << "Here in search 2" << "\n";
            
  				if(input.invest){
  					statistic.m_investigate=1;
  					statistic.m_estNMC=input.iestNMC;
  				}

				//std::cout << "after construct\n";
				statistic.verbose = true; //yell at us if something is wrong

				std::string path_v_cdf="./priorCDF/vperp.cdf";    //Define path to prior distributions
				std::string path_psi_cdf="./priorCDF/psiperp.cdf";

				NumericCdfInverse v_prior(path_v_cdf);    //Generate prior distributions
				NumericCdfInverse psi_prior(path_psi_cdf);

				// more data processing routines
				statistic.formRelativeCouplings();
            
            std::cout << "this is num clocks to acf: " << statistic.num_clocks << "\n";
            statistic.formAutoCorrelation(input.path_to_acf,input.acf_label);
//            statistic.formAutoCorrelation();		// use this for gpsdays before 12800

/*
		///****************************************************************************
        // ************************* inject signal test *******************************
        // ***************************************************************************
        // added on 12/6/19 to check injected signal
        //Inject Actual event into the clock data:
        std::string model = "twall";

//        double vDM = 100*30.;
//        double inT = 0.236;
//        double inP = 1.741;
//        double indTP = 0.210;
//        double d = 0.;
//        double Rdm = 0.;
//        double adm = 0.;
        //theta angle for DM:
//        double thetaDM = MFS_randDouble(inT*PI-indTP*PI,inT*PI+indTP*PI);
//        if(thetaDM<0) thetaDM = -thetaDM;
//        while(thetaDM>PI) thetaDM = 2*PI - thetaDM;
        // phi angle for DM:
//        double phiDM = MFS_randDouble(inP*PI-indTP*PI,inP*PI+indTP*PI);
        // Unit vector, in direction of DM velocity
//        double eDM[3];
//        eDM[0]=sin(thetaDM)*cos(phiDM);  //x
//        eDM[1]=sin(thetaDM)*sin(phiDM);  //y
//        eDM[2]=cos(thetaDM);             //z
      //  eDM[0]=-0.996654;  //x
      //  eDM[1]=-0.0235658;  //y
      //  eDM[2]=0.0782763;             //z
//        eDM[0]=0.46;  //x
//        eDM[1]=-0.49;  //y
//        eDM[2]=0.74;             //z

        double vDM, eDM[3], d, Rdm, adm;
        double amp = 0.2;    // nanosec
        double t0 = 1200.5;
        int iDGen=1;
        //Store the signal in 2D array:
        double s[statistic.num_clocks][MAXJW]; //zeroing not needed!

	    //Form  rsat and rref arrays:
        //Here, I just take the "slice" of positions that I need to calculate
        // the signals.
        // Note: this is the same as is done in the likelihoods program
        double rsat[statistic.num_clocks][3];
        double rref[3];
        int j0=(int)ceil(t0); //index to take
        int iClocks=statistic.num_clocks;
        //std::cout << "Here before rsat" << "\n";
        for(int i=0; i<iClocks; i++){
           for(int ix=0; ix<3; ix++)
               rsat[i][ix]=statistic.pos[i][j0][ix];
        }
        //std::cout << "Here after rsat" << "\n";
        for(int ix=0; ix<3; ix++)
        	rref[ix]=statistic.refpos[j0][ix];

        //std::cout << "Here after rref" << "\n";
        //generate the relatvant signal:
    	randomParameters(v_prior,psi_prior,model,j0,t0,vDM,eDM,d,Rdm,adm);
        int jmin=int(j0-30);
        int jmax=int(j0+30);
        int JW_gen=jmax-jmin+1; //new!
        generateSignal(model,s,rsat,rref,iClocks,jmin,JW_gen,iDGen,statistic.keff,t0,
                          vDM,eDM,d,Rdm,adm);                  // the problem lies here
        //std::cout << "Here after gensignal" << "\n";
        std::ofstream sigFile1("./results/injSignalTemplate"+std::to_string(10*iw+id)+".txt");
  		for(int a=0;a<statistic.num_clocks;a++){
  		//	if(statistic.clk[a]=="Rb"){
  				for(int b=0;b<MAXJW;b++){
		  			sigFile1<<s[a][b]<<" ";
     					std::cout << s[a][b] << " ";
				}
	  			sigFile1<<std::endl;
 				std::cout <<"\n";
  		//	}
  		}
  		sigFile1.close();
       std::ofstream sigFile2("./results/injSignalParams"+std::to_string(10*iw+id)+".txt");
       sigFile2 << eDM[0] << " "<< eDM[1] << " "<< eDM[2] <<" " << vDM << "\n";
       sigFile2.close();		

    	std::cout << "clocks " << iClocks << "\n";
    	//inject the signal into the data:
        for(int i=0;i<iClocks;i++){
        	for(int j=jmin;j<=jmax;j++){
            	statistic.bias[i][j]+=amp*s[i][j-jmin];
            }
        }

        std::cout << "Here after injsignal" << "\n";
        statistic.differenceData(0);
        //std::cout << "Here after diff3" << "\n";
    	// difference must be zero in JPL constructor
    	// end of check for signal 12/6/19
    	///****************************************************************************
        // ************************* inject signal test *******************************
        // ***************************************************************************
        
        statistic.calculateStdDev();
*/
        //std::cout<<statistic.calculateStdDev();

				statistic.calculateCrossClockCorrelation(61);

				//statistic.avgNetCCC(id,iw, path_to_avgccc);

				//std::ofstream b0val("b0_calc.txt", std::ios_base::app);
				//b0val<< iw+id <<"\t"<<statistic.b0[0]<<std::endl;
				//b0val.close();

				// Search trough the data!!
				std::cout<< "Day: " << iw << id << "-\n";
/*
				// this is just to store the number of clocks in network per day
				std::string clkPath = "./results3/num_clocks__" + wd + ".out";
    				std::ofstream clk(clkPath);
				clk << statistic.num_clocks << "\n";
				
    				clk.close();
				continue;
				// end of storing number of clocks per day
*/
				int JW = MAXJW;
				int jbeg;
				int jend;
				if(input.invest){
					jbeg = input.event - input.hJW;
					jend = input.event + input.hJW;
				}
				else{
					jbeg = JW/2+1;
					jend = statistic.num_epochs - JW/2-2;
				}

// commenting out this portion to make the biasdata file in methods.cpp parameter estimation
/*
				if(input.invest){
					std::string biasOut = "./results/"+wd+"/biasData"+wd+"-"+std::to_string(input.event)+".out";
					std::ofstream biasData(biasOut);
					if(input.invest == 1){
						for(int a=0; a<statistic.num_clocks; a++){
					//	if(statistic.clk[a]=="Rb"){     // commented out 11/3/20 since we already eliminated the Cs clocks
  							for(int b=jbeg; b<=jend; b++){
								biasData<<statistic.bias[a][b]<<" ";
							}
							biasData<<std::endl;
	  				//	}
					}
					biasData.close();
					}
					else if (input.invest == 2){	
						for(int a=0; a<statistic.num_clocks; a++){
							for(int b=jbeg; b<=jend; b++){
									biasData<<statistic.bias[a][b]<<" ";
							}
							biasData<<std::endl;
						}
						biasData.close();
					}
				}

*/
// end of commented out portion for biasdata file on 12/21/21 by tday
				
				statistic.MCSearchStats(input.model, statistic, input.icov, 
             	                	 	JW, jbeg, jend, input.NMC, 
										v_prior, psi_prior, input.iexact, iw, 
										id, input.ilimit_cores, 
										input.sanity,input.event);

            // commented this 11/06/20 for adding chi_thresh into flag tday
				//statistic.serialize(iw,id,input.sanity, input.NMC, jbeg, jend);
            
            // this is for the chi_squared filter
            // added here on 11/06/2020
            double alpha = pow(10,-10);
            double beta = 2*(1-alpha)-1;
            int chi_thresh = (statistic.num_clocks)*MAXJW + sqrt(2)*sqrt(2*(statistic.num_clocks)*MAXJW) * MFS_inverseErf(beta);
                statistic.serialize(iw,id,input.sanity, input.NMC, jbeg, jend, chi_thresh);

			if(iw == input.endweek && id == input.endday) break;
		}
	  }




    //Output the timing results:
  	time (&end);
  	double dif = difftime (end,start);

    return 0;
}

//func def
void readIn(options & obj)
{
	/********************* Read in the input .dat file **************************/
	std::ifstream fInput;
	std::string junk;
	fInput.open("search.dat");
	fInput >> obj.path_to_data;								getline(fInput, junk);
	fInput >> obj.path_to_psd;								getline(fInput, junk);
	fInput >> obj.path_to_acf;								getline(fInput, junk);
	fInput >> obj.path_to_cdf;								getline(fInput, junk);
	fInput >> obj.sanity >> obj.binary;						getline(fInput, junk);
	fInput >> obj.startweek >> obj.startday 
		   >> obj.endweek >> obj.endday;					getline(fInput, junk);
	fInput >> obj.invest >> obj.event >> obj.hJW;			getline(fInput, junk);
	fInput >> obj.acf_label;								getline(fInput, junk);
	fInput >> obj.iexact >>  obj.detrend >> obj.idiff;		getline(fInput, junk);
	fInput >> obj.model >> obj.iRgps >> obj.ialpha;			getline(fInput, junk);
	fInput >> obj.coupling;									getline(fInput, junk);
	fInput >> obj.NMC >> obj.iestNMC >> obj.hmin;			getline(fInput, junk);
	fInput >> obj.oddscut >> obj.snrcut;					getline(fInput, junk);
	fInput >> obj.ilimit_cores >> obj.iwhichpara;			getline(fInput, junk);
	fInput >> obj.label;									getline(fInput, junk);
	fInput.close();
}

void errorCheck(options & input)
{
	std::string start_day = MSC_wdtoymd(input.startweek, input.startday);
	std::string end_day =  MSC_wdtoymd(input.endweek, input.endday);

	std::string exact;
	std::string white_noise;
	std::string model;
	std::string couple;
	if(input.iexact == 1){
		exact = "Exact";
	}else{
		exact = "Perturbative";
	}
	std::string params;
	if(input.iestimate == 1){
		params= "yes";
	}else{
		params = "no";
	}

	std::cout << "\n\n ###########################################"
			  << "\n             Running search program            "
			  << "\n ###########################################\n\n";

	std::cout << "Search Parameters: \n\n"
			  << "Searching days \t\t-\t" << start_day << " through " << end_day;
			if(input.invest){
	std::cout << "\nSearch style \t\t-\tInvestigating epoch " << input.event;
			}
			else{
	std::cout << "\nSearch style \t\t-\tSearching full day";
			}
	std::cout << "\nWindow size \t\t-\t" << MAXJW
			  << "\nMCMC iterations \t-\t" << input.NMC
			  << "\nInverse method \t\t-\t" << exact;
			  

	if (input.iexact != 1) input.iexact = 0;
	//if (input.method != 1) input.method = 0;

	//Make sure input for input.models is ok. Will default back to thin-walls!
	if (input.model == "wall") input.model = "dwall";
	if (input.model == "mono") input.model = "monopole";
	if (input.model != "dwall" && input.model != "string" && input.model != "monopole") input.model = "twall";

	//Makes sure input input.coupling is ok.
	//If not alpha, m_e, or m_p, assumes all clocks the same.
	if (input.coupling == "a") input.coupling = "alpha";
	if (input.coupling == "ep" || input.coupling == "mep" || input.coupling == "e") input.coupling = "me";
	if (input.coupling == "q") input.coupling = "mq";
	if (input.coupling != "alpha" && input.coupling != "me" && input.coupling != "mq") input.coupling = "";

	if(input.model == "dwall"){
		model = "Domain wall";
	}else if(input.model == "twall"){
		model = "Thin wall";
	}else{
		model = input.model;
	}
	if(input.coupling == ""){
		couple = "all";
	}else{
		couple = input.coupling;
	}

	std::cout << "\nModel type \t\t-\t" << model
			  << "\nCoupling \t\t-\t" << couple
			  << "\nh-integral minimum \t-\t" << input.hmin
			  << "\nOdds-ratio cut-off \t-\t" << input.oddscut
			  << "\nSNR cut-off \t\t-\t" << input.snrcut << "\n\n";
			  


	//Use real, or "flat" priors for velocity/angles?
	bool use_flat_priors = true;

	if (input.acf_label == "na")input.acf_label = "";

}
