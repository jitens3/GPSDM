/*
PROGRAM:impTwWindow
-Improved "twWindow"

Program searches for "obvious" thin-wall events, using single-differenced data.
Looks for "windows" (in time) in which a certain proprtion of the clocks
experience jumps in both directions of a certain magnitude.
Loops over different time window lengths (J_W), and 'jump' height ranges.
Method described in: http://arxiv.org/abs/1704.06844

Outputs the total number of "events" it found for each window, and jump height.

It requires the 'ACF' files to exist, in order to know the standard deviation
for each clock. If they don't exist, will assume the worst s.d. for each clock.

At the moment, uses the averaged standard deviations that are averaged over
SVNs (from the ACF file).
This could easilt be changes to either
 a) use the actual s.d. from this day, or
 b) Use the averaged s.d. (from ACF file) for the specific SVN..

This program was used for the results in:
http://arxiv.org/abs/1704.06844

Actually, this is a slightly updated version. However, the results didn't
change.


-------------------------------------------------------------------------
CHANGE LOG: Older versions archived in ./versions/ subdirectory
====== CHANGE LOG ======
161115- new program. From "twWindow"
161122- new program again. Start from scratch
161123- fixing "jw" window - before, it assumed that 'ref jump' was in middle!
161231- Lots of fixes. Should work good now.
170103- There was an error - program was taking the s.d. from the H-maser
         reference file, instead of the Rb-Rb (or Cs-Cs) file!! :(
170104- Fixed small error that can occur if acf file doesn't exist
170316- small fixes
170414- Minor update: uses icv file (instead of acf) to read standard deviations
       - Found a "continue" that should have been a break. Actually, won't change
         anything, so not as bad as it sounds.
170505- Fixed input label/extension for ACF files
170630- Found error. Multiplying by the sign of the jump twice, meaning it
         only ever looked for 'positive' jumps, but then counted them each twice.
         Also, reqJDFrac was re-setting to 1
170702- Made arrays static
170731- Now, can be run for just a single block!
170830- Updated to work with class.

Note: should possibly use the actual (calculated) s.d. for given day!
--> or at least, the given svn?? XXX

XXX Would be nice to store +ve / -ve jumps seperately as well!
-------------------------------------------------------------------------

*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#include "JplGpsDataClass.h"

//******************************************************************************
//******************************************************************************
int main (void)
/*
Main Program: Loops through all input data in parallel. Calls other functions
*/
{


//Input Parameters:
  std::string path,PSDpath;           //="/home/ben/jpligsac/";
  std::string ACFlabext; //label+extension for input ACF files
  int minweek;           //=1060; first week with 30s data!
  int maxweek;           //=1879;
  int weekday;           //=7;   // 0-6 = Sun-Mon. 7=whole week.
  std::string label;          //label for the output file
  //Smallest fraction of clocks that must jump in window:
  double reqJmpFrac,reqJDFrac;
  std::string clockType,blockType;//="Rb";
   double xsig;
   double dh; //jump height step-size
   int minClocks;
   int iPoly,iOrder;
   int iMaxWindow;
   double minh,maxh;
   //int iuse_av_svn; //use the average svn (from ACF file), or calculated?
     // bool use_av_svn=true;

  //Read in the input file "impTwWindow.dat"
  std::string junk;
  std::ifstream fInput;
  fInput.open ("impTwWindow.dat");
    fInput >> path;                        getline(fInput,junk);
    fInput >> PSDpath;                     getline(fInput,junk);
    fInput >> ACFlabext;                   getline(fInput,junk);
    fInput >> minweek >> maxweek;          getline(fInput,junk);
    fInput >> weekday;                     getline(fInput,junk);
    fInput >> minClocks;                   getline(fInput,junk);
    fInput >> reqJmpFrac;                  getline(fInput,junk);
    fInput >> reqJDFrac;                   getline(fInput,junk);
    fInput >> dh;                          getline(fInput,junk);
    fInput >> minh >> maxh;                getline(fInput,junk);
    fInput >> iMaxWindow;                  getline(fInput,junk);
    fInput >> xsig;                        getline(fInput,junk);
    fInput >> iPoly >> iOrder;             getline(fInput,junk);
    fInput >> clockType >> blockType;      getline(fInput,junk);
    //fInput >> iuse_av_svn;                 getline(fInput,junk);
    fInput >> label;                       getline(fInput,junk);
  fInput.close();

  // Formats input, & reduces liklihood of errors
  int minday,maxday;
  if(weekday<0||weekday>7){weekday=7;}
  if(minweek<1060||minweek>2000){minweek=1060;}
  if(maxweek>2000||maxweek<1060){maxweek=2000;}
  if(weekday==7){
    minday=0;
    maxday=7;
  }else{
    minday=weekday;
    maxday=weekday+1;
  }
  if(reqJmpFrac>1||reqJmpFrac<=0)reqJmpFrac=0.6; //default % required to 'jump'


  std::cout<<"\n          ###########################################";
  std::cout<<"\n                 Looking for 'potential events'      ";
  std::cout<<"\n          ###########################################\n\n";

  {//Output params to screen:
    printf("Searching the 30s clock data, for weeks %i -> %i\n",minweek,
            maxweek);
    if(weekday!=7)printf(" for weekday %i\n",weekday);
    printf(" Considering days with >= %i %s clocks.\n",minClocks,
            clockType.c_str());
    printf(" Requiring >= %.2f ref jumps, and %.2f jump backs\n",reqJmpFrac,
           reqJDFrac);
    printf(" Jump range: %.3f -> %.3f ns (steps of: %.3f).\n",minh,maxh,dh);
    printf(" with h=h +/- %.4f sigma.\n",xsig);
    double vmin=double(50000/iMaxWindow)/30.;
    printf(" Max window size of %i; => vmin=%.1f km/s\n",iMaxWindow,vmin);
    if(iPoly==1){
      printf("Fitting and removing a %i-order polynomial.",iOrder);
    }else{
      std::cout<<"Not detrending the data. \n";
    }
    printf("~~~~~\n\n");
  }

  int iConly=2; //=0 skip Cs. =1 skip Rb. =2 do both!
  if(clockType=="Rb")iConly=0;
  if(clockType=="Cs")iConly=1;

//// Read the PRN<->SVN map file
// std::string sPRNmap[PRNLINES][PRNPTS];
// readSVNPRN(path,sPRNmap,1);

////////////////////////////////////////////////////////////////////////////////
// MAIN LOOP:
  time_t start,mid,end;
  time (&start);

  //count the clocks and days of data used
  long int daysUsed[2]; //0=Cs, 1=Rb
    daysUsed[0]=0;
    daysUsed[1]=0;
  long int clksUsed[2]; //0=Cs, 1=Rb
    clksUsed[0]=0;
    clksUsed[1]=0;

  //number of h (signal magnitude, S(1) steps)
  int ihmax=int((maxh-minh-dh)/dh)+1;


  //window size things:
  //Increases window size by 1 epoch at first, then by more..
  iMaxWindow++; //increment by one?
  int iWResCut=50;
  int iNumWindows= (int)ceil((iMaxWindow+24*iWResCut-1)/50);
  int iFineWin=(int)ceil(iWResCut/2); //place where window resolution drops (x2)

  const int MAXWINDOWS=100;
  const int MAXHMAX=500;
  if(iNumWindows>MAXWINDOWS||ihmax>MAXHMAX){
    std::cout<<"FAILURE 171: too many windows or h steps!"<<std::endl;
    return 1;
  }

  //array to store number of events!
  static long int nEvents[2][MAXWINDOWS][MAXHMAX]={0};

  //store h's and J_W's:
  std::vector<double> hlist(ihmax); //store actual h's used
  std::vector<int> jwlist(iNumWindows); //store the actual window lengths used

  //Define the sat blocks/clocks
  std::string sBlk[4]={"II","IIA","IIR","IIF"};
  std::string sClk[2]={"Cs","Rb"};

  //MAIN LOOP: loop over all data files
  //loop over weeks/days:
  for(int iW=minweek;iW<maxweek+1;iW++){  //loops through weeks
    for(int iD=minday;iD<maxday;iD++){//loops through each day for given week

      JplGpsData data;
      int jplok=data.readJplBiasData(path,iW,iD);
      if(jplok==2)return 1;//DIR doesn't exist
      if(jplok==1)continue;//this file doesn't exist

      std::cout<<"Reading: jpl"<<data.week+data.day<<" ("<<data.date
               <<". Ref: "<<data.refname<<")\n";

      //Do analysis for Rb and Cs clocks seperately:
      for(int iC=0;iC<2;iC++){//loop over Cs (iC=0) and Rb (iC=1)
        if(iC==iConly)continue; //skip Rb/Cs or not

        //Transfer data into temp "dataCopy" - so can change ref clock twice!
        //ineficient..but easier than re-writing.
        //from here down, only use tData, not dData!
        JplGpsData dataCopy;
        dataCopy.makeACopyOf(data);

        //remove the polynomial (only for testing)
        if(iPoly==1)dataCopy.polynomialDetrend(iOrder);

        //difference the data [s^(0) -> S^(1)]
        dataCopy.differenceData(1); //perform 1st-order differenceing.

        //Change the reference clock!
        bool by_svn=false; //false, so actually take best s.d.!
        if(!by_svn)dataCopy.calculateStdDev();
        dataCopy.swapReference(sClk[iC],blockType,by_svn);

        std::cout<<" Ref --> "<<dataCopy.refname<<"\n";

        //Open and read-in the various "standard deviations"
        //From the "ACF" file
        //XXX can easily update this to use SVN s.d.
        //XXX OR can use actual s.d. of this day!
        double clkSD[4];
        for(int iB=0;iB<4;iB++){//Block
          if((sBlk[iB]!=blockType)&&(blockType!="all"))continue;
          std::string ACFfileName=PSDpath+sClk[iC]+sBlk[iB]+"-av-"
                             +dataCopy.refclk
                             +dataCopy.refblk+"-"+ACFlabext;
          if(MSC_fexists(ACFfileName)){//"CsIIR" doesn't exist
            std::ifstream acffile;
            acffile.open (ACFfileName.c_str());
            std::string sLine;
            while(getline(acffile,sLine)){
              std::stringstream sTest(sLine);
              std::string sok;
              sTest>>sok;
              if(sok=="#")continue;
              double sd1,sd2;
              std::stringstream ssin(sLine);
              ssin>>sd1>>sd2;
              clkSD[iB]=sd1;
              break;
            }
          }else{
            if(sClk[iC]+sBlk[iB]!="CsIIR"){
              //std::cout<<"WARNING: "<<ACFfileName<<" Doesn't exist!\n";
              //if the file doesn't exist, still run, but assume worst error
              if(iC==1)clkSD[iB]=0.1056; //Rb IIR vs IIR (0.0973+0.0083)
              if(iC==0)clkSD[iB]=0.1333; //Cs IIA vs IIF (0.1296+0.0037)
            }
          }
        }//end loop over blocks


        //count the actual number of Cs/Rb clocks
        int numClocks=0;
        for(int i=dataCopy.num_receivers;i<dataCopy.num_clocks;i++){
          if((dataCopy.blk[i]!=blockType)&&(blockType!="all"))continue;
          if(dataCopy.clk[i]==sClk[iC])numClocks++;
        }
        if(numClocks<minClocks)continue; //not enough clocks, don't use
        daysUsed[iC]++; //this day is good, so count it.
        clksUsed[iC]+=numClocks; //"" ""


        // ******* ACTUAL ANALYSIS STARTS HERE *******
        //-Loop through each S1 (+ve and -ve)
        //-Look for "ref clock jumps" (verticle blue line)

        double hh=minh;//initialise
        for(int ih=0;ih<ihmax;ih++){//loop over jump heights [i.e. S^(1)]
          hh+=dh;
          hlist[ih]=hh; //store which h's we tested
          for(int pm=-1;pm<=1;pm+=2){//loop over +ve/-ve 'jumps'
            //nb: hh always +ve. pm*hh = actual h
            //look for "vertical" line:
            for(int j=1;j<dataCopy.num_epochs;j++){//loop over 'ref jump' epochs:
              int nCHJ=0; //number of clocks that "jumped" here
              //to remember which clocks jumped:
              std::vector<bool>bTCJ(dataCopy.num_clocks);
              //loop over only the satellite clocks
              for(int i=dataCopy.num_receivers;i<dataCopy.num_clocks;i++){
                bTCJ[i]=false;
                if(dataCopy.blk[i]!=blockType&&blockType!="all")continue;
                if(dataCopy.clk[i]!=sClk[iC])continue;//skip wrong clocks
                //identify the block:
                int iB=-1;
                for(int ib=0;ib<4;ib++){
                  if(dataCopy.blk[i]==sBlk[ib]){
                    iB=ib;
                    break;
                  }
                }//END identify the block:
                if(iB==-1)std::cout<<"FAILURE 303: didn't find block?\n";
                double xs=xsig*clkSD[iB];//"range" for jump
                // XXX here! or dataCopy.sdev??^^
                double dd=dataCopy.bias[i][j];
                if((dd>(-pm*hh-xs))&&(dd<(-pm*hh+xs))){
                  nCHJ++;
                  bTCJ[i]=true;//remember which clocks jumped (at ref jump)
                }
              }//END initial loop over clocks
              //check to see if enough of the clocks jumped:
              double fracHJ=double(nCHJ)/double(numClocks);//fraction that jumped
              if(fracHJ<reqJmpFrac)continue; //didn't find H jump. skip rest
              //If fracHJ>=reqJmpFrac, then found ref jump!
              //Found ref jump: Check if clocks jump "back down"
              #pragma omp parallel for
              for(int ijw=0;ijw<iNumWindows;ijw++){//loop over window lengths
                //'jw' is the length of the window before/after ref-clock epoch
                //i.e., half the window length!
                int jw;
                if(ijw<iFineWin){
                  //for small windows, increase 1 epoch at a time
                  jw=ijw+1;
                }else if(ijw>=iFineWin){
                  //for small windows, increase 25 epochs at a time
                  jw=iFineWin+(ijw-iFineWin+1)*25;
                }
                jwlist[ijw]=2*jw+1;  //store the actual window lengths
                int nCJB=0; //number of clocks that "jumped" here
                //step the window across
                for(int j0=j-jw; j0<=j+jw; j0++){//step the window across
                 //nb: 'j0' is center of window. 'j' is where ref clock jumped
                  //loop over clocks (satellite only) for the window
                  for(int i=dataCopy.num_receivers;i<dataCopy.num_clocks;i++){
                    if(!bTCJ[i])continue;//only check clocks for which ref jumped
                    if(dataCopy.blk[i]!=blockType&&blockType!="all")continue;
                    if(dataCopy.clk[i]!=sClk[iC])continue;//skip Rb/Cs clocks
                    //identify which 'block' this clock belongs to (for sigma!)
                    int iB=-1;
                    for(int ib=0;ib<4;ib++){
                      if(dataCopy.blk[i]==sBlk[ib]){
                        iB=ib;
                        break;
                      }
                    }//END identify the block:
                    if(iB==-1)std::cout<<"FAILURE 344: didn't find block?\n";
                    //j0 is center of window, j is ref jump, jw is half-window
                    for(int jj=j0-jw;jj<=j0+jw;jj++){//loop over j within window
                      if((jj<=1)||(jj>=dataCopy.num_epochs))continue;
                      if(jj==j)continue; //don't re-check the ref jump!
                      double xs=xsig*clkSD[iB];//"range" for jump
                      // XXX here! or dataCopy.sdev??^^
                      double dd=dataCopy.bias[i][jj];
                      if((dd>(pm*hh-xs))&&(dd<(pm*hh+xs))){//check if jumped
                        nCJB++;
                        break; //found one, don't re-count it, go to next clock!
                      }//END check if jumped
                    }//END loop over j within window!
                  }//END loop over clocks for the window
                  if(nCJB>=reqJDFrac*nCHJ){//these clocks jumped up/down!
                    nEvents[iC][ijw][ih]++;
                    break; //stop stepping window along
                  }
                }//END step window across
              }//std::endloop over windows
            }//END loop over 'ref jump' epochs:
          }//END loop over +ve/-ve 'jumps'
        }//END loop over jump heights

      }//END loop over Cs (iC=0) and Rb (iC=1)

    }//END loop over days[j]


    //print out time-remaining:
    double perc=(double(iW-minweek+1)/double(maxweek-minweek+1))*100.;
    time (&mid);
    double tsofar = difftime (mid,start);
    double tleft=tsofar*(100/perc)-tsofar;
    int dd=int(tleft/86400);
    int hh=int((tleft-dd*86400)/3600);
    int mm=int((tleft-dd*86400-hh*3600)/60);
    int ss=int(tleft-dd*86400-hh*3600-mm*60);
    printf("# %.1f%% done: ",perc);
    if(dd>0)printf("%i days ",dd);
    if(hh>0)printf("%i hours ",hh);
    if(mm>0)printf("%i mins ",mm);
    printf("%i secs remaining.\n",ss);

  }//END loop over weeks[iW]

  //write out results!
  std::cout<<"\nSnapshot of results:\n";
  for(int iC=0;iC<2;iC++){
    if(iC==iConly)continue; //skip Rb/Cs or not
    std::cout<<"For "<<sClk[iC]<<" ("<<blockType<<"):"<<std::endl<<"h\\jw: ";
    for(int ijw=0;ijw<20;ijw+=2){
      printf("%4i ",jwlist[ijw]);
    }
    std::cout<<"\n-------------------------------------------------------\n";
    for(int ih=0;ih<ihmax;ih+=5){
      printf("%.2f| ",hlist[ih]);
      for(int ijw=0;ijw<10;ijw++){
        if(nEvents[iC][ijw][ih]<=9999){
          printf("%4ld ",nEvents[iC][ijw][ih]);
        }else{
          std::cout<<"9999+";
        }
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  }

  for(int iC=0;iC<2;iC++){
    if(iC==iConly)continue; //skip Rb/Cs or not
    std::string filename="events-"+sClk[iC]+blockType+"-"
                         +std::to_string(minweek)+"-"
                         +std::to_string(maxweek)+"-"+label+".txt";
    std::ofstream oFile;
    oFile.open (filename.c_str());
    oFile<<minweek<<" "<<maxweek<<std::endl;
    oFile<<double(daysUsed[iC])/365.25<<" yrs, "<<clksUsed[iC]<<" clk-days"
         <<std::endl;
    oFile<<minClocks<<" "<<reqJmpFrac<<" "<<reqJDFrac<<" "<<dh<<" "<<xsig<<" "
         <<std::endl;
    for(int ijw=0;ijw<iNumWindows;ijw++){
      oFile<<jwlist[ijw]<<" ";
    }
    oFile<<std::endl;
    for(int ih=0;ih<ihmax;ih++){
      oFile<<hlist[ih]<<" ";
    }
    oFile<<std::endl;
    for(int ih=0;ih<ihmax;ih++){
      for(int ijw=0;ijw<iNumWindows;ijw++){
        oFile<<nEvents[iC][ijw][ih]<<" ";
      }
      oFile<<std::endl;
    }
    oFile.close();
  }

////////////////////////////////////////////////////////////////////////////////

  time (&end);
  double dif = difftime (end,start);
  printf ("Elasped time is %.0f seconds.\n", dif );
  return 0;
}
//******************************************************************************
//******************************************************************************
