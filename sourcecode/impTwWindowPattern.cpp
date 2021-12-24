/*
PROGRAM:impTwWindowPattern

Similar to impTwWindow, but this program has only a single window size
(50 epochs), and is designed to output 'pattern' files, for Conner/Wyatt
mathematica "pattern matching" program.
The actual window used by the program is 101 epochs long. But, since we always
give the output with j_r (when the reference clock jumped) in the middle, the
actual maximum window is ~50 epochs, which is still large.

It requires the 'ACF' files to exist, in order to know the standard deviation
for each clock. If they don't exist, will assume the worst s.d. for each clock.
As in "impTwWindow".
They can be downloaded from github (benroberts999/PowerSpectrums)

Format of output file is:
  * date jplFile Clock REF jmin jmax jR
  * S_cut
  * numClocks numEpochs
  * r.n(rproj) rx ry rz sig svn-prn-clk-blk ...data... 
  * [one line for each clock; the ref clock is always the last entry]

At the moment, it doesn't seem to work with openMP.
This probably can be fixed easily, however, it doesn't matter much, since
the bottle-neck is reading in the files, and it's fairly quick.

-------------------------------------------------------------------------
CHANGE LOG: Older versions archived in ./versions/ subdirectory
====== CHANGE LOG ======
170124- new program. From "impTwWindow"
170125- minor fixes, outputs full x,y,z position, AND projection 
170316- small fixes
170414- Minor update: uses icv file instead of acf to get standard deviations.
      - Found a "continue" that should have been a break. Actually, won't change
        anything, so not as bad as it sounds.
170505- Fixed input label for ACF files
170630- Found error. Multiplying by the sign of the jump twice, meaning it
        only ever looked for 'positive' jumps, but then counted them each twice.
        (Same as impTwWindow).
        Also, reqJDFrac was re-setting to 1
170702- Doesn't work with OMP! made arrays static
170726- removed static.
        Added csd (standard deviation) to the output file. NOTE: uses the sd
        that was averaged over all SVN. So only different per-block!
        Updated the output file! Better format!

XXX Note: should possibly use the actual (calculated) s.d. for given day!
--> or at least, the given svn?? XXX

XXX To Do XXX
-Make work with OMP
-Better clock sd's
-------------------------------------------------------------------------

*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//using namespace std;

#include "JplGpsDataClass.h"


//******************************************************************************
//******************************************************************************
int main (void)
/*
Main Program: Loops through all input data in parallel. Calls other functions
*/
{

 //srand(time(NULL));

//Input Parameters:
 std::string path,PSDpath,odir;   //input and output file directories
 std::string ACFlabext;      // Label + extension for input ACF files
 int minweek;           //=1060; first week with 30s data!
 int maxweek;           //=1879;
 int weekday;           //=7;   // 0-6 = Sun-Mon. 7=whole week.
 std::string label;          //label for the output file
 
 double reqJmpFrac,reqJDFrac;//=0.5; //Smallest fraction of clocks that jmp
 std::string clockType;//="Rb";
  double xsig;
  double dh; //jump height step-size
  int minClocks;
  double minh,maxh;
  
 std::string junk;
 std::ifstream fInput;
 fInput.open ("impTwWindowPattern.dat");
   fInput >> path;                        getline(fInput,junk);
   fInput >> PSDpath;                     getline(fInput,junk);
   fInput >> ACFlabext;                   getline(fInput,junk);
   fInput >> odir;                        getline(fInput,junk);
   fInput >> minweek >> maxweek;          getline(fInput,junk);
   fInput >> weekday;                     getline(fInput,junk);
   fInput >> minClocks;                   getline(fInput,junk);
   fInput >> reqJmpFrac;                  getline(fInput,junk);
   fInput >> reqJDFrac;                   getline(fInput,junk);
   fInput >> dh;                          getline(fInput,junk);
   fInput >> minh >> maxh;                getline(fInput,junk);
   fInput >> xsig;                        getline(fInput,junk);
   fInput >> clockType;                   getline(fInput,junk);
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
 if(reqJmpFrac>1||reqJmpFrac<=0)reqJmpFrac=0.6;
 //if(reqJDFrac>reqJmpFrac)reqJDFrac=1;

 std::cout<<"\n          ###########################################";
 std::cout<<"\n              Generating 'potential event' Patterns   ";
 std::cout<<"\n          ###########################################\n\n";
 
 std::string cmd="mkdir -p "+odir;
 MSC_execute(cmd);
 if(!MSC_direxists(odir)){
   std::cout<<"ERROR 136: Directory "<<odir<<" does not exist."
       <<"Please create it, and then try again."<<std::endl;
   return 1;
 }

 {//Output params to screen:
   printf("Searching the 30s clock data, for weeks %i -> %i\n",minweek,maxweek);
   if(weekday!=7)printf(" for weekday %i\n",weekday);
   printf(" Considering days with >= %i %s clocks.\n",minClocks,
           clockType.c_str());
   printf(" Requiring >= %.2f ref jumps, and %.2f jump backs\n",reqJmpFrac,
          reqJDFrac);
   printf(" Jump range: %.3f -> %.3f ns (steps of: %.3f).\n",minh,maxh,dh);
   printf(" with h=h +/- %.4f sigma.\n",xsig);
   printf("~~~~~\n\n");
 }
 
 int iConly=2; //=0 skip Cs. =1 skip Rb.
 if(clockType=="Rb"){
    iConly=0;
 }else if(clockType=="Cs"){
    iConly=1;
 }else{
   std::cout<<"Invalid clocktype: "<<clockType<<std::endl;
   return 1;
 }

//// Read the PRN<->SVN map file
// std::string sPRNmap[PRNLINES][PRNPTS];
// readSVNPRN(path,sPRNmap,1);
 
 double galvec[3];
 {
   const double galvecx=0.46332,galvecy=-0.49003,galvecz=0.73838;
   galvec[0]=galvecx;
   galvec[1]=galvecy;
   galvec[2]=galvecz;
 }
 
 //array to store 1D projections of sat. positions
 std::vector<double> projpos;

////////////////////////////////////////////////////////////////////////////////             
// MAIN LOOP:
 time_t start,mid,end;
 time (&start);
 
   
 int ihmax=int((maxh-minh-dh)/dh)+1;
 
 //Define the sat blocks/clocks
 std::string sBlk[4]={"II","IIA","IIR","IIF"};
 std::string sClk[2]={"Cs","Rb"};
       
 
 bool wpara=false;
 if((maxweek-minweek)<3)wpara=false;
 
 //loop over weeks/days:
 //#pragma omp parallel for if(wpara)
 for(int iW=minweek;iW<=maxweek;iW++){  //loops through weeks
   std::string week = std::to_string(iW);
   for(int iD=minday;iD<maxday;iD++){//loops through each day for given week 
     MSC_progressBar("Searching jpl"+week,35,iD,maxday);
     
     //create JPL data object:
     JplGpsData data;
     
     //read in the JPL bias and the ECI files:
     data.readJplBiasData(path,iW,iD);
     data.readEciPos("sat"); //only need sat. positions!

     //Start the analysis:
     for(int iC=0;iC<2;iC++){//loop over Cs (iC=0) and Rb (iC=1) XXX not needed!
       if(iC==iConly)continue; //skip Rb/Cs or not
       
       //difference the data [s^(0) -> S^(1)]
       data.differenceData(1); //perform 1st-order differenceing.
       
       //calculate the clock s.d.'s 
       //These are used for output only; "average" sds used in search.
       data.calculateStdDev();
       
       //Change the reference clock:
       data.swapReference(sClk[iC]); //uses "best s.d.", USE ANY BLOCK:

       //Open and read-in the various "standard deviations" from the ACF files
       std::string acfFN=PSDpath;
       double clkSD[4]; //holds the s.d. XXX For now, just each block!
       for(int iB=0;iB<4;iB++){//Block
         std::string ACFfileName=acfFN+sClk[iC]+sBlk[iB]+"-av-"+data.refclk
                            +data.refblk+"-"+ACFlabext;
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
             //std::cout<<ACFfileName<<" "<<clkSD[iB]<<std::endl;
             break;
           }
         }else{
           if(sClk[iC]+sBlk[iB]!="CsIIR"){
             std::cout<<"WARNING 213: "<<ACFfileName<<" Doesn't exist!\n";
             clkSD[iB]=0.125; //"worst" possible error..
           }
         }
       }
       
       //count number of Cs/Rb clocks
       int numClocks=0;
       for(int i=data.num_receivers;i<data.num_clocks;i++){
         if(data.clk[i]==sClk[iC])numClocks++;
       }
       if(numClocks<minClocks)continue; //not enough clocks, don't use
       
       // *** Do the actual data analysis: look for potential events *** 
       double hh=minh;//initialise
       for(int ih=0;ih<ihmax;ih++){//loop over jump heights
         hh=minh+ih*dh;
         for(int pm=-1;pm<=1;pm+=2){//loop over +ve/-ve 'jumps'
           //double h=pm*hh;
           //look for "horizontal line!
           for(int j=1;j<data.num_epochs;j++){//loop over 'ref jump' epochs:
             int nCHJ=0; //number of clocks that "jumped" here
             //to remember which clocks jumped:
             std::vector<bool>bTCJ(data.num_clocks);
             for(int i=data.num_receivers;i<data.num_clocks;i++){//loop over clocks
               bTCJ[i]=false;
               if(data.clk[i]!=sClk[iC])continue;//skip wrong clocks
               //identify the block:
               int iB=-1;
               for(int ib=0;ib<4;ib++){
                 if(data.blk[i]==sBlk[ib]){
                   iB=ib;
                   //continue;
                   break;
                 }
               }//END identify the block:
               if(iB==-1)std::cout<<"FAILURE 276: didn't find block?\n";
               double xs=xsig*clkSD[iB];//"range" for jump
               double dd=data.bias[i][j];
               if((dd>(-pm*hh-xs))&&(dd<(-pm*hh+xs))){
                 nCHJ++;
                 bTCJ[i]=true;//remember which clocks jumped (at ref jump)
               }
             }//END initial loop over clocks
             double fracHJ=double(nCHJ)/double(numClocks);//fraction that jumped
             if(fracHJ<reqJmpFrac)continue; //didn't find H jump. skip rest
             //Found ref jump: Check if clocks jump "back down"
             //'jw' is the length of the window before/after ref-clock epoch
             //i.e., half the window length!
             int jw=50; 
             int nCJB=0; //number of clocks that "jumped" here
             //int j0=j;
             for(int j0=j-jw; j0<=j+jw; j0++){//step window across
              //nb: 'j0' is center of window. 'j' is where ref clock jumped
              //loop over clocks for the window:
               for(int i=data.num_receivers;i<data.num_clocks;i++){
                 if(!bTCJ[i])continue;//only check clocks for which ref jumped
                 if(data.clk[i]!=sClk[iC])continue;//skip Rb/Cs clocks 
                 
                 //identify which 'block' this clock belongs to
                 int iB=-1;
                 for(int ib=0;ib<4;ib++){
                   if(data.blk[i]==sBlk[ib]){
                     iB=ib;
                     //continue;
                     break;
                   }
                 }//END identify the block:
                 
                 if(iB==-1)std::cout<<"FAILURE 305: didn't find block?\n";
                 //j0 is center of window, j is ref jump, jw is half-window
                 for(int jj=j0-jw;jj<=j0+jw;jj++){//loop over j within window!
                   if((jj<=1)||(jj>=data.num_epochs)||(jj==j))continue;
                   double xs=xsig*clkSD[iB];//"range" for jump
                   double dd=data.bias[i][jj];
                   if((dd>(pm*hh-xs))&&(dd<(pm*hh+xs))){//check if jumped
                     nCJB++;
                     break; //found one, don't re-count it, go to next clock!
                   }//END check if jumped
                 }//END loop over j within window!
                 
               }//END loop over clocks for the window
               if(nCJB>=reqJDFrac*nCHJ){//these clocks jumped up/down!
                 //HERE! FOUND a "potential event" but only once!? 
                 
                 //Project each clock position onto ECI gal. vector
                 projpos.resize(data.num_clocks);
                 for(int i=data.num_receivers;i<data.num_clocks;i++){
                   double temp=0;
                   for(int x=0;x<3;x++){
                     //temp+=galvec[x]*clkpos[i][j][x];
                     temp+=galvec[x]*data.pos[i][j][x];
                   }
                   projpos[i]=temp;
                 }
                 
                 //write out window to file
                 std::ofstream ofile;
                 std::string shh=std::to_string(int(hh*100));
                 std::string fname=odir+"/pat"+data.week+data.day+"-"
                    +std::to_string(j)+"-"+shh+"-"+label+".txt";
                 ofile.open (fname.c_str());
                 int iBeg=j-jw;
                 int iEnd=j+jw;
                 std::string sREF;
                 if(data.refblk=="AR"){
                   sREF=data.refprn;
                 }else{
                     sREF=data.refsvn+"-"+data.refprn+"-"+data.refclk
                         +"-"+data.refblk;
                 }
                 ofile<<"# date jplFile Clock REF jmin jmax jR;";
                 ofile<<"# S_cut"<<std::endl;
                 ofile<<"# numClocks numEpochs"<<std::endl;
                 ofile<<"# r.n(rproj) rx ry rz sig svn-prn-clk-blk ...data..."
                      <<std::endl;
                 ofile<<"# ref-clock is last entry"<<std::endl;
                 // actual output:
                 ofile<<data.date<<" jpl"<<data.week<<data.day<<" "<<clockType
                      <<" "<<sREF
                      <<" "<<iBeg<<" "<<iEnd<<" "<<j<<std::endl;
                 ofile<<hh<<std::endl;
                 ofile<<numClocks<<" "<<2*jw+1<<std::endl;
                 //loop through all clocks (sat only)
                 for(int i=data.num_receivers;i<data.num_clocks;i++){
                   if(data.svn[i]==data.refsvn)continue; //skip ref!
                   std::string tempCLK=data.clk[i];
                   //int numclocks=0; //stores number of clocks
                   if((tempCLK==clockType)||(clockType=="all")){
                     ofile<<projpos[i]  //proj of pos onto "galactic vector"
                          <<" "<<data.pos[i][j][0]  //ECI x-position
                          <<" "<<data.pos[i][j][1]  //y ""
                          <<" "<<data.pos[i][j][2]  //z ""
                          <<" "<<data.sdev[i]           // the s.d. of this clock
                          <<" "<<data.svn[i]<<"-"<<data.prn[i]
                          <<"-"<<data.clk[i]<<"-"<<data.blk[i]<<" ";
                     for(int jj=iBeg;jj<=iEnd;jj++){
                       if((jj>=0)&&(jj<data.num_epochs)){
                         double tmpB=data.bias[i][jj];
                         ofile<<tmpB<<" ";
                       }else{
                         ofile<<0<<" ";
                       }
                     }
                     ofile<<std::endl;
                   }
                 }
                 //print out ref-clock details:
                 for(int i=data.num_receivers;i<data.num_clocks;i++){
                   if(data.svn[i]!=data.refsvn)continue; //only ref!
                   ofile<<projpos[i]  //proj of pos onto "galactic vector"
                        <<" "<<data.pos[i][j][0]  //ECI x-position
                        <<" "<<data.pos[i][j][1]  //y ""
                        <<" "<<data.pos[i][j][2]  //z ""
                        <<" "<<data.sdev[i]           // the s.d. of this clock
                        <<" "<<data.svn[i]<<"-"<<data.prn[i]
                        <<"-"<<data.clk[i]<<"-"<<data.blk[i]<<" ";
                   for(int jj=iBeg;jj<=iEnd;jj++){
                     ofile<<0<<" ";
                   }
                   ofile<<std::endl;
                 }//end output ref-clock
                 ofile.close();
               
                 break; //stop stepping window along
               }
             }//END step window across
           }//END loop over 'ref jump' epochs:
         }//END loop over +ve/-ve 'jumps'
       }//END loop over jump heights
       
     }//END loop over Cs (iC=0) and Rb (iC=1)
     
   }//END loop over days[j]
   MSC_progressBar("Searching jpl"+week,35);//reset counter 

   //print out time-remaining:
   if(!wpara){
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
   }

 }//END loop over weeks[i]


////////////////////////////////////////////////////////////////////////////////

 time (&end);
 double dif = difftime (end,start);
 printf ("Elasped time is %.0f seconds.\n", dif );
 return 0;
}
//******************************************************************************
//******************************************************************************





























