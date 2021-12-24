/*
PROGRAM:

Plots the clock solutions for a given day.
Can plot all the clocks, or only a particular set, or a particular PRN/station.

Creates gnuplot files, and places them in a new directory called:
output-plotClockData/
[unless that directory doesn't exist and couldn't be created, in which case
it puts them in the current directory!]
-Can be deleted immediately if you no longer need them.

Either plots "to screen" (using gnuplot),
or to a pdf file (using gnuplot and pdflatex)
If you don't have gnuplot/latex installed, it will still work, it will just
create the .gnu/.tex files without actually plotting them.

Uses 'system' to call gnuplot/latex, so will only fully work on linux machines.
However, you can just call gnuplot/latex sepperately to plot to produced files.

Can plot raw data, polynomial-reduced data, or differenced data.
Plots either all clocks, or a specific clock (by PRN/CLK/Block)



=== Change Log ===
160512- Added wget lookup
160515- Added option to just print 1 specific clock
160517- Now has option to single-difference, double-difference, or neither
160811- Differencing with input of "0" no longer does "nothing", it removes the
         average. This program doesn't run the difference routine for iDD=-1:
         iDD=0 means don't difference, but do subtract mean.
         iDD=-1 means don't difference OR subtract mean.
160815- Looks up clock types
160830- Changed sRefClk to array.
160915- New "plot" routine. Much nicer.
170316- Checks if directory exists (i.e. non-unix friendly) + checkJPLfiles
170414- Minor changes.
170702- made static arrays (store on heap). no ulimit problem!
170827- Works with class
=== TO DO ===




*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#include "plotClockData.h"
#include <stdio.h>

//******************************************************************************
//******************************************************************************
int main (void){

  //input parameters:
  std::string path;//="/home/ben/jpligsac/";   //should be input
  int minweek;//=1060;//1060 - first week with 30s data!;
  int maxweek;//=1879;//1879;
  int weekday;//=7;
  int whichclocks;//=0;   //which clocks to examine: 0:Sat, 1:Station, 2:both
  int iPoly,iOrder;//=2;    //order of the fit poly
  int iWeight;//=1;   //weighted fit(1)? or not(0)
  int iPlot,iStyle;
  double dOffset;
  int iInit,iFin;
  double R2min;
  int iDD;               //Difference data: Signle, double, or neither
  std::string oneclock;  //Print 1 specific clock. Put "all" to print all

  //read in the input .dat file:
  std::ifstream fInput;
  fInput.open ("plotClockData.dat");
  {
    std::string junk;
    fInput >> path;                        getline(fInput,junk);
    fInput >> minweek >> maxweek ;         getline(fInput,junk);
    fInput >> weekday;                     getline(fInput,junk);
    fInput >> whichclocks;                 getline(fInput,junk);
    fInput >> iInit >> iFin;               getline(fInput,junk);
    fInput >> R2min;                       getline(fInput,junk);
    fInput >> iPoly >> iOrder >> iWeight;  getline(fInput,junk);
    fInput >> iDD;                         getline(fInput,junk);
    fInput >> iPlot >> iStyle >> dOffset;  getline(fInput,junk);
    fInput >> oneclock;                    getline(fInput,junk);
  }
  fInput.close();

  //Format the inputs (avoid some errors)
  int minday,maxday;
  if(whichclocks<0||whichclocks>2){whichclocks=0;}
  if(weekday<0||weekday>7){weekday=7;}
  if(maxweek<1060){maxweek=minweek;}
  if(iDD<0||iDD>2)iDD=-1;
  if(weekday==7){
    minday=0;
    maxday=7;
  }else{
    minday=weekday;
    maxday=weekday+1;
  }
  int max_epochs=2880; //just a "safe-guard" against bad parameters
  if(iWeight<0||iWeight>1){iWeight=1;}
  if(iPoly<0||iPoly>1){iPoly=1;}
  if(iOrder<0){iOrder=2;}
  if(iPlot<0||iPlot>2){iPlot=1;}
  if(iStyle<0||iStyle>2){iStyle=1;}
  if(iInit<0||iInit>max_epochs){iInit=0;}
  if(iFin<=0||iFin>max_epochs){iFin=max_epochs;}
  if(R2min<0){R2min=0;}
  if(oneclock=="0")oneclock="all";

  std::cout<<"\n          ###########################################";
  std::cout<<"\n                    Plotting 30s Clock Data          ";
  std::cout<<"\n          ###########################################\n\n";

  {//Output params to screen:
    printf("Plotting 30s clock data, for weeks %i -> %i\n",minweek,maxweek);
    if(weekday!=7)printf(" for weekday %i\n",weekday);
    printf(" %s clocks \n",oneclock.c_str());
    if(whichclocks==0){
      printf(" AS: Satellite clocks\n");
    }else if(whichclocks==1){
      printf(" AR: Base station clocks\n");
    }else if(whichclocks==0){
      printf(" Ax: Satellite AND base station clocks\n");
    }
    if(iPoly==1){
      printf("Fitting and removing a %i-order polynomial",iOrder);
      if(iWeight==1){
        printf(", using a weighted fit.\n");
      }else{
        printf(" (not weighted).\n");
      }
    }else{
      std::cout<<"Not detrending the data. \n";
    }
    if(iDD==0){
      printf("Not Differenceing the data.\n");
    }else if(iDD==1){
      printf("Single Differenceing the data.\n");
    }else if(iDD==2){
      printf("Double Differenceing the data.\n");
    }
    printf("~~~~~\n\n");
  }

  // *** Main loop over weeks ***
  for(int iw=minweek;iw<=maxweek;iw++){  //loops through weeks
    for(int id=minday;id<maxday;id++){//loops through each day for given week

      //create the data object:
      JplGpsData data;

      //set verbose=true so any messages will print:
      data.verbose=true;

      int jplok=data.readJplBiasData(path,iw,id);
      //int jplok=data.binaryJpl30s(path,iw,id);
      if(jplok==2)return 2;//DIR doesn't exist
      if(jplok==1)continue;//this file doesn't exist

      if(iPoly==1)data.polynomialDetrend(iOrder,iWeight);
      //Difference the data:
      if(iDD!=-1)data.differenceData(iDD);
      //Note, "zero differenceing" doesn't do nothing, it removes the average!

      //Work out which clocks to plot:
      int iCi=0; //index of intial clock to consider
      int iCf=0; // "" final ""
      if(whichclocks==0){
        iCi=data.num_receivers;
        iCf=data.num_clocks;
      }else if(whichclocks==1){
        iCi=0;
        iCf=data.num_receivers;
      }else if(whichclocks==2){
        iCi=0;
        iCf=data.num_clocks;
      }

      //make directory for gnuplot files, prepare the gnuplot files:
       //dir. to place gnuplot files:
      std::string outbase="./output-plotClockData/";
      MSC_execute("mkdir -p "+outbase);
      if(!MSC_direxists(outbase)){
        //For non-unix machines (new directory won't exist):
        std::cout<<"!! I couldn't make directory "<<outbase<<".\n";
        std::cout<<"!! Placing gnuplot files in current directory instead.\n";
        outbase="./";
      }

      std::string filename="jpl"+data.week+data.day;
      std::cout<<"Plotting file: "<<filename<<" ["<<data.refname<<"]\n";
      std::string title=filename+" ["+data.refname+"]. "+data.date+".";
      if(oneclock!="all")title+=" Only "+oneclock+" clocks.";

      plotClocks(data,oneclock,outbase,filename,R2min,iInit,iFin,
                 title,iCi,iCf,iPlot,iStyle,dOffset);
      //Use either gnuplot, or gnuplot+LaTeX to plot the files.
      //following only works on unix machines. Otherwise, just have to do
      //gnuplot/latex sepperately. Won't give error, will still output correct
      // .gnu files.

      if(iPlot==1){
        std::string cmd="gnuplot -persist "+outbase+filename+".gnu";
        std::cout<<cmd<<std::endl;
        MSC_execute(cmd.c_str());
      }else if(iPlot==2){
        std::string cmd="gnuplot "+outbase+filename+".gnu";
        std::cout<<cmd<<std::endl;
        MSC_execute(cmd.c_str());
        std::string cmd1="pdflatex "+outbase+filename+".tex";
        std::string cmd2="pdflatex -halt-on-error "+outbase+filename+".tex"+
                    " > /dev/null 2>&1";
        std::cout<<cmd1<<std::endl;
        MSC_execute(cmd2.c_str());
        MSC_execute("mv -f "+filename+".pdf "+outbase);
        MSC_execute("rm -f "+filename+".log");
        MSC_execute("rm -f "+filename+".aux");
        MSC_execute("rm -f "+outbase+filename+"-inc*");
      }

    }//END for(int j=minday;j<maxday;j++)
  }//END loop over weeks

  std::cout<<"Done :)"<<std::endl;

  return 0;
}
//******************************************************************************


//******************************************************************************
int plotClocks(JplGpsData data, std::string oneclock, std::string base,
                std::string filename, double R2min, int jmin, int jmax,
                std::string title,
                int iCi, int iCf, int iPlot, int iStyle, double dOfs)
/*
This functions uses gnuplot to plot the data. Actually, it just produces the
gnuplot file, which gnuplot can then use to make the plot.
Can either plot all clocks, or a particular group (by clock type and/or block)
or just a single clock (by PRN).

INPUT:
   data        :: Clock data object (JplGpsData)
   oneclock    :: Which clocks: e.g. "all", Cs", "IIR", "CsIIF", or "G12"
   base        :: Directory for output files
   filename    :: name for output file
   jmin, jmax  :: Initial/final epochs to plot. Default is 0->num_epochs
   title       :: A title for the plot. Can be blank.
   iCi, iCf    :: Initial/final clocks (usually to exclude the base stations)
   iPlot       :: How to format output file? 1=for screen, 2=pdf (latex)
   iStyle      :: Line style: 0=points, 1=lines, 2=lines+points
   dOfs        :: Amount (in ns) to offset each subsequent clock solution

=== Change Log ===
160915- New program, nicer than old one.
170321- Minor change to plot legends
170827- Works with class
===== TO DO =====

*/
{

  bool latex=false;
  if(iPlot==2){latex=true;}
  if(jmax==0)jmax=data.num_epochs;

  std::ifstream gnuHeader;
  std::ofstream gnuPlot;
  std::string gnuPlotfilename=base+filename+".gnu";
  std::string latexfilename=base+filename+".tex";
  gnuPlot.open (gnuPlotfilename.c_str());
  if(latex){
   gnuPlot<<"set terminal epslatex size 7.5,4.5 standalone color linewidth 3\n";
   gnuPlot<<"set output \'"+latexfilename+"\'\n";
  }else{
   gnuPlot<<"#set terminal epslatex size 7.5,4.5 standalone color linewidth 3\n";
   gnuPlot<<"#set output \'"+latexfilename+"\'\n";
  }


  if(title!=""){
    gnuPlot<<"set title \'"+title+"\'"+"\n";
  }

  //Decide where to place the "legend"
  if(oneclock=="all"||oneclock=="Rb"||oneclock=="Cs"||
     oneclock=="II"||oneclock=="IIA"||oneclock=="IIR"||oneclock=="IIF"||
     oneclock=="RbII"||oneclock=="RbIIA"||oneclock=="RbIIR"||oneclock=="RbIIF"||
     oneclock=="CsII"||oneclock=="CsIIA"||oneclock=="CsIIR"||oneclock=="CsIIF")
  {
     gnuPlot<<"set key outside\n";
  }

  //Writes the "plot xyz" command
    gnuPlot<<"plot ";
    bool first = true;
    for(int i=iCi;i<iCf;i++){
      if(data.blk[i]==oneclock||
         data.clk[i]==oneclock||
         data.prn[i]==oneclock||
         data.clk[i]+data.blk[i]==oneclock||
         oneclock=="all"){
        std::string sTitle=data.prn[i]+"-"+data.svn[i];
        if(oneclock=="all")sTitle+="-"+data.clk[i]+data.blk[i];
        if(data.clk[i]==oneclock)sTitle+="-"+data.blk[i];
        if(data.blk[i]==oneclock)sTitle+="-"+data.clk[i];
        if(data.prn[i]==oneclock)sTitle+="-"+data.clk[i]+data.blk[i];
        {
          if(data.dR2.size()>0){
            if(data.dR2[i]<R2min) continue;
          }
          if(!first){gnuPlot<<", ";}
          if(first) first=false;
          if(iStyle==0){
            gnuPlot<<"\"-\" title \'"<<sTitle<<"\' with errorbars ls "<<i;
          }else if(iStyle==1){
            gnuPlot<<"\"-\" title \'"<<sTitle<<"\' with lines ls "<<i;
          }else if(iStyle==2){
            gnuPlot<<"\"-\" title \'"<<sTitle<<"\' ls "<<i;
          }
          // if(i!=iCf-1){gnuPlot<<", ";}
        }
      }
    }
    gnuPlot<<"\n";

  double offset=0;
  for(int i=iCi;i<iCf;i++){
    if(data.blk[i]==oneclock||
       data.clk[i]==oneclock||
       data.prn[i]==oneclock||
       data.clk[i]+data.blk[i]==oneclock||
       oneclock=="all"){
      {
        if(data.dR2.size()>0){
          if(data.dR2[i]<R2min) continue;
        }
        for(int j=jmin;j<=jmax;j++){
        if((j>=0)&&(j<data.num_epochs))
          gnuPlot<<j<<" "<<data.bias[i][j]+offset<<" "<<data.ferr[i][j]<<"\n";
        }
        gnuPlot<<"e\n";
        offset=offset+dOfs;
      }
    }
  }
  gnuPlot.close();

  return 0;
}
