/*
PROGRAM:convertmout.cpp

-Converts JPL 30s clock data into a form that can be easily read by mathematica.
Two output formats:
  a) ".mout" format, suitable for input to Andrei's MMA notebooks
  b) A simpler ".mout2" format, that is a simple table


Format for ".mout2"/simple version:
1: date (human readable)
2: PRN
3: SVN
4: Block
5: Clock type
6+: Bias data
-Each column is a different clock, each row an epoch.
-Does not output the formal error, but can be 'normalised' by it.
-Good for just checking one clock. The students prefer it this way.

Uses openMP for parallelisation. Parallel by week.
The location of the input JPL files, and all other options are given inside
the input file: "convertmout.dat".
Detailed info on each of the options is also given inside the .dat file.




-------------------------------------------------------------------------
CHANGE LOG: Older versions archived in ./versions/ subdirectory
160401- New program. Converts JPL 30s clock files to ".mout" format for Mathem.
160404- Works.
160516- Now has option for either full "Andrei" format, or simple student format
160517- Now has option to signle-difference, double-difference, or neither.
160830- Changed sRefClk to array. NOT tested. Tested. Works.
161005- Now reads in ECI files, and can change the ref clock!
170830- Updating (NOT FINISHED) to work with JplGpsDataClass class.
        Done. (not tested)
-------------------------------------------------------------------------


*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//using namespace std;

#include "convertmout.h"


//******************************************************************************
//******************************************************************************
int main (void)
/*
Main Program: Loops through all input data in parallel. Calls other functions
*/
{

//Input Parameters:
 std::string path;           //="/home/ben/jpligsac/";
 std::string outbase;        //="/home/ben/jpligsac/";
 int minweek;           //=1060; first week with 30s data!
 int maxweek;           //=1879;
 int weekday;           //=7;   // 0-6 = Sun-Mon. 7=whole week.
 int iOverwrite;        //=1;   If file exists, skip, or overwrite?
 //int remunzip;
 //int zipout;
 int iPoly;             //=2;   //order of the fit poly
 int iOrder;            //=2;   //order of the fit poly
 int iWeight;           //=1;   //weighted fit(1)? or not(0)
 int iDD;               //Double difference, Signle Difference, or neither.
 int iFormat,iNorm;
 std::string label;          //Label for output files. "na" means no label!
 int icRef; //,irPos;
 std::string clockType,sRefBLK;

 std::string junk;
 std::ifstream fInput;
 fInput.open ("convertmout.dat");
   fInput >> path;                           getline(fInput,junk);
   fInput >> outbase;                        getline(fInput,junk);
   fInput >> minweek >> maxweek;             getline(fInput,junk);
   fInput >> weekday;                        getline(fInput,junk);
   fInput >> iOverwrite;                     getline(fInput,junk);
   fInput >> iPoly >> iOrder >> iWeight;     getline(fInput,junk);
   fInput >> icRef >> clockType >> sRefBLK;  getline(fInput,junk);
   fInput >> iDD;                            getline(fInput,junk);
   fInput >> iFormat;                        getline(fInput,junk);
   fInput >> iNorm;                          getline(fInput,junk);
   fInput >> label;                          getline(fInput,junk);
   //"getline" is a clumsy way of going to next line..
 fInput.close();

// Formats input, & reduces liklihood of errors
 int minday,maxday;
 if(weekday<0||weekday>7){weekday=7;}
 if(minweek<1060||minweek>2000){minweek=1060;}
 if(maxweek>2000){maxweek=2000;}
 if(maxweek<minweek)maxweek=minweek;
 if(weekday==7){
   minday=0;
   maxday=7;
 }else{
   minday=weekday;
   maxday=weekday+1;
 }
 if(iPoly<0||iPoly>1){iPoly=0;}
 if(iWeight<0||iWeight>1){iWeight=1;}
 if(iOrder<0){iOrder=2;}
 if(iDD<-1||iDD>2){iDD=0;}
 if(iOverwrite<0||iOverwrite>1){iOverwrite=0;}
 if(iFormat>1||iFormat<0){
   std::cout<<"iFormat="<<iFormat<<" Not valid."<<std::endl;
   return 1;
 }
 if(iNorm==3&&iFormat==1&&iDD<1){
   std::cout<<"Doesn't make sense to normalise by s.d., if you don't"
            <<" difference the data! Try again.\n";
   return 1;
 }

 {//Output params to screen:
   printf("\n~~~~~\n");
   printf("convertmout: Creating mout files for Mathematica.\n");
   if(iFormat==0){
     printf("output files will be '.mout' - Andrei's MMA format.\n");
   }else{
     printf("output files will be '.mout2' - simpler format.\n");
   }
   printf("Reading JPL 30s clock data, for weeks %i -> %i\n",minweek,maxweek);
   if(weekday!=7){printf(" for weekday %i\n",weekday);}
   printf(" Only for GPS/Sattelite clocks\n");
   if(iPoly==1){
     printf("Fitting and removing a %i-order polynomial",iOrder);
     if(iWeight==1){
       printf(", using a weighted fit.\n");
     }else{
       printf(" (not weighted).\n");
     }
   }else{
     printf("Not removing a polynomial.\n");
   }
   if(iDD==-1){
     printf("Not Differenceing the data.\n");
   }else if(iDD==0){
     printf("Not Differenceing the data, but centering it.\n");
   }else if(iDD==1){
     printf("Single Differenceing the data.\n");
   }else if(iDD==2){
     printf("Double Differenceing the data.\n");
   }
   printf("~~~~~\n\n");
 }



//// Read the PRN<->SVN map file
// string sPRNmap[PRNLINES][PRNPTS];
// readSVNPRN(path,sPRNmap,1);

//if omp is enabled, sets number of cores (only used for "% done")
 int numCores=1;
#if defined(_OPENMP)
 numCores=omp_get_num_procs();
#endif


 if(outbase!="0"){
   MSC_execute("mkdir -p "+outbase);
   if(!MSC_direxists(outbase)){
     std::cout<<"Directory "<<outbase<<" DNE. Create it, and try again.\n";
     return 1;
   }
 }

// MAIN LOOP:
 time_t start,mid,end;
 time (&start);
#pragma omp parallel for //runs through the weeks in parallel
 for(int i=minweek;i<maxweek+1;i++){  //loops through weeks
   for(int j=minday;j<maxday;j++){//loops through each day for given week

       JplGpsData data;
       int jplok=data.readJplBiasData(path,i,j);
       if(jplok==2)continue;//DIR doesn't exist
       if(jplok==1)continue;//this file doesn't exist

       //name for the output files:
       std::string outname;
       if(label=="na")label="";
       std::string filename="jpl"+data.week+data.day;
       if(outbase=="0"){
         outname=path+filename+"-"+label+".mout";
       }else{
         outname=outbase+filename+"-"+label+".mout";
       }

      //only do if outputfile doesn't exist, or if we should overwrite it
      if(MSC_fexists(outname)&&(iOverwrite==0))continue;

       std::cout<<"Reading: "<<filename<<" ("<<data.date<<". Ref: "
                <<data.refprn<<")\n";

       //De-trend (remove polynomial) from data
       if(iPoly==1){
         data.polynomialDetrend(iOrder,iWeight);
       }

       //Difference (or not) the data
       if(iDD!=-1) data.differenceData(iDD);

       bool use_svn=true;
       if(iDD>0){
         use_svn=false;
         data.calculateStdDev();
       }
       if(icRef==1){
         data.swapReference(clockType,sRefBLK,use_svn);
         std::cout<<" Ref --> "<<data.refname<<"\n";
       }

       writeMout(data,10,iFormat,iNorm);


   }//END for(int j=minday;j<maxday;j++)

   double perc=(double(i-minweek+1)/double(maxweek-minweek+1))*100.;
   if(perc<(100./numCores)){
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

 }//END loop over weeks

 time (&end);
 double dif = difftime (end,start);
 printf ("Elasped time is %.1f seconds.\n", dif );
 return 0;
}
//******************************************************************************
//******************************************************************************





//******************************************************************************
int writeMout(JplGpsData data, int prec, int iFormat, int iNorm)
/*
Writes the output data to file, Formatted for Andrei's MMA notebook 'mout'
Outputs a text file in specified directory called filename."mout"

INPUT:
  data      ::  clock data object [JplGpsData class]
  prec      ::  number of digits to write to output file
  iFormat   ::  which output format to use (see below)
  iNorm     ::  Normalise the data to sigma? See below

If iFormat=0: write with Andrei's format
if iFormat=1: write with simple "student" format
--
iNorm: "Normalises" the data (devides by 'sigma'). Only for "simple" format!
       0 => do nothing.
       1 => Normalise by formal error for each data point.
       2 => Normalise by average formal error for each clock
       3 => Normalise by standard deviation of each clock
       4 => Normalise by overall average of formal error

NB: only problem is "begin/end" time - Andrei had in in GPS time (epochs),
while this routine just does it day-by day.. [does output the data though]
--> no, actually fixed?

=== Change Log ===
160516- Now has option for either full "Andrei" format, or simple student format
160830- Updated sRefClk to array
170321- Removed "orb"
*/
{


 std::ofstream outFile;
 std::string filename="./mout/jpl"+data.week+data.day+".mout";
 if(iFormat==1)filename=filename+"2";//for "simple" format
 outFile.open (filename.c_str());

 //extract the year/month/day (integers)
 //...should just add this to the class
 int iY, iM, iD;
 std::stringstream ssin(data.date);
 ssin>>iY>>iM>>iD;
 iM=abs(iM); iD=abs(iD);//removes the "-", interpreted as a negative


 int week=std::stoi(data.week);
 int day=std::stoi(data.day);
 int iBeginTime=(week*7+day)*24*60*60;
 int iEndTime=iBeginTime+30*2879; //not including second midnight!

 std::string sREF;
 sREF=data.refname;

  if(iFormat==0){//Write with Andrei's format

   outFile<<"{0, "<<"\"Derived from jplwwwwd.clk_30s file. "
          <<sREF<<" clock fixed\", "
          <<"{"<<iY<<", "<<iM<<", "<<iD<<"}, \n"
          <<"{"<<iBeginTime<<"., "<<iEndTime<<"., "<<"30.}, \"GPS\", \n"
          <<"{";
   //loop over satellite clocks only:
   for(int i=data.num_receivers;i<data.num_clocks;i++){
     outFile<<data.svn[i];
     if(i!=data.num_clocks-1){
       outFile<<", ";
     }else{
       outFile<<"}, \n";
     }
   }
   outFile<<"{";
   for(int i=data.num_receivers;i<data.num_clocks;i++){
     outFile<<"{"
            <<data.svn[i]<<", "
            <<data.prn[i]<<", "
            <<"\""<<data.blk[i]<<"\", "
            <<"\""<<"ORB"<<"\", "
//            <<"\""<<sClock[i][ORB]<<"\", "
            <<"\""<<data.clk[i]<<"\""
            <<"}";
     if(i!=data.num_clocks-1){
       outFile<<", \n";
     }else{
       outFile<<"}, \n";
     }
   }
   outFile<<"{";
   for(int i=data.num_receivers;i<data.num_clocks;i++){
     outFile<<"{";
     for(int j=0;j<data.num_epochs;j++){
       outFile.precision(prec);
       outFile.setf(std::ios::fixed, std::ios::floatfield);
       outFile<<"{"
              <<data.bias[i][j]
              <<", "
              <<data.ferr[i][j]
              <<"}";
       if(j!=data.num_epochs-1){
         outFile<<", \n";
       }else{
         outFile<<"}";
       }
     }
     if(i!=data.num_clocks-1){
       outFile<<", \n";
     }else{
       outFile<<"}\n";
     }
   }
   outFile<<"}";

  }else{ //Write with simple format for students
  //This routine outputs the data in the form that was easy for Alex to use.

   //ref clock name:
   outFile<<data.date<<" Ref:"<<data.refname<<std::endl;

   //all other clock names (prn, svn, block, clock)
   for(int i=data.num_receivers;i<data.num_clocks;i++)
         outFile<<data.prn[i]<<" ";
   outFile<<std::endl;
   for(int i=data.num_receivers;i<data.num_clocks;i++)
         outFile<<data.svn[i]<<" ";
   outFile<<std::endl;
   for(int i=data.num_receivers;i<data.num_clocks;i++)
         outFile<<data.blk[i]<<" ";
   outFile<<std::endl;
   for(int i=data.num_receivers;i<data.num_clocks;i++)
         outFile<<data.clk[i]<<" ";
   outFile<<std::endl;

   //Calculate quantities for "normalising"
   std::vector<double>avfer(data.num_clocks);
   double oaavfer=0;
   if(iNorm==2||iNorm==4){//Calc average FER for each clock
     for(int i=data.num_receivers;i<data.num_clocks;i++){
       double tavg=0;
       for(int j=0;j<data.num_epochs;j++){
         tavg+=data.ferr[i][j];
       }
       avfer[i]=tavg/data.num_epochs;
       oaavfer+=avfer[i]/(data.num_satellites);
     }
   }

   //Output the clock data:
   for(int j=0;j<data.num_epochs;j++){
     for(int i=data.num_receivers;i<data.num_clocks;i++){
       outFile.precision(prec);
       outFile.setf(std::ios::fixed, std::ios::floatfield);
       if(iNorm==0)outFile<<data.bias[i][j]<<" ";
       if(iNorm==1){
         if(data.ferr[i][j]!=0){
           outFile<<data.bias[i][j]/data.ferr[i][j]<<" ";
         }else{
           outFile<<0.0<<" ";
         }
       }
       if(iNorm==2)outFile<<data.bias[i][j]/avfer[i]<<" ";
       if(iNorm==3)outFile<<data.bias[i][j]/data.sdev[i]<<" ";
       if(iNorm==4)outFile<<data.bias[i][j]/oaavfer<<" ";
     }
     outFile<<std::endl;
   }

  }

 outFile.close();

 return 0;
}
