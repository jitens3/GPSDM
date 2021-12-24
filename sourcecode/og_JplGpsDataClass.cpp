//class JplGpsData::
#include "JplGpsDataClass.h"

/*
170826.
Contains member functions for the class JplGpsData

Each function returns 0 on success, 1 if there is a warning/error but the
program should continue, and 2 if a "fatal" error occurred, and the program
should stop.

Has member called "verbose", a bool, which is false by default.
Most 'cout's will not print unless verbose is true.
Some important error messages still will, and all the "fatal" error messages
still print.

At the moment, it is set up to read the PRN_GPS file every day (readSVNPRN).
This is not a big deal, but isn't necisary. Not sure the best way to deal with
this.

*/

//******************************************************************************
JplGpsData::JplGpsData(std::string path, int in_jplweek, int in_jplday,
                       bool detrend, bool difference)
/*
Constructor. All parameters optional.
the bools: detrend and difference are false(?) by default.
detrend=true    => weighted 2-order poly (the default)
difference=true => difference (1-order), AND calc s.d.
If difference=false, still 'zeros' the data (subtracts mean)
*/
{

  //set some default values:
  num_epochs=0, num_receivers=0, num_satellites=0, num_clocks=0;
  verbose=false;
  acf_points=5; // Here??

  //if no input given, don't do anything else
  if(in_jplweek==0) return; //(in_jplweek=0 by default)
    std::cout << "Here 1 JPL" << "\n";

  //otherwise: initialise arrays with JPL data:
  readJplBiasData(path, in_jplweek, in_jplday);
    std::cout << "Here 2 JPL" << "\n";
  readEciPos();
    std::cout << "Here 3 JPL" << "\n";
  if(detrend)polynomialDetrend();
  if(difference){
    differenceData(1);          // changed from 3, 1 is for network of only Rb clocks
    calculateStdDev();
    //swapReference("all");
  }else{
    differenceData(0);
  }

}

//******************************************************************************
template <class T>
int JplGpsData::reSizeVec(std::vector< std::vector<T> > &invec, int rows,
                          int cols)
/*
170823.
Private.
Small function to easily re-size the vector arrays.
*/
{
  invec.resize(rows, std::vector<T>(cols));
  return 0;
}

//----Overloaded version (1D vector):
template <class T>
int JplGpsData::reSizeVec(std::vector<T> &invec, int rows)
{
  invec.resize(rows);
  return 0;
}

//----Overloaded version (3D vector):
template <class T>
int JplGpsData::reSizeVec(std::vector< std::vector< std::vector<T> > > &invec,
                          int rows, int cols, int slcs)
{
  invec.resize(rows, std::vector< std::vector<T> >(cols, std::vector<T>(slcs)));
  return 0;
}



//******************************************************************************
int JplGpsData::readJplBiasData(std::string path, int in_jplweek, int in_jplday,
                                bool download)
/*
170823.
"New" function, based off readclk30s. Minor changes to work with new class.
Reads in the (already unzipped) clk_30s JPL clock files, and puts the data into
the relevant arrays.
Written for RINEX version 3.00.
Files from: ftp://sideshow.jpl.nasa.gov/pub/jpligsac/
Files located in directory "path/", and are called jpl+week+day+ext
It searches the files for the string "END OF HEADER" - before this is header,
after this is data, and then reads in the clock data, sorts it, and puts it in
the relevant arrays.
The input files are in seconds. This is converted into nanoseconds.
Every 10th epoch, the formal error (sigma) is wrong. This is fixed by
extrapolation.
Once the data has been read, this called the mapSVNPRN function, to map each
clock (here given by the PRN) to its SVN, clock type, block etc.

bool download is False by default. If true, it will download the JPL files
if they aren't already on hard-drive.

PRN/SVN index for satellite clocks comes from
ftp://sideshow.jpl.nasa.gov/pub/gipsy_products/gipsy_params/PRN_GPS.gz

INPUT
  path           ::  directory that holds JPL 30s data files
  in_jplweek     ::  jpl week (as an integer)
  in_jplday      ::  jpl day [0=sunday], integer

OUTPUT
  bias, ferr     ::  array for clock bias data (and formal error)
  prn, svn, clk,
  blk            ::  Array of prn/svn/clk/blk for each clock.
                    NB: blk, clk, and svn not filled until after map PRN_SVN
                    Same for ref clock refprn etc.
                    NB: PRN and SVN for ref clock are the 4char codes (eg, USN3)
                    NO! svn for all station clocks now 0!
  week           ::  string. jpl week
  day            ::  string. jpl day [0=sunday]
  date           ::  human readable date. e.g., 2016-02-26
  num_clocks,
  num_receivers,
  num_satellites ::  Number of: total clocks, Base stations(recievers),
                    satellites. num_receivers+num_satellites=num_clocks.


====== CHANGE LOG ======
160815- There was a bug that caused the first epoch of some clocks to be written
        as zero, even though it existed in the file. It occured when the
        previous clock had missing data points at the end of the day, which
        resulting in me "moving to the next line" too many times. Fixed with the
        use of the "nextclock" bool. Seems to work now.
        ...Now introduced another problem...incorrect num_receivers sometimes!
        -Fixed. [(k!=0) && (gpsEpoch==0)]
160822- No, last time, i didn't fix it. Now it is fixed I think:
        with: (gpsDay!=gpsDayline) && (gpsEpoch==0)
160830- Updated sRefClk to array
170321- Removed "Orb" reference; not relaible, never used.
170403- Had small error. If one of the 'every 10th' epochs was a skipped point,
        it would overwrite the 'zero' tag in FER. If one or both of the
        neighbouring epochs were non-zero, this would cause the 'missed' point
        to be treated as an actual point. Likely extremely rare case anyway.
170702- Changed fileLine to std::vector
170823- Updated to work with JplGpsData class.
170830- Updated to reflect JPL's file naming convention change: clk_30s -> clk
        Note: New ".clk" files do not include the 'second' midnight.
        This shouldn't matter. They also have WAY more recievers.
170912- Added "download" option, which is FALSE by default.
====== To Do ======

*/
{
    //

    
    std::cout << "Here 1 JPLreadata" << "\n";
  //sets the jplweek and day
  week=std::to_string(in_jplweek);
  day=std::to_string(in_jplday);
  //store the location of the data files:
  jpl_file_dir=path;
  std::cout << "this is jpl_file_dir: " << jpl_file_dir << "\n";
  //check if the CLK_30s file exists, and if not, attempt to download it
  int ijplok=checkJPLfiles(jpl_file_dir, download);
  std::cout << "this is ijplok: " << ijplok << "\n";
  if(ijplok!=0)return ijplok;

    std::cout << "Here 2 JPLreadata" << "\n";
  //Set the number of epochs: the JPL files have 2880 epochs
  num_epochs=2880;

  //Opens the input file.
  std::string pref="jpl";
  std::string ext=".clk_30s";
  if(in_jplweek>=1934){
//  if(in_jplweek>=1640){
    //After week 1934, JPL's naming convetnion changed
    ext=".clk";
  }
  std::string sJPLclk=path+pref+week+day+ext;
  std::ifstream inFile;
  inFile.open (sJPLclk.c_str());
    std::cout << "Here 3 JPLreadata" << "\n";

  //This counts the number of data lines.
  //It also find the "ref clock" name, and stores it in sRefName
  //Searches for  str:"END OF HEADER" - before this is header, after is data
  //#header lines (starting from 1, so actual number) stored in 'iNumJunkLines'
  //#data lines (from 0, for the array) stored in iTempLines->iLines
  int iTempLines=0, iNumJunkLines=0;
  int iLineNum=0;
  bool bData=false; //true once we get past the junk/header lines
  std::string sTempLine;
  std::string sRefName;
  while ( getline (inFile,sTempLine) ) {
    //Find the reference clock name:
    if (sTempLine.find("ANALYSIS CLK REF") != std::string::npos){
      std::stringstream srefin(sTempLine);
      std::string sRefClkCode;//Not used... but may want later (Ref Clock code)
      srefin>>sRefName>>sRefClkCode;
    }
    iLineNum++;
    if (!bData){
      //Find the start of the data lines (end of header/junk lines)
      if (sTempLine.find("END OF HEADER") != std::string::npos) {
        bData=true;
        iNumJunkLines=iLineNum;
      }
    }else{
      //count the # of data lines:
      iTempLines++;	//Just lines of clock DATA (not total lines)
    }
  }
//    std::cout << "Here 4 JPLreadata" << "\n";
  int iLines=iTempLines; //why a new variable? meh.

  //"rewind" the file, ready to be 'getlined' again:
  inFile.clear();
  inFile.seekg(0, inFile.beg);

  //writes Ref Clock info:
  refprn=sRefName;
  refsvn="00";
  refname=sRefName;
  refblk="AR"; // Means "reciever", i.e. eaerth-based clock
  refclk="H";  // The reference clocks for JPL files always H-masers (??)

  //Loads the each line of the input data file into string array
  //NOTE: we could actually do this in the above loop, using pushback()!
  // Doesn't matter. This routine is fast enough, and called seldomly.
  std::vector<std::string> fileLine(iLines);
  std::string sLine;
  int iLine=0;
  iLineNum=0;
  while ( getline (inFile,sLine) ) {
    iLineNum++;
    if((iLineNum>iNumJunkLines)){
      fileLine[iLine]=sLine;
      iLine++;
    }
  }
 //   std::cout << "Here 5 JPLreadata" << "\n";

  inFile.close(); //Don't need the file any more

  //Sorts the array (of lines) by clock (PRN) and time
  std::sort(fileLine.begin(), fileLine.end());

 //   std::cout << "Here 6 JPLreadata" << "\n";
  //Counts the total number of clocks:
  {
    std::string sCurrentClock,sCurrentType,sPrevClock,sPrevType;
    std::stringstream ssin(fileLine[0]);
    ssin>>sPrevType>>sPrevClock;
    int iTClocks=1,iTAR=0,iTAS=0;
    for(int i=1;i<iLines;i++){
      std::stringstream ssin(fileLine[i]);
      ssin>>sCurrentType>>sCurrentClock;
      if(sCurrentClock!=sPrevClock){
        iTClocks++;
      }
      if(sCurrentType!=sPrevType&&iTAR==0){
        iTAR=iTClocks-1;
      }
      sPrevClock=sCurrentClock;
    }
    iTAS=iTClocks-iTAR;
    num_clocks=iTClocks;  //Number of clocks (sats + stations)
    num_receivers=iTAR;   //Number of stations/recievers (earth-based clocks)
    num_satellites=iTAS;  //Number of satellite clocks
  }
 //   std::cout << "Here 7 JPLreadata" << "\n";

  num_clocks = num_satellites; //We don't care about station clocks!
  num_receivers = 0;

  //Now we know the number of clocks, we can "re-size" the vector arrays:
  reSizeVec(bias,num_clocks,num_epochs);  //clock bias (2d vector array)
  reSizeVec(ferr,num_clocks,num_epochs);  //clock formal error
  reSizeVec(prn, num_clocks); //prn ("G01" etc)
  reSizeVec(svn, num_clocks); //SVN (space vehicle number)
  reSizeVec(clk, num_clocks); //Clock type (Cs, Rb, H)
  reSizeVec(blk, num_clocks); //Block: II, IIA, IIR, IIF, or AR (for recievers)
  //  std::cout << "Here 8 JPLreadata" << "\n";

  // XXX reSize here?? Or when calc'd!!
  //reSizeVec(dR2,num_clocks);//r^2 value (not used yet, but don't want seg-flts)
  //reSizeVec(sdev,num_clocks); //Stores the std deviation for each clock.

  //Puts clock data into the data (bias/ferr) arrays
  //The complicating factor here, is that not all epochs are always recorded!
  // i.e. some clocks have less than 2880 (2881) data points!
  std::string Ax, sPrn;
  int yy,mm,dd,hh,min,idl,gpsEpoch,gpsDay=-1,gpsDayline;
  double sec,dbias,sigma;
  int iYear=0,iMonth=0,iDay=0; //For final date
  double nano=pow(10.,9.);      //sec to nano-seconds
  int i=0,j=0;		//i: line number; j: clock #
  bool nextclock=false;
  std::string temp = "";
  while(i<iLines){
    for(int k=0;k<num_epochs+1;k++){//includes the unwanted midnight
      Ax="0";sPrn="0";dbias=0;sigma=0;
      if(i<iLines){
        std::stringstream ssin(fileLine[i]);
        ssin>>Ax>>sPrn>>yy>>mm>>dd>>hh>>min>>sec>>idl>>dbias>>sigma;
      }

      gpsEpoch=(hh*60*60+min*60+int(sec))/30;
      if(i==0)gpsDay=dd; //initial day
      gpsDayline=dd; //day for this line (only different for 'second' midnight)
      dbias=dbias*nano; //from sec -> ns
      sigma=sigma*nano; //from sec -> ns

      if(Ax == "AR"){
        k= num_epochs+1;
        i++;
        temp = Ax;
        continue;
      }
    //    std::cout << "Here 8.5 JPLreadata" << "\n";

      temp = Ax;

      if(k==0){
        //If first line of new clock, write clock info to clock array
        //SVN and block not known untill read-in PRN_GPS file.
        prn[j]=sPrn;
        svn[j]="00";
        clk[j]="**";
        blk[j]=Ax;   //For station clocks, this field will be "AR". XXX type?
        iYear=yy;
        iMonth=mm;
        iDay=dd;
      }
      if( (gpsDay!=gpsDayline) && (gpsEpoch==0) ){//ignore the second midnight!
        if(!nextclock)i++;
        //Ignore this line, go to next line (and clock!)
        //but not if the data file was missing the last
        //points, because in that case, we're already on the
        //next line!
      }else if(gpsEpoch==k){
        nextclock=false;
        bias[j][k]=dbias;
        ferr[j][k]=sigma;
        i++; //go to next line in file
      }else if (k<gpsEpoch){
        bias[j][k]=0;	//This epoch missing from file; write zeros
        ferr[j][k]=0;	//They should be excluded from analysis
        //nextclock=true;
      }else if (k<num_epochs){
        //Data stopped pre-midnight. NO i++. Fill w/ zeros
        //(note: k>gpsEpoch here!)
        bias[j][k]=0;
        ferr[j][k]=0;
        nextclock=true;//see above { if(!nextclock) }
      }//ENDIF
    }//END for(int k=0;k<num_epochs;k++)
    if(temp == "AS"){
      j++;
    }// goes to next clock
  }//END while(i<iLines)
  //  std::cout << "Here 9 JPLreadata" << "\n";

  //Get the date (in human-readable format).
  date=MSC_padIntString(iYear)+"-"+MSC_padIntString(iMonth)+"-"
       +MSC_padIntString(iDay);

  //once data is read in, use the PRN_GPS file to map the SVNs etc.
  mapSVNPRN(jpl_file_dir);
  // and the station clocks:
  mapStations(jpl_file_dir);

  //The formal errors given by JPL are incorrect every 10th epoch.
  //This is a known issue that will likely be fixed in the future.
  //The are understeimated by ~1000x.
  //We use the formal errors for 'weights' when removeing a polynomial, so
  // without fixing this, the weighting would essentially mean only these (every
  // 10th epoch) would be used in the fitting.
  //This section 'fixes' the incorrect errors by extrapolation (+/-1 epoch)
  //The formal errors are not extremely important, and they typically vary
  // slowly and only by a small fraction, so this is not a big deal.
  // After (and including) week 1934, JPL fixed this!
  if(in_jplweek>=1934) return 0; //don't need to fix!
  for(int i=0;i<num_clocks;i++){
    for(int j=0;j<num_epochs;j++){
      if(j%10==0){//This extrapolates the errors for the "wrong" points
        if(j==0){
          double newsig=ferr[i][j+1];
          //don't overwrite if this is a 'missed' epoch!
          if(ferr[i][j]!=0)ferr[i][j]=newsig;
        }else if(j<num_epochs-1){
          double newsig=(ferr[i][j-1]+ferr[i][j+1])/2.;
          if((ferr[i][j-1]==0)||(ferr[i][j+1]==0))newsig*=2;
          //(If one of these is zero, will devide by two when not averaging!)
          //if they are both zero, will write a zero, effectively deleting this
          //epoch.
          if(ferr[i][j]!=0)ferr[i][j]=newsig;
          if(newsig==0)bias[i][j]=0; //kill this epoch.
        }
      }
    }
  }
    
    num_satellites = 0;         // made this on 09/30/20 tday
    
    for(int i=0;i<num_clocks;i++){
        
        std::cout << "this is prn in array order: " << prn[i] << "\n";
        
    }
    
    // for some reason not finding svn 40?
    for(int i=0;i<num_clocks;i++){
        
        std::cout << "this is svn in array order: " << svn[i] << "\n";
        std::cout << "this is clk type in array order: " << clk[i] << "\n";
        
    }
  
    // added 09/30/20 for Cs clock elimination
    int cs_count;
    cs_count = 0;
    //std::vector<int> tmpsvn(num_clocks);
    reSizeVec(tmpsvn, num_clocks);
    reSizeVec(tmpsvnord, num_clocks);
    for(int i=0;i<num_clocks;i++){
        tmpsvn[i] = 99;
    }
    for(i=0;i<num_clocks;i++){
       // std::cout << "Hello: " << i << "\n";
        if(clk[i]=="Cs"){
            tmpsvn[i] = i;
            tmpsvnord[i] = svn[i];
            cs_count++;
            std::cout << "this is array order of Cs clock: " << tmpsvn[i] << "\n";
            std::cout << "this is array sorted order of Cs clock: " << tmpsvnord[i] << "\n";
        }
        //std::cout << "this is array order of Cs clock: " << tmpsvn[i] << "\n";
    }
    std::cout << "this is number of Cs clocks: " << cs_count << "\n";
    std::cout << "this is number of all clocks: " << num_clocks << "\n";
    std::cout << "this is number of satellites: " << num_satellites << "\n";
    num_satellites = cs_count;      // num_satellites is now number of Cs clocks
    std::cout << "this is number of satellites: " << num_satellites << "\n";
    
    for(int i=0;i<num_clocks;i++){
        std::cout << "this is tmpsvn[i]: " << tmpsvn[i] << "\n";
    }
    
    for(int i=0;i<num_clocks;i++){
        std::cout << "this is tmpsvnord[i]: " << tmpsvnord[i] << "\n";
    }
    
    // end of 09/30/20 additions
    
    // added on 09/02/20 for eliminating a clock from network
    
    std::vector< std::vector<double> > tbis(num_clocks,std::vector<double>(num_epochs));
    std::vector< std::vector<double> > tfer(num_clocks,std::vector<double>(num_epochs));
    
    for(int i=0;i<num_clocks;i++){
        
        for(int j=0;j<num_epochs;j++){
            
            tbis[i][j] = bias[i][j];
            tfer[i][j] = ferr[i][j];
            
        }
        
    }
    
    // end of new on 09/02/20
    
    // 08/17/20 new addtions for the skipping of one clock G#
 //   num_clocks = num_satellites - 1; //We don't care about station clocks!
    num_receivers = 0;
    
    // 09/16/20 new addtions for the skipping of all Cs clocks
    num_clocks = num_clocks - num_satellites;
    std::cout << "this is number of all clocks after elim: " << num_clocks << "\n";
    
    std::string gvalPath = "./results/gval13101.out";
    std::ofstream g2val(gvalPath);

    //Now we know the number of clocks, we can "re-size" the vector arrays:
    reSizeVec(bias,num_clocks,num_epochs);  //clock bias (2d vector array)
    reSizeVec(ferr,num_clocks,num_epochs);  //clock formal error
    reSizeVec(prn, num_clocks); //prn ("G01" etc)
    reSizeVec(svn, num_clocks); //SVN (space vehicle number)
    reSizeVec(clk, num_clocks); //Clock type (Cs, Rb, H)
    reSizeVec(blk, num_clocks); //Block: II, IIA, IIR, IIF, or AR (for recievers)
    //std::cout << "Here 8 JPLreadata" << "\n";

    int nit;
    nit = 0;
    
    
    // this is for a one clock elimination (phase slip filter)
//    for(int i=0;i<num_clocks+1;i++){
     
    // this is for elimination of entire Cs clocks from network
    //for(int i=0;i<num_clocks+9;i++){
    for(int i=0;i<num_clocks+num_satellites;i++){
     
     // this if statement is only for getting rid of a single clock (phase slip)
  //      if(i==24){
  //          continue;
  //      }
     // end of phase slip if statement
     
     // this if statement is to get rid of all Cs clocks from network
//      if(i==0 || i==2 || i==4 || i==7 || i==8 || i==9 || i==13 || i==22 || i==23){
//          continue;
//      }
        
      if(i==tmpsvn[i]){
          std::cout << "Hello: " << i << "\n";
          continue;
      }
      // this if statement does the same thing as the if statement on lines 512-514
     // end of elimination of Cs clocks from network 09/16/20
     
      for(int j=0;j<num_epochs;j++){
          bias[nit][j] = tbis[i][j];
          g2val <<  bias[nit][j] << "\n";
          ferr[nit][j] = tfer[i][j];
        }
          nit = nit + 1;
      }
        
    //}

    // this is for eliminating a single clock (phase slip filter)
//    for(int i=0;i<num_clocks+1;i++){
    
    // this is for elimination of the all Cs clocks from network
    for(int i=0;i<num_clocks+num_satellites;i++){
        for(int j=0;j<num_epochs;j++){
            
            if(j==0){
             std::cout << bias[i][j] << "\n";
            }
            
        }
        
    }
    
    g2val.close();
    
    // end of new additions on 08/17/20
    
    
  return 0;
}


//******************************************************************************
int JplGpsData::checkJPLfiles(std::string path, bool download)
/*
170316.
Checks to see if the JPL files exist. If not, uses wget to fetch them from
ftp://sideshow.jpl.nasa.gov/pub/jpligsac/
Will only work on linux machines.
====== CHANGE LOG ======
170830- Updated to work with JplGpsData class.
      - Updated to reflect JPL's file naming convention change: clk_30s -> clk
170912-
*/
{
  std::string ext=".clk_30s";
  std::string zip=".Z";

  if(std::stoi(week)>=1934){
//  if(std::stoi(week)>=1640){
    //After week 1934, JPL's naming convetnion changed
    ext=".clk";
    zip=".gz";
  }

  std::string filename="jpl"+week+day;

  if((!MSC_fexists(path+filename+ext+zip))&&(!MSC_fexists(path+filename+ext))){
    //If file not in directory, try to wget it from web
    std::string fileURL="ftp://sideshow.jpl.nasa.gov/pub/jpligsac/"+week+"/"
                   +filename+ext+zip;
    MSC_execute("mkdir -p "+path);
    if(!MSC_direxists(path)){
      //The "create directory (mkdir)" command above only works on linux pcs:
      std::cout<<"ERROR 116 in JplGpsData::checkJPLfiles: Directory: "<<path
               <<" doesn't exist. You have to create it manually.\n";
      return 2;
    }
    if(download){
      if(verbose)std::cout<<"\nAttempting to find "<<filename
                          <<" on jpl.nasa.gov: ..."<<std::flush;
      MSC_execute("wget -q -P "+path+" "+fileURL);
//      //uses wget to download file.
//      //Downloads first to sub-folder (wgettemp), then copies back to main dir
//      //Does this to be "safe", so doesn't ever 'half' download file ???
//      MSC_execute("mkdir -p "+path+"/wgettemp/");
//      MSC_execute("wget -q -P "+path+"/wgettemp/"+" "+fileURL);
//      MSC_execute("mv -f "+path+"/wgettemp/"+filename+ext+zip+" "+path
//                  +" 2>/dev/null"); //2>/dev/null suppresses output msgs
// Nope: doesn't work.
    }
    if(MSC_fexists(path+filename+ext+zip)){
      if(verbose)std::cout<<"found! Continuing.\n";
    }else{
      if(verbose)std::cout<<"\n"<<filename<<" NOT found! Skipping.\n"
                          <<std::flush;
      return 1;
    }
  }
  //unzips
  if(!MSC_fexists(path+filename+ext)){
    MSC_execute("gunzip -f -q "+path+filename+ext+zip);
  }
  // removes left-over unziped
  if(MSC_fexists(path+filename+ext+zip)){
    MSC_execute("rm -f "+path+filename+ext+zip);
  }

  return 0;
}



//******************************************************************************
int JplGpsData::mapSVNPRN(std::string DIR, bool gps_dm)
/*
2016-04-01
Maps the PRNs (from the input JPL 30s clock files) to the SVN/clock/block info.
==
The sPRNmap array is created by readSVNPRN function.
Fills info for the SVN, Block, and clock type.
-Should work the same with JPL/APD version of GPS_PRN! (yep, checked!)

INPUTS:
  DIR       ::  location of the PRN_GPS file
  gps_dm    ::  true (dafult): use "our" updated version. Fale=use orig. ver

====== CHANGE LOG ======
160720- Now writes out the failed PRN-GPS lookups to a txt file
170315- Also writes out if SVN is found, but CLK not known.
170321- Removed "Orb" reference; not relaible, never used.
170824- Works with class.
170905- An error occured with std::stoi when there was a missing PRN-SVN
        assignment. i.e., we went passed the end of the PRN_GPS file without
        finding the correct mapping, and the argument of stoi() was not a string
        that could be converted to an integer. Fixed this by changing the length
        of sPRNmap such that it is exactly as long as it needs to be.
        Also added if(tprn==" ")break; that should also fix it.
170925- Updated. Uses better date handling.
*/
{
  if(prn.size()==0){
    std::cout<<"ERROR 442 in JplGpsData::mapSVNPRN: no JPL data?\n";
    return 2;
  }

  //Array to hold the PRN<->SVN mappings.
  //XXX Make this a private array! DON"T need to read the file each time!!
  // This vector is re-sized inside 'readSVNPRN'
  std::vector< std::vector<std::string> > sPRNmap;

  //call the function that reads in the PRN_GPS file:
  readSVNPRN(DIR,sPRNmap,gps_dm);

  //work out the date for this JPL file (in integer wwwwd format):
  int today = std::stoi(week)*10 + std::stoi(day);
    
    //std::cout << "num_clocks in mapSVN: " << num_clocks << "\n";
    //std:: cout << "this is prn[10]: " << prn[9] << "\n";

  //Loop through the clock array and the PRN<->SVN mapping, do the mapping
  // It tries to find the entry in the PRN_GPS file where the current date
  // (from the JPL 30s file) is inbetween the start/end dates in PRN_GPS.
  // This is probably possible with nicer code using ctime features.
  for(int i=num_receivers;i<num_clocks;i++){//only do sattelites.
  
    blk[i]="**"; //default value if not found (svn, clk already set)

    for(std::size_t j=0; j<sPRNmap.size(); j++){
      //Goes through each line of sPRNmap, and checks to see
      //if it corrsponds to the clock we want.
      std::string tprn = sPRNmap[j][3];
        if(tprn == "G010"){
            tprn = "G10";
        }
        //std::cout << "tprn: " << tprn << "\n";
      if(tprn!=prn[i]) continue; //skip non-correct lines
        //std::cout << "tprn: " << tprn << "\n";
      std::string ini_date = sPRNmap[j][0];
      std::string fin_date = sPRNmap[j][1];
      std::string tsvn = sPRNmap[j][2];
      std::string tblk = sPRNmap[j][4];
      std::string tclk = sPRNmap[j][5];
      int beg, end;
      //convert the dates from the map array to jpl wwwwd (int) format:
      beg = MSC_ymdtowd(ini_date);
      if(fin_date=="0000")
        end = 999999; //just really big number. Assignment still current.
      else
        end = MSC_ymdtowd(fin_date);
      //check if todays date falls in range of this line:
      if(today>=beg && today <=end){
        //Found the right entry!
        svn[i]=tsvn;
        blk[i]=tblk;
        clk[i]=tclk;
        break;
      }
    }//goes through each line until found. (or eof)

    //check if not found:
    if((clk[i]!="Rb")&&(clk[i]!="Cs")){
      //Will write which PRN lookups failed to a log file.
      //found SVN, but don't know which clock used!
      //ALSO write this out
      std::ofstream oFile;
      std::string filename="PRN_GPS-faillog.txt";
      oFile.open (filename.c_str(),std::ios_base::app);
      oFile<<date<<" PRN:"<<prn[i]<<" SVN:"<<svn[i]<<" "
           <<clk[i]<<"\n";
      oFile.close();
    }

      
      //std::cout << "this is svn in order: " << svn[i] << "\n";        // length 30 (total) for now
  }//End loop over clocks

  return 0;
}

//******************************************************************************
int JplGpsData::readSVNPRN(std::string DIR,
        std::vector< std::vector<std::string> > &sPRNmap, bool gps_dm)
/*
Private function.
2016-04-01.
Reads in the "PRN_GPS" file that maps the PRN to the SVN,
sat Block, clock type etc., and puts into an array ordered by line number!
PRN/SVN index for satellite clocks comes from
ftp://sideshow.jpl.nasa.gov/pub/gipsy_products/gipsy_params/PRN_GPS.gz
If the file is not on the hd, uses wget to get it.
This means, if you want to use updated files, just delete the existing one!

"DIR" is the base directory in which the file goes.
outputs an array sPRNmap[PRNLINES][PRNPTS], that includes the dates and all info

Can either use the JPL version of the file (but note this has some incorrect
clock assignments), or the corrected "GPSDM" version:
gps_dm=false => Use JPL version
gps_dm=true => Use our updated version

INPUTS:
   DIR      ::  Base directory where the PRN_GPS file lives
   gps_dm   ::  Which PRN_GPS file? false=> JPL, true=> GPSDM (fixed) version
OUTPUTS:
   sPRNmap  ::  Array that holds the mapping info PRNs->SVNs etc.

---
Change Log:
160627- Now works for either Andrei's version, or JPL version!
170321- Removed "Orb" reference; not relaible, never used.
170824- Works with class
170925- Updated. Uses better date handling.
*/
{

  std::string fileURL=
    "ftp://sideshow.jpl.nasa.gov/pub/gipsy_products/gipsy_params/PRN_GPS.gz";
  std::string gpsdm_ver="_GPSDM.txt"; //For our updated version
  std::string sPRNGPS=DIR+"PRN_GPS";

  if(gps_dm){
    sPRNGPS=sPRNGPS+gpsdm_ver;
  }else{
    if((!MSC_fexists(sPRNGPS))&&(MSC_fexists(sPRNGPS+".gz"))){
      //if the zipped, but not unzipped file exists, unzip it.
      MSC_execute("gunzip -f -q "+sPRNGPS+".gz");//unzips
    }
    if(!MSC_fexists(sPRNGPS)){
      //if file doesn't exist, use wget to get it, and unzips it
      MSC_execute("wget -q -P "+DIR+" "+fileURL);
      MSC_execute("gunzip -f -q "+sPRNGPS+".gz");
    }
  }

  if(!MSC_fexists(sPRNGPS)){
    //failure if doesn't exist
    if(verbose){
      std::cout<<"FAILURE 629 in readSVNPRN: PRN lookup file "<<sPRNGPS
               <<" does not exist!\n";
    }
    return 1;
  }

  //Open the PRN_GPS file
  std::ifstream inFile;
  inFile.open (sPRNGPS.c_str());

  //Read the file
  std::string sLine;
  bool firstline=true;
  while(getline(inFile,sLine)){
    if(firstline){
      firstline=false;//skips the first junk line
      continue; //skip the first junk line
    }
    std::stringstream ssTest(sLine);
    std::string sTest;
    ssTest>>sTest;
    if(sTest=="!")continue;//skips "!" comment lines (GPSDM version)
    std::stringstream ssin(sLine);
    std::string ini_date, fin_date, l_blk, l_orb, l_clk;
    int l_svn, l_prn; //ints so can "pad" easily.
    ssin >> ini_date >> fin_date >> l_svn >> l_prn >> l_blk >> l_orb >> l_clk;
    std::vector<std::string> temp_vec(6);
    temp_vec[0]=ini_date;
    temp_vec[1]=fin_date;
    temp_vec[2]=MSC_padIntString(l_svn);
    temp_vec[3]="G"+MSC_padIntString(l_prn);
    temp_vec[4]=l_blk;
    temp_vec[5]=l_clk;
    sPRNmap.push_back(temp_vec);
  }//END WHILE (getline(inFile,sLine))
  inFile.close();

  return 0;
}


//******************************************************************************
int JplGpsData::mapStations(std::string path)
/*
170925.
Uses the station_map array to work out which station employes which clock.

Format of station_map:
  Name Clock intial_date final_date
Final date of 'c' mean "current" assignment.
Final date of '?' means not sure; if no future line, assume current assignment.

NOTE:
The start/end dates overlap! I include the start date, but not the end date in
the mappings. I don't know if this is correct/ok. (Probably, clocks not used for
at least a day or so when they are swapped??)
====== CHANGE LOG ======

*/
{
  if(prn.size()==0){
    std::cout<<"ERROR 761 in JplGpsData::mapStations: no JPL data?\n";
    return 2;
  }

  //Array to hold the PRN<->SVN mappings.
  //XXX Make this a private array! DON"T need to read the file each time!!
  // This vector is re-sized inside 'readStaMap'
  std::vector< std::vector<std::string> > station_map;

  //call the function that reads in the PRN_GPS file:
  readStaMap(path,station_map);

  //work out the current year/month/day for this JPL file:
  // in the integer wwwwd format!
  int today = std::stoi(week)*10 + std::stoi(day);

  // Loop through each station clock, work out what clock is used.
  for(int i=0;i<num_receivers;i++){//only do recievers/stations
    std::string name = prn[i]; //name of this clock
    for(std::size_t l=0; l<station_map.size(); l++){
      std::string sta_name, sta_clock, ini_date, fin_date;
      sta_name = station_map[l][0]; //data from the map
      sta_clock =station_map[l][1];
      ini_date = station_map[l][2];
      fin_date = station_map[l][3];
      if(name!=sta_name) continue; //skip wring files
      //intial/final dates for this assignment:
      int beg, end;
      //convert the dates from the map array to jpl wwwwd (int) format:
      beg = MSC_ymdtowd(ini_date);
      if(fin_date=="c"){
        end = 999999; //just really big number. Assignment still current.
      }else if(fin_date=="?"){
        //Ambigues notation. Check if this is the last entry for this station.
        // if so, assume assignment current. Else, go to next entry.
        if(l+1<station_map.size()){
          if(station_map[l+1][0]==name) continue; //go to next line
        }
        //otherwise, assume '?' means 'c'
        end = 999999; //just really big number. Assignment still current.
      }else{
        end = MSC_ymdtowd(fin_date);
      }

      if(today>=beg && today<end){
        //Found correct assignment!
        //Store info
        svn[i]="00"; //All station clocks have 0?
        clk[i]=sta_clock;
        //blk and prn already stored.
        break;
      }

    }//goes through each line until found. (or eof)

  }//End loop over clocks

   return 0;
}


//******************************************************************************
int JplGpsData::readStaMap(std::string path,
        std::vector< std::vector<std::string> > &station_map)
/*
170925.
Reads in the station clock map file (created by me from the IGS station logs).
Forms the vector map: station_map, which is used to map each station to which
clock was used.

Format of station_map:
  Name Clock intial_date final_date
Final date of 'c' mean "current" assignment.
Final date of '?' means not sure; if no future line, assume current assignment.

*/
{
  std::string map_file_name = path+"/STA_CLOCK_GPSDM.txt";

  if(!MSC_fexists(map_file_name)){
    //failure if doesn't exist
    if(verbose){
      std::cout<<"FAILURE 765 in readStaMap: Station lookup file "
               <<map_file_name<<" does not exist!\n";
    }
    return 1;
  }

  //Open the Station map file
  std::ifstream in_file;
  in_file.open (map_file_name.c_str());

  //number of entries per line:
  // Name Clock intial_date final_date
  int number_entries=4;

  //Read the file
  std::string sline;
  bool firstline=true;
  while(getline(in_file,sline)){
    if(firstline){
      firstline=false;//skips the first junk line
      continue; //skip the first junk line
    }
    //Check for "comment" lines (marked with "!")
    std::stringstream sstest(sline);
    std::string stest;
    sstest >> stest;
    if(stest=="!")continue;//skips "!" comment lines (GPSDM version)
    //Read in this line of the station map file
    std::stringstream ssin(sline);
    std::string sta_name, sta_type, sta_clock, ini_date, fin_date;
    ssin >> sta_name >> sta_type >> sta_clock >> ini_date >> fin_date;
    std::vector<std::string> map_line(number_entries);
    map_line[0] = sta_name;
    //Parse/format the clock type:
    if(sta_type=="EXTERNAL"){
      if(sta_clock=="H-MASER"||sta_clock=="H_MASER"||sta_clock=="MASER"
         ||sta_clock=="H-maser"||sta_clock=="18"||sta_clock=="39"
         ||sta_clock=="1"||sta_clock=="13")
        map_line[1]="H";
      else if(sta_clock=="CESIUM") map_line[1]="Cs";
      else if(sta_clock=="RUBIDIUM") map_line[1]="Rb";
      else map_line[1]="??";
    }else{
      //not external clock => internal (?)
      map_line[1] = "int";
    }
    map_line[2]=ini_date;
    map_line[3]=fin_date;
    station_map.push_back(map_line);
  }

  return 0;
}

//******************************************************************************
int JplGpsData::readEciPos(std::string which)
/*
170823.
Reads the ECI position files for the satellites and/or base-stations (receivers)
and stores them in the arrays pos and refpos.
Satellite and station positions files have name formats ecisatWWWWD.pos and
ecistaWWWWD.pos, respectively.
Note: for now, it doesn't store the velocities as these are not yet used;
this can easily be updated.
Written for ECI files generated by Geoff using GYPSY software.
This function can only be called {\em after} the JPL files are read in
(so it knows which satellite positions to put into which array positions).

Public function that calls one of two (or both) private functions that
read in the ECI position files for the satelites and reciever stations.
Note: 'which' is an optional parameter, which defaults to 'both'

Note: Must have "sta" or "both" in order to read the position of the reference
clock (otherwise, they will be zero).
which="both" is dafault.

NOTE: The mapSVNPRN function MUST be called prior to this one, since this
function uses the SVNs to map which orbit belongs to which clock.
XXX - add bool!

INPUTS:
   which    :: Which files to read in: sat or sta or both (default)
*/
{
    std::cout << "Here 1 ECIpos" << "\n";
    std::cout << "num_clocks: " << num_clocks << "\n";
  //if(which!="sat"&&which!="sta") return 2;
  if(!(which=="sat"||which=="sta"||which=="both")) return 2;
    std::cout << "Here 2 ECIpos" << "\n";
  reSizeVec(pos,num_clocks,num_epochs,iXYZ);
  reSizeVec(refpos,num_epochs,iXYZ);
  int iret=2;
    std::cout << "Here 3 ECIpos" << "\n";
  std::string ecipath = jpl_file_dir; //assume same directory
//    std::string ecipath = "/home/tday/work/new_jpl/"; //new directory for later than gpsday 16400 on 11/24/20
    std::cout << "Here 4 ECIpos" << "\n";
    std::cout << jpl_file_dir << "\n";
    std::cout << ecipath << "\n";
  if(which=="sat"||which=="both") iret=readEciSatPos(ecipath);
  if(which=="sta"||which=="both") iret=readEciStaPos(ecipath);
    std::cout << "Here 5 ECIpos" << "\n";
  return iret;
}



//******************************************************************************
int JplGpsData::readEciSatPos(std::string ecipath)
/*
170823. "new" function, work with class
2016-07-20.
Reads the ECI SATELLITE files (ecisatWWWWD.pos) for the positions/velocities.
Puts the positions into the array pos.
Note: doesn't store the velocities as these are not yet used.
Doesn't fail if the file isn't found(?)

NOTE: Relies on the fact that the ECI files have exactly 2881 entries! (i.e. no
"missed" data points) - this seems to be always true.

Units: km

INPUTS:
   ecipath  :: directory of ECI files

====== CHANGE LOG ======
160720- Started. Works.
160815- Commented out STA read-in..not needed!
160822- Reads in SAT and/or STA files. Also better error handling.
160831- Reads in the ref clock positions/velocities too.
170221- Commented out reference to vx, vy, vz. Not used now, but adds large
        amount of memory!
170320- Added "faillog" - writes out cases where program failed.
170702- Changed tempX to static!
170824- Updated to work with class.

*/
{
  if(svn.size()==0){
    std::cout<<"ERROR 757 in JplGpsData::readEciSatPos: no JPL data?\n";
    return 2;
  }
    
    std::cout << "Here 1 ECIsatpos" << "\n";

  //Read in the sat file (ecisatWWWWD.pos):
  //ECI sat files are ordered by SVN.. jpl30s files ordered by PRN..
  std::string sLine;
  std::ifstream inFile;
  std::string filename=ecipath+"ecisat"+week+day+".pos";
  std::cout << "this is filename: " << filename << "\n\n";
  if(MSC_fexists(filename)){
    inFile.open (filename.c_str());
  }else{
    if(verbose)std::cout<<"~ WARNING: ECI satellite file: ecisat"<<week<<day
                        <<".pos NOT found!\n";
    std::ofstream oFile;
    std::string filename="ECI-faillog.txt";
    oFile.open (filename.c_str(),std::ios_base::app);
    oFile<<"ecisatNotFound: "<<week<<day<<"\n";
    oFile.close();
    //goto gtENDSAT;//nb: this could be a return..(but NOT the above one!)
    return 1;
  }
    
    std::cout << "Here 2 ECIsatpos" << "\n";
    std::cout << "number of clock" << num_clocks << "\n";

  //Array to store the unsorted positions from the input ECI file:
  std::vector< std::vector< std::vector<float> > > tempX;
    std::vector< std::vector< std::vector<float> > > tempX2;
  //reSizeVec(tempX,num_clocks,num_epochs,iXYZ);
  //reSizeVec(tempX,num_clocks+1,num_epochs,iXYZ);    // added 08/20/20 for phase slip filter
  reSizeVec(tempX,num_clocks+num_satellites,num_epochs,iXYZ);    // added 08/20/20 for all Cs clock elim
    reSizeVec(tempX2,num_clocks,num_epochs,iXYZ);    // added 08/20/20

  //array to store SVN from ECI files (to match w/ jpl):
  //std::vector<int> testSVN(num_clocks);
    //std::vector<int> testSVN(num_clocks+1);       // added 08/20/20 for phase slip filter
    std::vector<int> testSVN(num_clocks+num_satellites);       // added 08/20/20 for all Cs clock elim
    std::vector<int> testSVN2(num_clocks);       // added 08/20/20
    std::cout << "this is XX: " << XX << "\n";
    std::cout << "this is YY: " << YY << "\n";
    std::cout << "this is ZZ: " << ZZ << "\n";
    
    
  std::cout << "here" << "\n\n";    
  //Loop through the file, store position data (initially sorted by SVN)
  int eciSat=0; //number of satelites in eci file
  int iEp=0;
  while(getline(inFile,sLine)){
    std::stringstream ssLine(sLine);
    std::string junk1,junk2,junk3;
    std::string stempSVN;
    double sx,sy,sz,svx,svy,svz;
    ssLine>>junk1>>stempSVN>>junk2>>junk3>>sx>>sy>>sz>>svx>>svy>>svz;
    stempSVN = stempSVN.substr(3); //remove "GPS"
    if(iEp!=num_epochs){//ignore ep2880 (second midnight)
      //populate pos array:
        
        
      tempX[eciSat][iEp][XX]=sx;
      tempX[eciSat][iEp][YY]=sy;
      tempX[eciSat][iEp][ZZ]=sz;
      iEp++;         //increment epoch
    }else{
      testSVN[eciSat]=std::stoi(stempSVN); //only need to do this once..
      iEp=0;         //re-set epoch for next sat.!
      eciSat++;      //increment # of sats.
        std::cout << "eciSat: " << eciSat <<"\n";
    }
  }//end while get-line
  inFile.close();
    
    std::cout << "total eciSat: " << eciSat << "\n";
    
    for(int i=0;i<eciSat;i++){
        std::cout <<"testSVN: " << testSVN[i] << "\n";
    }
    
    for(int i=0;i<eciSat;i++){
        std::cout << "tmpsvn[i]: " << tmpsvn[i] << "\n";
    }
    
    // added 09/30/20 to elim cs clocks properly since resorted svn here
    for(int i=0;i<eciSat;i++){
        tmpsvn[i] = 99;
    }
    
    for(int i=0;i<eciSat;i++){
        if(clk[i]=="Cs"){
            tmpsvn[i] = i;
        }
    }
    
    for(int i=0;i<eciSat;i++){
        std::cout << "new tmpsvnord[i]: " << tmpsvnord[i] << "\n";
    }
    std::sort(tmpsvnord.begin(),tmpsvnord.end());
    for(int i=0;i<eciSat - num_satellites;i++){
        tmpsvnord[i] = "0";
        //std::cout << "new2 tmpsvnord[i]: " << tmpsvnord[i] << "\n";
    }
    
    for(int i=0;i<eciSat;i++){
        std::cout << "new2 tmpsvnord[i]: " << tmpsvnord[i] << "\n";
    }
    
    double tmpord[num_satellites];
    int ixt;
    ixt = 0;
    for(int i=0;i<eciSat;i++){
        if(std::stoi(tmpsvnord[i])==0){
            continue;
        }
        tmpord[ixt] = std::stoi(tmpsvnord[i]);
        ixt++;
    }
    
    for(int i=0;i<num_satellites;i++){
        std::cout << "new3 tmpord[ixt]: " << tmpord[i] << "\n";
    }
    
    
    // end of 09/30/20
    
    
    std::cout << "Here 3 ECIsatpos" << "\n";
    
    // added on 08/20/20 to try and figure out how to eliminate a single clock
    int ix;
    ix = 0;
    ixt = 0;
    
    for(int i=0;i<num_clocks;i++){      // length 29 right now 08/20/20
        
        // this is for eliminating a single clock for the phase slip analysis
 //       if(ix == 3){
 //           ix = 4;
 //       }
        // end of phase slip clock elimination
        
        // this if statement is to get rid of all Cs clocks from network
/*
        if(ix == 0){
            ix++;
        }
        if(ix == 2){
            ix++;
        }
        if(ix == 3){
            ix++;
        }
        if(ix == 9){
            ix++;
        }
        if(ix == 10){
            ix++;
        }
        if(ix == 12){
            ix++;
        }
        if(ix == 15){
            ix++;
        }
        if(ix == 16){
            ix++;
        }
        if(ix == 17){
            ix++;
        }
*/
    // the number of if statements here equals the max number of Cs clocks in a network
 
        if(testSVN[ix] == tmpord[ixt]){
            ix++;
            ixt++;
        }
        if(testSVN[ix] == tmpord[ixt]){
            ix++;
            ixt++;
        }
        if(testSVN[ix] == tmpord[ixt]){
            ix++;
            ixt++;
        }
        if(testSVN[ix] == tmpord[ixt]){
            ix++;
            ixt++;
        }
        if(testSVN[ix] == tmpord[ixt]){
            ix++;
            ixt++;
        }
  
        // end of elimination of Cs clocks from network 09/16/20
        
        std::cout << "Here 1 in loop" << ix <<  "\n";
        for(int j=0;j<num_epochs;j++){
            
            
            
            tempX2[i][j][XX]=tempX[ix][j][XX];
            tempX2[i][j][YY]=tempX[ix][j][YY];
            tempX2[i][j][ZZ]=tempX[ix][j][ZZ];
            
            
        }
        testSVN2[i]=testSVN[ix];
        std::cout << "testSNV2: " << testSVN2[i] << "\n";
        ix++;
    }
    
    //end of new additions on 08/20/20 to kill a clock
    
    std::cout << "Here 3.5 ECIsatpos" << "\n";

  //print out warning if any problems are found.
  //if(eciSat!=num_satellites){
  if(eciSat!=num_clocks+num_satellites){
    if(verbose){
      std::cout<<"Warning 732 in readECI: jpl"<<week+day
               <<": eciSat=/=num_satellites !!!"
               <<" -- ("<<eciSat<<" =/= "<<num_satellites<<")\n";
    }
    std::ofstream oFile;
    std::string filename="ECI-faillog.txt";
    oFile.open (filename.c_str(),std::ios_base::app);
    oFile<<"eciSat=/=num_satellites: "<<week+day<<". eciSat="<<eciSat
         <<", num_satellites="<<num_satellites<<"\n";
    oFile.close();
  }
    
    std::cout << "Here 4 ECIsatpos" << "\n";
    std::cout << "eciSat: " << eciSat << "\n";

  //Order the ECI sat files into the pos array by PRN:
  //Uses the SVN to match the positions
    
    int n;
    n = 0;
  for(int i=num_receivers;i<num_clocks;i++){            // length 20 for now 08/20/20
    bool found=false;
    //for(int k=0;k<eciSat;k++){
    for(int k=0;k<eciSat;k++){        // added on 08/20/20 (length 30)
       // if(k==24){
       //     continue;       // added on 08/20/20
       // }
      //if(testSVN[k]==std::stoi(svn[i])){
      if(testSVN2[i]==std::stoi(svn[k])){       // added on 08/20/20
          std::cout << "SVN2[i]: " << testSVN2[i] << "  SVN[k]: " << svn[k] << "\n";
        found=true;
        for(int j=0;j<num_epochs;j++){
          for(int l=0;l<iXYZ;l++){
            //pos[i][j][l]=tempX[k][j][l];
        //    pos[i][j][l]=tempX[k][j][l];       // added on 08/20/20
            pos[i][j][l]=tempX2[i][j][l];       // added on 08/20/20
          }
        }
      }
        //n++;
    }
    if(!found){
      //print errors/warnings if somethine goes wrong
      if(verbose){
        std::cout<<"WARNING 749 in read ECI: jpl"<<week<<day
                 <<" missing SVN?\n";
        std::cout<<"! ---> Clock:"<<prn[i]<<"-"<<svn[i]
                 <<"("<<clk[i]<<blk[i]<<") "
                 <<"not found. Writing zeros to the positions.\n";
      }
      std::ofstream oFile;
      std::string filename="ECI-faillog.txt";
      oFile.open (filename.c_str(),std::ios_base::app);
      oFile<<"SVNnotFound: "<<week+day<<" "<<prn[i]<<"-"
           <<svn[i]<<" ("<<clk[i]<<blk[i]<<")\n";
      oFile.close();
      //in this case, just populate w/ zeros... (just to by-pass a crash)
      for(int j=0;j<num_epochs;j++){
        for(int l=0;l<iXYZ;l++){
          pos[i][j][l]=0;
        }
      }
    }
  }
    
    std::string psPath = "./results/ps13100.out";
    std::ofstream ps(psPath);
    
    for(int i=0;i<num_clocks;i++){
        
        for(int j=0;j<num_epochs;j++){
            
            for(int l=0;l<iXYZ;l++){
                
                ps << pos[i][j][l] << "\n";
                
            }
            
        }
        
    }
    
    ps.close();
    
    std::cout << "Here 5 ECIsatpos" << "\n";

  return 0;
}




//******************************************************************************
int JplGpsData::readEciStaPos(std::string ecipath)
/*
170823. "new" function, work with class
2016-07-20.
Reads in the ECI STATION files (ecistaWWWWD.pos) for the positions/velocities.
Puts the positions into the array pos.
Note: doesn't store the velocities as these are not yet used.
Doesn't fail if the file isn't found(?)

NOTE: It relies on the fact that every base station for which 30s data is
available appears in the ECI files (probably always true).
Also, relies on the fact that the ECI files have exactly 2881 entries! (i.e. no
"missed" data points) - this seems to be always true.
Also, relies on the fact that the station names are in the same (aplhabetical)
order in both the jpl and eci files..true for now, but could potentially change.

Units: km

INPUTS:
   ecipath  :: directory of ECI files

====== CHANGE LOG ======
160720- Started. Works.
160815- Commented out STA read-in..not needed!
160822- Reads in SAT and/or STA files. Also better error handling.
160831- Reads in the ref clock positions/velocities too.
170221- Commented out reference to vx, vy, vz. Not used now, but adds large
        amount of memory!
170320- Added "faillog" - writes out cases where program failed.
170702- Changed tempX to static!
170824- Updated to work with class.

====== To Do ======
*/
{
  int iret=0;

  if(prn.size()==0){
    std::cout<<"ERROR 905 in JplGpsData::readEciStaPos: no JPL data?\n";
    return 2;
  }

  //Read in station ECI files (ecistaWWWWD.pos)
  //They are ordered aplhabetically, same as jpl30s..
  //Perhaps dangerous - this might change!??
  std::string sLine;
  int eciSta=0;
  int iEp=0;
  std::ifstream inStaFile;
  std::string stafname=ecipath+"ecista"+week+day+".pos";
  if(MSC_fexists(stafname)){
    inStaFile.open (stafname.c_str());
  }else{
    if(verbose)std::cout<<"~ WARNING: ECI station files: ecista"<<week+day
                        <<".pos NOT found!\n";
    //write out to file
    std::ofstream oFile;
    std::string filename="ECI-faillog.txt";
    oFile.open (filename.c_str(),std::ios_base::app);
    oFile<<"ecistaNotFound: "<<week+day<<"\n";
    oFile.close();
    return 1;
  }

  //Put the positions into the position array, including the reference clock
  bool foundRef=false;
  while ( getline (inStaFile,sLine) ) {
    std::stringstream ssLine(sLine);
    std::string junk1,junk2,junk3;
    std::string stempSTA;
    double sx,sy,sz,svx,svy,svz;
    ssLine>>junk1>>stempSTA>>junk2>>junk3>>sx>>sy>>sz>>svx>>svy>>svz;
    if(prn[eciSta]==stempSTA){
      //only read in stations that have 30s data
      //note: this only works because the jpl30s and ECI files list the
      //station in the same order!
      if(iEp!=num_epochs){//ignore ep2880 (second midnight), populate pos array
        pos[eciSta][iEp][XX]=sx;
        pos[eciSta][iEp][YY]=sy;
        pos[eciSta][iEp][ZZ]=sz;
        iEp++;         //increment epoch
      }else{
        iEp=0;         //re-set!
        eciSta++;      //increment # of sats.
      }
    }else if(stempSTA==refprn){//Also read the REF clock!
      if(iEp!=num_epochs){//ignore ep2880 (second midnight), populate pos array
        refpos[iEp][XX]=sx;
        refpos[iEp][YY]=sy;
        refpos[iEp][ZZ]=sz;
        iEp++;         //increment epoch
      }else{
        iEp=0;         //re-set!
        foundRef=true;
      }
    }
  }
  inStaFile.close();

  if(eciSta!=num_receivers){
    //check if we read in the same # of stations are from jpl 30s files!
    //Needs to be the same, so the 'i' indices match up!
    if(verbose)std::cout
        <<"! WARNING 748: Error 655 in readECI: eciSta=/=num_receivers !!!\n";
    std::ofstream oFile;
    std::string filename="ECI-faillog.txt";
    oFile.open (filename.c_str(),std::ios_base::app);
    oFile<<"eciSta=/=num_receivers: "<<week+day<<". eciSta="<<eciSta
         <<", num_receivers="<<num_receivers<<"\n";
    oFile.close();
    iret++;
  }
  if(!foundRef){
    if(verbose){
      std::cout<<"! ERROR 757: Error 659 in readECI: Did not find Ref clock "
               <<refprn<<" !!!\n";
    }
    std::ofstream oFile;
    std::string filename="ECI-faillog.txt";
    oFile.open (filename.c_str(),std::ios_base::app);
    oFile<<"RefNotFound: "<<week<<day<<". "<<refprn<<"\n";
    oFile.close();
    iret++;
  }

  return iret;
}


//******************************************************************************
int JplGpsData::differenceData(int in_dif)
/*
170706. New function (from old "doubleDiff" function)
Differences (and zeros) the data: Single, Double, or just zero..
  in_dif=0 - Don't difference, but DO subtract the "average" (center the data)
  in_dif=1 - Single difference
  in_dif=2 - Double difference
  in_dif=3 - Mixed difference
Note that the first (one or two) data point(s) will be lost! (replaced w/ zero)
The "zeros" (from missing data) are replaced with zero.
--Note, that each "zero" in the data will affect 3 pts in the DD data
---ALL replaced with zero.

If in_dif = 3, means use "mixed" differencing: 1-order for Rb, CsIIF (and all
"others"), but 2-order for CsII and CsIIA.

INPUTS:
   in_dif           :: Order for differenceing

=== Change Log ===
170914- Uses new definition of d(2) [d(2) = d_j-2 - 2d_j-1 + d_j]
171021- Uses "mixed" differencing.
171115- Re-written to allow calling multiple times. in_dif given the 'target'
        will only difference if target>current level of differencing. And will
        only apply the appropriate level.
180319- Doesn't use weighted mean to 'zero' data

-- TO DO --
Check: I have fixed the j=0 outlier problem, but there are (less frequently)
major outliers at other positions...these still need to be deleted (?)
--at least they should be excluded from the "mean" etc.

*/
{
  if(bias.size()==0){
    std::cout<<"\nERROR 1025 in JplGpsData::differenceData: no JPL data?\n";
    return 2;
  }

  if(in_dif>3){
    if(verbose)std::cout<<"\nWARNING 1032 in differenceData: in_dif>3. "
                        <<"Data NOT differenced (just centred; in_dif=0)\n";
    in_dif=0; //just zero the data!
    differenceData(0);
    return 1;
  }
    
    std::cout << "Here in difference data" << "\n";
    std::cout << "this is num_clocks in diff data: " << num_clocks << "\n";

  //vdiff is a vector that stores level of differening.
  //NOTE: vector, that isn't re-sized untill diff is called.
  //Therefore, vdiff.size()=0 until it's called!
  if((int)vdiff.size()!=num_clocks) reSizeVec(vdiff,num_clocks);

  for(int i=0;i<num_clocks;i++){
    if(in_dif==0) break; //skip this, just 'zero' the data
    int current_dif = vdiff[i];
    //Figure out what order differening should be applied:
    //e.g., if data has already been 1-order differenced, but we want 2-order,
    // then, we need to apply 1-order differencing again!
    int target_dif = 0;
    if(in_dif==1){
      target_dif = 1;
    }else if(in_dif==2){
      target_dif = 2;
    }else if(in_dif==3){
      //default is 1-order differencing (Rb + CsIIF)*
      target_dif = 1;
      //But for CsII and CsIIA, use 2-order
      if(clk[i]=="Cs"&&(blk[i]=="IIA"||blk[i]=="II"))
        target_dif = 2;
    }
    //Already correct differening, so skip:
    if(target_dif==current_dif) continue;
    int clk_dif=0; //order of differening to apply (on top of existing)
    //Data has not been differenced before:
    if(current_dif==0) clk_dif = target_dif;
    //Data _has_ been differenced before, but not `enough'
    else if(current_dif==1 && target_dif==2) clk_dif = 1;
    if(clk_dif==0) continue; //Safety. Should never happen
    //store correct difference level
    vdiff[i] = target_dif;

    //Temporary arrays to hold data
    std::vector<double> temp(num_epochs);
    std::vector<double> tempFER(num_epochs);
    //Perform differenceing:
    if(clk_dif==2){
      //Double Difference:
      for(int j=2; j<num_epochs; j++){//do the differencing
        if(ferr[i][j-2]*ferr[i][j-1]*ferr[i][j]!=0){
          temp[j]=bias[i][j-2]-2*bias[i][j-1]+bias[i][j];
          tempFER[j]=sqrt(pow(ferr[i][j-2],2) + 4*pow(ferr[i][j-1],2)
                          + pow(ferr[i][j],2) );
        }else{
          //If any of the three points are "bad", just make it zero
          temp[j]=0;
          tempFER[j]=0;
        }//End if pts bad
      }//END epoch loop
    }else if(clk_dif==1){
      //Single-Difference
      for(int j=1; j<num_epochs; j++){//do the differencing
        if(ferr[i][j-1]*ferr[i][j]!=0){
          temp[j]=bias[i][j]-bias[i][j-1];
          tempFER[j]=sqrt(pow(ferr[i][j],2) + pow(ferr[i][j-1],2));
        }else{
           //If any of the two points are "bad", just make it zero
          temp[j]=0;
          tempFER[j]=0;
        }//End if pts bad
      }//END epoch loop
    }//END If single or double

    //replace data with differenced data:
    for(int j=0; j<num_epochs; j++){
      bias[i][j]=temp[j];
      ferr[i][j]=tempFER[j];
    }

  }//END first loop over each clock

  //"zeros" the data (subtracts the weigthed mean)
  // ?? Better outlier exlusion!! ??
  for(int i=0;i<num_clocks;i++){
    double avg=0;
    double ww=0;//sum of "weights"
    for(int j=0;j<num_epochs;j++){//Calculate the weighted mean
      if(ferr[i][j]!=0){
        //double wj=1/pow(ferr[i][j],2);
        double wj=1.; //don't use weigts!
        avg+=bias[i][j]*wj;
        ww+=wj;
      }
    }
    avg=avg/ww;
    for(int j=0;j<num_epochs;j++){//Subtract the mean!
      if(ferr[i][j]!=0){
        bias[i][j]=bias[i][j]-avg;
      }
    }
  }
    
    
        std::string gvalPath = "./results/g10val13100.out";
        std::ofstream gval(gvalPath);
    
        std::string gval2Path = "./results/g10ferr13100.out";
        std::ofstream gval2(gval2Path);
        
        for(int i=0;i<num_clocks;i++){
          for(int j=0;j<num_epochs;j++){
              //if(i==10){
    //              continue;
    //          }
                  gval <<  bias[i][j] << "\n";
                  gval2 <<  ferr[i][j] << "\n";
    //          ferr[nit][j] = ferr[i][j];
              //}
              
          }
    //        nit = nit + 1;
        }
        
        gval.close();

  return 0;
}



//******************************************************************************
int JplGpsData::polynomialDetrend(int poly_order, int weighted)
/*
Uses the GSL library to perform an any-order weighted polynomial fit (and
subtraction) to given data.
Weighting can be turned off (weighted=0 => not weighted)
Does not include any 'missing' data points (marked with a sigma=0) in the
fit [gives them zero weight, even if weigting is switched off).
::polywfit
Is a wrapper contained in a sepperate file, for the GSL library
INPUT---
  poly_order:  order of the polynomial
  weighted: weighting switch
  dData:   Input clk data array: dData[clock][epoch][type]. type: 0=bias; 1=sig
OUTPUT---
  dR2:   R^2 value for each clock: dR2[clock]

NB: Only removes ploynomial from Satellite clocks, not base-stations!

====== CHANGE LOG ======
160226- Replaced my 2-order poly fit with that of GSL:
160809- Gave first data point a weight of zero, since it is often an outlier!
        (note, this is true even for non-weighted fit! Good)
160810- Writes 0's to the bias of clocks that have a "zero" flag for the FER..
        ..I don't understand why this wouldn't already be the case though!
160811- Ignores (assigns 0 weight) to first/last 5 points!
160822- does NOT ignore first/last points. They are fine
170702- Changed coef to std::vector
170827- Works with class.
170914- Only removes poly from sat. clocks, not base-stations
*/
{
  if(bias.size()==0){
    //nb: bias.size should = num_clocks
    std::cout<<"ERROR 1144 in JplGpsData::polynomialDetrend: no JPL data?\n";
    return 2;
  }

  //before this, has .size()=0
  reSizeVec(dR2,num_clocks); //r^2 value

  for(int i=0;i<num_clocks;i++){//Base-station clocks?
  //for(int i=num_receivers;i<num_clocks;i++){//Only satelite clocks
    std::vector<double> dx(num_epochs); //"data points" (just epochs)
    std::vector<double> dy(num_epochs); //Data (observations)
    std::vector<double> dw(num_epochs); //Weights
    //chi^2 and R^2 values:
    //(don't use chi^2 for now, may later?)
    double chisqi,R2i;
    // one more coef than the order:
    const int coefSize=poly_order+1;
    std::vector<double> coef(coefSize);
    //prep the data for polywfit
    for(int j=0;j<num_epochs;j++){
      dx[j]=j;
      dy[j]=bias[i][j];
      double sig=ferr[i][j];
      if(sig!=0){
        if(weighted!=0){
          dw[j]=1/pow(sig,2);
        }else{
          dw[j]=1;
        }
      }else{
        //Ignore the "missed" data points
        dw[j]=0;
      }
    }

    //perform the actual fit:
    //GSL wrapper, contianed in "mathematicsFunctions.cpp"
    MFS_polywfit(num_epochs,poly_order,dx,dy,dw,coef,chisqi,R2i);

    //Output fit-parameters R^2 and chi^2 of polywfit
    dR2[i]=R2i;
    //dChi2[i]=chisqi;

    //detrends the data
    for(int j=0;j<num_epochs;j++){
      double sig=ferr[i][j];
      double funj=0;
      for(int k=0;k<poly_order+1;k++){
        funj+=pow(dx[j],k)*coef[k]; //value of the best-fit function
      }
      if(sig!=0){//don't de-trend the missing data!
        bias[i][j]=dy[j]-funj;
      }else{
        if(bias[i][j]!=0){
          bias[i][j]=0;
          //This happens if data points on either side is missing, but this one
          //isn't! When I try to extrapolate the errors, it still ends up zero.
          //This data point is useless..so remove it?
        }
      }
    }
  }

  return 0;
}


//******************************************************************************
int JplGpsData::subtractWeightedMean(int dif)
/*
171115.
Function that forms a single time-series that is the weighted mean of all
available first-order differenced time series'. The weights are the inverse
standard deviations. This time-series should then approximate the "raw" time
series (first-order diff'd) of the reference clock.
This is then subtracted from each of the clocks, in order to "remove" the
influence of the reference clock and hence remove any cross-clock correlations.
Actually, what is subtracted is the weighted mean of all OTHER clocks

INPUT:
  dif   ::  the level the data should be differenced at (1=1, 2=2, 3=mixed)
            0: not allowed (will be the same as '1')

This should be called _instead_ of differenceData and calculateStdDev.
Will NOT work if differenceData has already been called.
*/
{

  if(vdiff.size()>0){
    std::cout<<"\nERROR 1395 in JplGpsData::subtractWeightedMean. This must"
             <<" be called instead of differenceData!\n";
    return 2;
  }
  if(dif==0){
    std::cout<<"\nWARNING 1418 in JplGpsData::subtractWeightedMean. dif=0. "
    <<" Changed to 1.\n";
    subtractWeightedMean(1);
    return 1;
  }

  //First, perform first-order differencing on the data
  differenceData(1);
  //then, calculate the sds
  calculateStdDev();

  //Form the 'weighted norm constant'
  double W=0;
  for(int i=0; i<num_clocks; i++){
    if(sdev[i]==0) continue; //? just avoid a crash?
    W+=(1/sdev[i]);
  }
  if(W==0){
    std::cout<<"\nERROR 1435 in JplGpsData::subtractWeightedMean. W=0\n";
    return 2;
  }

  //Form the weighted average
  std::vector<double> dbar;
  for(int j=0; j<num_epochs; j++){
    double temp_dbar=0;
    double temp_W = W;
    for(int i=0; i<num_clocks; i++){
      if(ferr[i][j]==0){
         temp_W-=(1/pow(sdev[i],2)); //get rid of this from the 'weight'
         continue;
      }
      if(sdev[i]==0) continue; //just to avoid crash..
      temp_dbar += bias[i][j]/pow(sdev[i],2);
    }
    dbar.push_back(temp_dbar/temp_W);
  }

  //Subtract the weighted mean from each clock
  for(int i=0; i<num_clocks; i++){
    for(int j=0; j<num_epochs; j++){
      if(ferr[i][j]==0) continue; //skip bad points
      //bias[i][j]-=dbar[j];
      //Subtract the weighted mean of all OTHER clocks:
      bias[i][j] -= dbar[j] + (dbar[j] - bias[i][j])/(W*pow(sdev[i],2)-1);
      //XXX DOUBLE CHECK! XXX
    }
  }

  //Difference the data to appropriate level
  // and re-calculate std deviation.
  differenceData(dif);
  calculateStdDev();

  return 0;
}


//******************************************************************************
double JplGpsData::calculateStdDev(void)
/*
170727.
Calculates the standard deviations for each satellite and receiver clock,
and stores them (in the sdev array).
Note: it doesn't make much sense to apply this function to non-differenced data,
If called after just 0th order differencing, it will throw a warning.
If called before any differencing, it will fail.
Function also returns the average of this s.d.
--> NB: Returns avg sd of ONLY the sat clocks! (not recievers!)
  *Note: ignores "bad" epochs, as in missing data points.
  Does NOT skip any outliers.
  ALSO: assumes data is centered on zero!
XXX return thing doesn't work with 'double'...
*/
{
    
    std::cout << "Here in calc std" << "\n";

  if(bias.size()==0){
    std::cout<<"ERROR 1222 in JplGpsData::calculateStdDev: no JPL data?\n";
    return 2;
  }

  // If "differenceData" hasn't been called, the data
  if(vdiff.size()==0){
    std::cout<<"ERROR 1239 in JplGpsData::calculateStdDev: "
             <<"Data has not been differenced at all. s.d. is meaningless\n"
             <<"--Fatal error, calculateStdDev exited\n";
    return 2;
  }

  //If data has been zero'd, but not differenced, also not good.
  if(vdiff[0]==0){
    // XXX also check for polyremoval? Without, complete junk.
    // But, with, maybe useful sometimes (for tests etc.)
    if(verbose){
      std::cout<<"WARNING 1246 in JplGpsData::calculateStdDev: "
               <<"Data has been zerod, but not differenced.\n"
               <<"--Program will continue, but results likely meaningless!\n";
    }
  }

  //Size the array that stroes the clock s.d.
  //The size of this array should be zero unless filled!
  reSizeVec(sdev,num_clocks);

  int good_sat_clocks=0;
  double avg_sat_sd=0;
  for(int i=0; i<num_clocks; i++){
    double sd=0;  //hold the intermediate sd calculation
    int good_epochs=0; //denominator for sd calculation (ignore 'bad' epochs)*
    for(int j=0; j<num_epochs; j++){
      if(ferr[i][j]==0) continue;
      good_epochs++;
      sd+=pow(bias[i][j],2);
    }
    if(good_epochs!=0){
      sd=sqrt(sd/good_epochs);
      if(blk[i]!="AR"){
        //only include the satellite clocks in the average sd! (XXX NOT checked)
        good_sat_clocks++;
        avg_sat_sd+=sd;
      }
    }
    sdev[i]=sd;
  }
    
    for(int i=0; i<num_clocks; i++){
        
        std::cout << "i: " << i << "sdev[i]: " << sdev[i] << "\n";
        
    }

    // this is included on 10/16/20 tday
    std::cout << "This is avg_sat_sd: " << avg_sat_sd/good_sat_clocks << "\n";

  return avg_sat_sd/good_sat_clocks;
}



//******************************************************************************
int JplGpsData::formAutoCorrelation(std::string path_to_acf,
                                    std::string acf_label)
/*
171013.
Public function.
Reads in the ACF functions for each clock.
If it cannot find a suitable file, will first try 'av' (the average SVN) and/or
Hmas ref. (average of the Hmaser reference clocks).
Failing that, it will just calculate the ACF using todays data.
input: path_to_acf is optional. It is "calculateACF" by default.
Therefore, if this function is called without any arguments, the ACF will be
calculated, instead of read from the input file.
*/
{
  reSizeVec(acf,num_clocks,acf_points);
    
    std::cout << "number of clocks in read ACF: " << num_clocks << "\n";
    /* need to make sure the correct SVN clocks are being read
     should be a simple thing as was done with the previous
     Cs clock removals as above
     */
    
    // added for Cs elim
    int nt;
    nt = 0;
    //

  //Loop through each clock, get/form the ACF:
// commented out on 09/18/20 for elimination of Cs clocks from network
//  for (int i=0; i<num_clocks; i++){
     

  // if we are doing ACF from file we have to add the extra 9 to num_clocks 09/20/20
  // if we are doing ACF from data we do not include the extra nine since already eliminated!
  // num_satellites represents number of Cs clocks
  for (int i=0; i<num_clocks+num_satellites; i++){
//    acf[i][0]=1; //Always true by default // commented out on 09/18/20
    acf[nt][0]=1; //Always true by default
      
//    std::cout << "this is SVN in read ACF : " << svn[i] << "\n";
      std::cout << "this is nt: " << nt <<" in read ACF this is acf[nt] : " << acf[nt][0] << "\n";

    //Calculate the ACF, or read it in from file?
    if(path_to_acf=="calculateACF"){
      calculateAcfFromData(i);
      continue;
    }

    // if we are doing ACF from data or not eliminating Cs clocks we use this 09/20/20
    //If not Rb, Cs, or H clock, assume white:
 //   if(clk[i]!="Cs" && clk[i]!="Rb" && clk[i]!="H") continue;  // commented out for Cs clock elim ACF from file
      
      // taking out the Cs clocks from network 09/18/20
      // if we are doing ACF from file we do this 09/18/20
      if(clk[i]!="Rb" && clk[i]!="H") continue;
    //if(clk[i]!="Cs" || clk[i]!="Rb" || clk[i]!="H") continue;
      
      std::cout << "this is SVN in read ACF : " << svn[i] << "\n";

    //Work out ACF file name:
    //std::string prefix = path_to_acf+"/"+clk[i]+blk[i]+"-";       // commented out 09/24/20 tday
    std::string prefix = path_to_acf+clk[i]+blk[i]+"-";             // commented in 09/24/20 tday
    std::string suffix = "-"+acf_label+".acf";
    if(acf_label=="na"||acf_label=="") suffix=".acf";
    // File name for the ACF file:
    std::string acf_fname = prefix+svn[i]+"-"+refname+suffix;

    //Check if ACF file exists. If not, try the average SVN
    if(!MSC_fexists(acf_fname)){
      acf_fname = prefix+"av-"+refname+suffix;
      //Check if /that/ file exists. If not try real SVN, with 'Hmas' reference
      if(!MSC_fexists(acf_fname)){
        // If a sat clock is reference, can't use Hmas! Calculate ACF!
        if(refblk!="AR"){
          calculateAcfFromData(i);
          continue;
        }
        //Try with Hmas:
        acf_fname = prefix+svn[i]+"-Hmas"+suffix;
        //Check that file. If not, try averags SVN, and 'Hmas' ref.
        if(!MSC_fexists(acf_fname)){
          acf_fname = prefix+"av-Hmas"+suffix;
          // If still can't find file, then just assume white noise:
          if(!MSC_fexists(acf_fname)){
            calculateAcfFromData(i);
            continue;
          }
        }
      }
    }//End if full acf path exists
      
      std::cout << "this is acf_fname in read ACF: " << acf_fname << "\n";

    int id = vdiff[i];
    //readAcfFile(acf_fname,i,id);      // commented out 09/22/20 by tday
    // this is added 09/24/20 by tday
    id = 1;   // we are getting junk numbers, this represents just an Rb network.
    //readAcfFile(acf_fname,nt,id);
    readAcfFile(acf_fname,i,id,nt);
      
      std::cout << "this is acf[nt][j] after read: " << acf[nt][1] << "\n\n";
	std::cout << "here nt: " << nt << "\n";      
      nt++;
	std::cout << "here 2 nt: " << nt << "\n"; 
// some addtion on 10/12/20
	std::cout << "here n_clocks: " << num_clocks << "\n";
	std::cout << "here is i: " << i << "\n"; 
    if(nt+1>num_clocks){
	break;
    }	
// end of addition on 10/12/20
  }//END loop over clocks
    
    for(int i=0;i<num_clocks;i++){
        for(int j=0;j<acf_points;j++){
            std::cout << "i: " << i <<"  this is all acf[i][j]: " << acf[i][j] << "\n";
        }
        std::cout <<"\n";
    }

  return 0;
}

//******************************************************************************
int JplGpsData::calculateAcfFromData(int ic)
/*
171017.
Quick routine to calculate ACF using the clock data for given day.
Calculates for a single clock [with index ic]
This is not averaged.
*/
{
  //The standard deviations need to exist.
  if(sdev.size()==0) calculateStdDev();

  //Starting point (1 for 1st order diff'd data, 2 for 2nd)
  int jmin=vdiff[ic];

  //Calculate the ACF using todays data
  acf[ic][0]=1; //should already be set, but may as well.
  for(int t=1; t<acf_points; t++){
    double temp_acf=0;
    int jmax = num_epochs - t;
    for(int j=jmin; j<jmax; j++){
      temp_acf += bias[ic][j]*bias[ic][j+t];
    }
    acf[ic][t] = temp_acf / ( pow(sdev[ic],2)*(num_epochs-t) );
  }

  return 0;
}


//******************************************************************************
int JplGpsData::readAcfFile(std::string acf_fname, int ic, int id, int nt)
/*
171013.
Private function.
Reads in a single ACF from file.
Separated from "formAutoCorrelation" just for clarity.
  ic: clock index (which array position to write ACF to)
  id: differencing level for this clock/ACF
  nt: correct clock index when removing the Cs clocks from network 09/30/20 tday
*/
{
    
/*
 need to change the ic value in this script to match what the nt value is from the formAutocorellation function above
 so we should probably pass that value into this function and then change the acf at bottom of this function to be acf[nt][jj] and not acf[ic][jj]
 */
    
    std::cout << "this is snv[i] in readAcfFile: " << svn[ic] << "\n";
    std::cout << "this is nt in readAcfFile: " << nt << "\n\n";
    
    
  //open the ACF file, and read in the relevant info
  std::ifstream ifs_acf;
  ifs_acf.open (acf_fname.c_str());
  std::string s_line;

    
  int jj=-1; //numerates position of ACF. -1=>sd
  while(getline(ifs_acf,s_line)){

    std::stringstream ss_test(s_line);
    std::string s_temp;
    ss_test>>s_temp;
    if(s_temp=="#")continue; //"#" marks header/junk lines

    if(jj==-1){
      //read in the sd's from the file:
      double sd1,sd2;
      std::stringstream ssin(s_line);
      ssin>>sd1>>sd2;
      // NB: don't use these now, but might in future?
    }else{
      //read in the ACF
      double acf0,acf1,acf2;
      std::stringstream ssin(s_line);
      ssin>>acf0>>acf1>>acf2;
      //if(id==1) acf[ic][jj]=acf1;
      //if(id==2) acf[ic][jj]=acf2;
      if(id==1) acf[nt][jj]=acf1;
      if(id==2) acf[nt][jj]=acf2;
        std::cout << "this is id: " << id << "  " << "this is acf1: " << acf1 <<"\n";
        std::cout << "this is acf[nt][jj]: " << acf[nt][jj] << "\n";

    }
          jj++; //increment which ACF point we are up to
      
      std::cout << "this is jj: " << jj << "\n";
      

    //only read in "num_acf_points" points from ACF file
    // If acf_points is larger than file size, rest just zeros.
    if(jj>=acf_points)break;

  }//END while getline

  ifs_acf.close();
  return 0;
}


//******************************************************************************
int JplGpsData::calculateCrossClockCorrelation(int max_lag)
/*
171017.
Calculates the cross-clock correlations, defined:
  ccc_ik = < d^i d^k >
The CCC should be the same for each pair of clocks! Therefore, only need to
calculate it once!
-Averages over all clock pairs, for better numerics
-Can calculate non-zero lag terms [up to max_lag].
Note: Includes DOES NOT include the negative sign, or standard deviations!
These are needed for the covariance matrix.

XXX NOTE: Since this depends ONLY on reference clock, should probably include
this instead in 'processNoise', and just read it in??
-> maybye not: better averageing that way, but the clock may change over time!

=== Change Log ===
171112- Starting to add non-zero lag terms...XXX not tested!
180118- Updating. CCC is same for each clock pair - so, can just store
        one number, and average over all clocks! Better efficiency
*/
{
  //ensure array is clear, and make it correct size.
  //nb: max_lag is typically 0; we neglect non-zero lag terms here
  b0.clear();
  reSizeVec(b0,max_lag+1);

  //Loop over each (unique) pair of clocks.
  //Calculate the CCC between each pair, and average them.
  int num_clock_combos=0;
  for(int i=0; i<num_clocks; i++){
    for(int k=i+1; k<num_clocks; k++){
    //for(int k=i+1; k<=i+1; k++){
      num_clock_combos++;
      for(int l=0; l<=max_lag; l++){ //loop over lag.
        //std::cout<<i<<" "<<k<<" HELLO\n\n";
        double temp_ccc=0;
        int good_epochs=0;
        //average over epochs:
        for(int j=0; j<num_epochs-l; j++){
          if(ferr[i][j]*ferr[k][j+l]==0) continue; //skip 'bad' points
          good_epochs++;
          temp_ccc += bias[i][j]*bias[k][j+l];
          //std::cout<< bias[i][j]<<" "<<bias[k][j+l]<<" "<<temp_ccc<<"\n";
        }// j (epoch)
        //std::cout<<"\nGood epochs:"<<good_epochs<<"\n";
        if(good_epochs!=0)
          temp_ccc = temp_ccc/(2*good_epochs); /// XXX ???
          //XXX Need to have the '2' here, but I don't know why!????? XXX
          //temp_ccc = temp_ccc/(good_epochs);
        b0[l] += temp_ccc;
      }// l (lag)
    }//k
  }//i
  //devide by number of clock combos (to get average)
  for(int l=0; l<=max_lag; l++){
    //std::cout<<"\n XX1: "<<sqrt(fabs(b0[l]))<<" "<<num_clock_combos<<"\n";
    b0[l]/=num_clock_combos;
    std::cout << "this is b0: " << b0[l] << "\n";
    //std::cout<<"\n XX2: "<<sqrt(fabs(b0[l]))<<" "<<num_clock_combos<<"\n";
  }
    
    

  //b0[0] = sdev[0]*sdev[0];

  return 0;
}

//******************************************************************************
int JplGpsData::avgNetCCC( int day, int week, std::string path_to_avgccc)//for now adding day and week for file.
/* FUNCTION: avgNetCCC
 this function now see's if our day in is which network and we should use this instead of using the calculate_crossclockcorrelation function
    03/11/2020, tdaykin
*/
{
    
    
    
    // need to make a variable for the week + day so I can use to check
    int GPSday;
    std::string tempGPSday;
    
    int num_networks_minus1 = 102;
    
    //ensure array is clear, and make it correct size.
    b0.clear();
    reSizeVec(b0,num_networks_minus1);
    
    //CCC.clear();
    //netEnd.clear();
    reSizeVec(CCC,num_networks_minus1);
    reSizeVec(netEnd,num_networks_minus1);
    
    avgNetCCCread(num_networks_minus1, path_to_avgccc);
    
    tempGPSday = std::to_string(week) + std::to_string(day);
    GPSday = std::stoi(tempGPSday);
    
    //std::cout << GPSday << "\n\n";
    //std::cout << netEnd[0] << "\n\n";
    //std::cout << CCC[0] << "\n\n";
    //std::ofstream cccmat("./e16results/avgcrossclocksearch.txt");
    
    for(int i=1; i<num_networks_minus1; i++){
        
        if(GPSday <= netEnd[i-1]){
            b0[0] = CCC[i-1];
            break;
        }
        if(GPSday > netEnd[i-1] && GPSday <= netEnd[i]){
            b0[0] = CCC[i];
            break;
        }
        
        
    }
    //cccmat << b0[0] << "\n";
    
    //cccmat.close();
    
    return 0;
}


//******************************************************************************
int JplGpsData::avgNetCCCread(int num_networks, std::string path_to_avgccc)//for now adding day and week for file.
/* FUNCTION: avgNetCCCread
 this function is to read in each networks CCC value and end date
 and so that you may determine which network is needed for the day that we are currently on
 03/11/2020
*/
{

    double tempCCC[num_networks];
    double tempenddate[num_networks];


    std::ifstream indata(path_to_avgccc);

    for(int i=0;i<num_networks;i++){
            indata >> tempenddate[i] >> tempCCC[i];
            netEnd[i] = tempenddate[i];
            CCC[i] = tempCCC[i];
    }
    
    return 0;
    
}

//******************************************************************************
int JplGpsData::formRelativeCouplings(std::string interaction)
/*
171024.
Function that forms an array of "effective" (relative) coupling strengths for
each clock. Depends on which interaction is assumed (i.e. alpha, m_e etc.).
For now: Normalised to Rb (i.e.,
  keff(Rb) = 1
  keff(Cs) = k_Cs / k_Rb, etc.
Note: the reference clock entry is stored last - in the [num_clocks] position
NOTE: two ways to parameterise m_p variation. Either, absorb into m_e term (as
in our Nat.Comm. paper), or say dmp/mp = +0.05 dmq/mq [Stadnik/Flambaum].
Probably the latter is more correct!

Input:
  interaction ::  Can be: "alpha", "me", "mq", or ""
                  --blank by defult (which assumes each clock the same)
*/
{
  //Re-size the array:
  reSizeVec(keff,num_clocks);

  //index identifiers (to avoid mistakes)
  int ih=0, ics=1, irb=2; // clocks: H, Cs, Rb
  int ia=0, ime=1, imq=2; // interactions: alpha, m_e, m_q

  //Which interactions? "-1" means assume all clocks react same way!
  int x=-1;
  if(interaction=="alpha")    x=ia;   //Coupling to alpha
  else if(interaction=="me")  x=ime;  //Coupling to electron mass*
  else if(interaction=="mq")  x=imq;  //Coupling to quark mass

  //Array that holds the sensitivity coeficients.
  //A constant.
  //Order: H, Cs, Rb ; alpha, m_e, m_q
  // Note: there are two possibilities in regard to m_p..
  // // Like in our Nat.Comm. paper: 2m_e/m_p:
  // double kx[3][3]={4.0,   1., -0.1,     //H
  //                  4.83,  1., 0.002,    //Cs
  //                  4.34,  1., -0.019};  //Rb
  ////// ~~ OR ~~
  // Assume dmp/mp = +0.05 dmq/mq. Middle entry mow just m_e !
  // I think this is more valid!
  double kx[3][3]={4.0,   2., -0.15,    //H
                  4.83,  2., -0.048,   //Cs
                  4.34,  2., -0.069};  //Rb

  for(int i=0; i<num_clocks; i++){
    if(x<0){
      //assuming all clocks interact the same way. Just fill with 1.
      keff[i]=1.;
      continue;
    }
    if(clk[i]=="Rb")      keff[i] = 1.;
    else if(clk[i]=="Cs") keff[i] = kx[ics][x]/kx[irb][x];
    else if(clk[i]=="H")  keff[i] = kx[ih][x]/kx[irb][x];
    else keff[i] = 1.; // Or zero??
  }

  //For the reference clock:
  float refkeff=1;
  if(x>=0){
    if(refclk=="Rb")      refkeff = 1.;
    else if(refclk=="Cs") refkeff = kx[ics][x]/kx[irb][x];
    else if(refclk=="H")  refkeff = kx[ih][x]/kx[irb][x];
    else if(refclk=="USN3")  refkeff = kx[ih][x]/kx[irb][x];
    else refkeff = 1.; // Or zero??
  }
  keff.push_back(refkeff);

  return 0;
}

//******************************************************************************
bool JplGpsData::swapReference(std::string clock, std::string block,
                               bool by_svn)
/*
170827.
20160229 New program.
Sets one of the satellite clocks as the reference, instead of the base station
(by subtracting the time-series').
Updates the reference clock name data, and position data to reflect this.
At first, will look for the ``best'' clock of the correct type (Cs/Rb)
[and, optionally, block], that have no missed data points.
[If all clocks have missing points, it drops this condition].
By default, it decides the ``best'' clock by the best standard deviation.
Alternatively, can do it by highest SVN (newest satellite).
If called with `best s.d.' option, but s.d. not yet calculated, will revert to
doing by SVN.
Finally, if still no suitable clocks are found as reference, it will write zeros
to all the clocks - they should not be used in any analysis!

INPUT:
  clock     :: which clock type to use as reference? (Rb, Cs, or 'all')
  block     :: which block to use? II,A,R,F,'all'
  by_svn    :: False by dafault. If true, picks largest svn. If false, picks
               smallest standard-deviation

=== Change Log ===
160824- Writes zeros to the 'best' clock FER (to mark it not for use)
160825- Fixed rare error that occurs when no clocks meet criteria to be ref.
160830- New way to store "new" ref clock info. & added sRefBLK
160831- Large changes. Changed name of program.
        Killed abilty to subtract average.
        Finds best clock in a smarter way - optimized for swapping reference clk
        Now, also updates the ref clock positions!
160920- Also updates the formal errors! XXX Are they dependent??
161122- Changed dR2[i]>R2min to dR2[i]>=R2min (for case where both zero)
170727- Works with class
===
TO-DO:

Should the new formal error be from quadrature??


*/
{
  int best_index=-1;
  bool found=false;

  if(clock=="any")clock="all"; //common error
  if(block=="any")block="all"; //common error

  //if calculateStdDev hasn't been called, "sd" method won't work.
  if(sdev.size()==0&&!by_svn){
    std::cout<<"WARNING 1340 in JplGpsData::swapReference: "
             <<"Tried to identify best clock by sd, but sd not yet calcd!\n"
             <<"--Will instead identify by best SVN, but there is probably"
             <<" a mistake somewhere!\n";
    by_svn=true;
  }

  std::vector<bool> any_missing(num_clocks);
  bool all_missing=true;
  for(int i=num_receivers; i<num_clocks; i++){ // XXX only sat clocks?
    if(clk[i]!=clock && clock!="all") continue;
    if(blk[i]!=block && block!="all") continue;
    for(int j=2;j<num_epochs-2;j++){
      if(ferr[i][j]==0){
        any_missing[i]=true;
        break;
      }else{
        any_missing[i]=false;
      }
    }
    if(!any_missing[i])all_missing=false;
  }

  if(by_svn){
    int best_svn=0;
    for(int i=num_receivers; i<num_clocks; i++){ // XXX only sat clocks?
      //check if clock and block are correct:
      if(clk[i]!=clock && clock!="all") continue;
      if(blk[i]!=block && block!="all") continue;
      if(any_missing[i] && !all_missing) continue;
      int this_svn=std::stoi(svn[i]);
      if(this_svn>best_svn){
        best_svn=this_svn;
        best_index=i;
        found=true;
      }
    }
  }else{
    double best_sd=1000.;
    for(int i=num_receivers; i<num_clocks; i++){ // XXX only sat clocks?
      //check if clock and block are correct:
      if(clk[i]!=clock && clock!="all") continue;
      if(blk[i]!=block && block!="all") continue;
      if(sdev[i]==0)continue;
      if(any_missing[i] && !all_missing) continue;
      if(sdev[i]<best_sd){
        best_sd=sdev[i];
        best_index=i;
        found=true;
      }
    }
  }

  // Subtract the 'best' clock from the rest, in order to swap the reference
  if(found){

    //loop to store "best" clock...otherwise overwite it first..
    std::vector<double>temp(num_epochs);
    std::vector<double>tempFER(num_epochs);
    for(int j=0;j<num_epochs;j++){
      temp[j]=bias[best_index][j];
      tempFER[j]=ferr[best_index][j];
    }
    //Subtract best clock from all satellite clocks!
    for(int i=num_receivers;i<num_clocks;i++){
      for(int j=0;j<num_epochs;j++){
        if(ferr[i][j]==0)continue;
        if(ferr[best_index][j]!=0){
          bias[i][j]=bias[i][j]-temp[j];
          ferr[i][j]=sqrt(pow(ferr[i][j],2)+pow(tempFER[j],2));
        }else{
          //If "best" clock has zeros - all data useless->zero
          bias[i][j]=0;
          ferr[i][j]=0;
        }
      }
    }
    //"Mark" the "best clock" with zeros for FER
    for(int j=0;j<num_epochs;j++){
      //should already be the case, but may as well:
      bias[best_index][j]=0;
      ferr[best_index][j]=0;
    }
    //update ref clock info
    refprn=prn[best_index];
    refsvn=svn[best_index];
    refclk=clk[best_index];
    refblk=blk[best_index];
    refname=refclk+refblk; //+"-"+refsvn+refprn;

    //update Ref clock position/velocity values
    //but only if they exist/have been initialised!
    if(refpos.size()*pos.size()>0){
      for(int j=0;j<num_epochs;j++){
        for(int k=0;k<iXYZ;k++){
          refpos[j][k]=pos[best_index][j][k];
        }
      }
    }

   //Call differenceData function with 0, this just "zero's" the data
   // (i.e., subtracts the mean). Note: technically, this isn't needed here,
   // the data already has a very small mean. But can't hurt?
   differenceData(0);

   //Now that I have changed the reference clock, I need to re-calculate
   // the standard deviations of each clock:
   calculateStdDev();

  }else{
    //std::cout<<"Didn't find a good ref!\n";
    //if no suitable reference clock, kill this entire day!
    for(int i=num_receivers;i<num_clocks;i++){
      for(int j=0;j<num_epochs;j++){
        bias[i][j]=0;
        ferr[i][j]=0;
      }
    }
    refprn="**";
    refsvn="00";
    refclk="**";
    refblk="none";
    refname="none";
  }

  return found;
}



//******************************************************************************
/*
Helper functions for the binary read in/out
*/
template<typename T>
int binary_write(std::fstream& stream, const T& value){
    stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return 0;
}
template<typename T>
int binary_read(std::fstream& stream, T& value){
    stream.read(reinterpret_cast<char*>(&value), sizeof(T));
    return 0;
}
int binary_str_write(std::fstream& stream, std::string& value){
    binary_write(stream,value.length());
    stream.write(value.c_str(), value.length());
    return 0;
}
int binary_str_read(std::fstream& stream, std::string& ovalue){
    size_t temp_len;
    binary_read(stream,temp_len);
    char* value = new char[temp_len+1];
    stream.read(value, temp_len);
    value[temp_len] = '\0'; //null 'end of string' character
    ovalue = value;
    delete [] value;
    return 0;
}
template<typename T>
int binary_rw(std::fstream& stream, T& value, bool write){
    if(write) binary_write(stream, value);
    else binary_read(stream, value);
    return 0;
}
template<typename T>
int binary_str_rw(std::fstream& stream, T& value, bool write){
    if(write) binary_str_write(stream,value);
    else binary_str_read(stream,value);
    return 0;
}

//******************************************************************************
int JplGpsData::binaryJpl30s(std::string path, int in_jplweek, int in_jplday,
  bool pos, bool write)
/*
180316.
Program to read/write the JPL30s files in/out of binary file.
Will return 1 if asked to read a data file that doesn't exist.
Note: fill NOT give error message if asked to read position file that doesn't
exist
*/
{

  //attempt to make the directory, and check it exists
  MSC_execute("mkdir -p "+path);
  if(!MSC_direxists(path)){
    //The "create directory (mkdir)" command above only works on linux pcs:
    std::cout<<"ERROR 2060 in JplGpsData::binaryJpl30s: Directory: "<<path
             <<" doesn't exist. You have to create it manually.\n";
    return 2;
  }

  int wwwd=in_jplweek*10+in_jplday;
  std::string fname1=path+"/bin"+MSC_padIntString(wwwd,5)+".clk";

  if((!write)&&(!MSC_fexists(fname1))) return 1;
  binaryReadWriteData(fname1,write);
  if(pos){
    std::string fname2=path+"/bin"+MSC_padIntString(wwwd,5)+".pos";
    binaryReadWritePositions(fname2,write);
  }

  return 0;
}

//******************************************************************************
int JplGpsData::binaryReadWriteData(std::string fname, bool write)
/*
180316.
Private function.
Reads or writes clock data (and names etc) to binary file
*/
{

  std::fstream iof;
  if(write) iof.open(fname,std::ios_base::out|std::ios_base::binary);
  else      iof.open(fname,std::ios_base::in|std::ios_base::binary);
  //Intitial sizes etc.
  binary_rw(iof,num_receivers,write);
  binary_rw(iof,num_satellites,write);
  binary_rw(iof,num_clocks,write);
  binary_rw(iof,num_epochs,write);
  binary_str_rw(iof,date,write);
  binary_str_rw(iof,week,write);
  binary_str_rw(iof,day,write);
  //Reference clock names:
  {
    binary_str_rw(iof,refprn,write);
    binary_str_rw(iof,refsvn,write);
    binary_str_rw(iof,refclk,write);
    binary_str_rw(iof,refblk,write);
    binary_str_rw(iof,refname,write);
  }
  if(!write){
    //if reading data, re-size necisary arrays
    reSizeVec(bias,num_clocks,num_epochs);
    reSizeVec(ferr,num_clocks,num_epochs);
    reSizeVec(prn,num_clocks);
    reSizeVec(svn,num_clocks);
    reSizeVec(clk,num_clocks);
    reSizeVec(blk,num_clocks);
  }
  //Clock names and data
  for(size_t i=0; i<bias.size(); i++){
    {
      binary_str_rw(iof,prn[i],write);
      binary_str_rw(iof,svn[i],write);
      binary_str_rw(iof,clk[i],write);
      binary_str_rw(iof,blk[i],write);
    }
    for(size_t j=0; j<bias[i].size(); j++){
      binary_rw(iof,bias[i][j],write);
      binary_rw(iof,ferr[i][j],write);
    }
  }
  iof.close();

  return 0;
}

//******************************************************************************
int JplGpsData::binaryReadWritePositions(std::string fname, bool write)
/*
180316.
Private function.
Reads or writes clock positions to binary file
*/
{
  std::fstream iof;
  if(write) iof.open(fname,std::ios_base::out|std::ios_base::binary);
  else      iof.open(fname,std::ios_base::in|std::ios_base::binary);

  if(!write){
    reSizeVec(refpos,num_epochs,iXYZ);
    reSizeVec(pos,num_clocks,num_epochs,iXYZ);
  }
  //read/write reference clock position:
  for(size_t j=0; j<refpos.size(); j++){
    for(size_t k=0; k<refpos[j].size(); k++){
      binary_rw(iof,refpos[j][k],write);
    }
  }
  //read/write clock positions
  for(size_t i=0; i<pos.size(); i++){
    for(size_t j=0; j<pos[i].size(); j++){
      for(size_t k=0; k<pos[i][j].size(); k++){
        binary_rw(iof,pos[i][j][k],write);
      }
    }
  }
  iof.close();
  return 0;
}

//******************************************************************************
int JplGpsData::makeACopyOf(JplGpsData inObject)
/*
170830.
Makes a complete identical copy of an input JplGpsData object.
*/
{
  num_receivers=inObject.num_receivers;
  num_satellites=inObject.num_satellites;
  num_clocks=inObject.num_clocks;
  num_epochs=inObject.num_epochs;

  jpl_file_dir=inObject.jpl_file_dir;

  date=inObject.date;
  week=inObject.week;
  day=inObject.day;

  verbose=inObject.verbose;

  //Copy accross clock bias/ferr data, and clock names
  if(inObject.bias.size()!=0){
    reSizeVec(bias,inObject.bias.size(),inObject.bias[0].size());
    reSizeVec(ferr,inObject.ferr.size(),inObject.ferr[0].size());
    // If bias.size()=0, bias[0] will cause seg-fault
  }
  prn.clear();//in case already exists!
  svn.clear();
  clk.clear();
  blk.clear();
  for(std::size_t i=0; i<bias.size(); i++){
    for(std::size_t j=0; j<bias[0].size(); j++){
      bias[i][j]=inObject.bias[i][j];
      ferr[i][j]=inObject.ferr[i][j];
    }
    prn.push_back(inObject.prn[i]);
    svn.push_back(inObject.svn[i]);
    clk.push_back(inObject.clk[i]);
    blk.push_back(inObject.blk[i]);
  }
  refprn=inObject.refprn;
  refsvn=inObject.refsvn;
  refclk=inObject.refclk;
  refblk=inObject.refblk;
  refname=inObject.refname;

  //Copy across clock positions (if they exist):
  if(inObject.pos.size()!=0){
    reSizeVec(pos,inObject.pos.size(),inObject.pos[0].size(),
              inObject.pos[0][0].size());
  }
  for(std::size_t i=0; i<pos.size(); i++){
    for(std::size_t j=0; j<pos[0].size(); j++){
      for(std::size_t k=0; k<pos[0][0].size(); k++){
        pos[i][j][k]=inObject.pos[i][j][k];
      }
    }
  }

  //Copy across reference clock positions (if they exist):
  if(inObject.refpos.size()!=0){
    reSizeVec(refpos,inObject.refpos.size(),inObject.refpos[0].size());
  }
  for(std::size_t j=0; j<refpos.size(); j++){
    for(std::size_t k=0; k<refpos[0].size(); k++){
        refpos[j][k]=inObject.refpos[j][k];
    }
  }

  //Copy R^2 values (if they exist):
  dR2.clear();
  for(std::size_t i=0; i<inObject.dR2.size(); i++){
    dR2.push_back(inObject.dR2[i]);
  }

  //Copy std dev values (if they exist):
  sdev.clear();
  for(std::size_t i=0; i<inObject.sdev.size(); i++){
    sdev.push_back(inObject.sdev[i]);
  }

  //Copy vdiff values (if they exist):
  vdiff.clear();
  for(std::size_t i=0; i<inObject.vdiff.size(); i++){
    vdiff.push_back(inObject.vdiff[i]);
  }

  //Copy acf values (if they exist):
  if(inObject.acf.size()!=0){
    reSizeVec(acf,inObject.acf.size(),inObject.acf[0].size());
  }
  for(std::size_t i=0; i<inObject.acf.size(); i++){
    for(std::size_t j=0; j<inObject.acf[0].size(); j++){
      acf[i][j]=inObject.acf[i][j];
    }
  }

  //Copy cross-clock-correlation values (if they exist):
  b0.clear();
  for(std::size_t i=0; i<inObject.b0.size(); i++)
    b0.push_back(inObject.b0[i]);

  //The "effective" (relative) couplings between clocks
  keff.clear();
  for(std::size_t i=0; i<inObject.keff.size(); i++){
    keff.push_back(inObject.keff[i]);
  }

  return 0;
}


//******************************************************************************
//******************************************************************************
//******************************************************************************
int JplGpsData::gpsSimulator(
  std::string pathToPSD, //Path to the PSD files
  std::string inLabel,   //Label for the PSD files
  std::vector<int> &numClks,   //How many of each clock to simulate
  double sig,       //Standard deviation for the white clocks.
  int in_receivers, //number of recievers (base stations) to simulate
  double sta_sd,    //s.d. of recievers/base-stations
  std::string sRef,      //which ref clock to simulate
  int in_epochs,
  bool useRandSVN,  //Use randomised SVNs? (or use 'av')
  bool swap_white_ref,
  double ref_sig, // for simulating ccc
  bool print       //Print? Or be quiet
)
/*
170726.

The GPS simulator.
Generates a specified number of random clock solutions, and satelite positions.
The clock noise can either be white or can simulate a real clock type.
The real satelite clocks are simulated by "colouring" white noise with the known
GPS satelite clock Power Spectral Desnity (PSD) - as in Alex Rollings' thesis.
How many of each type of clock is an input (numClks).
You can use either 100% random satelite positions, or can semi-uniformly
place them (with just some random element); this more accurately portrays actual
satellite distribution.

The program optionally makes use of power spectrums (for simulating clocks) and
the inverse covariance functions (for covariance matrix). These have been
calculated already, and can be downloaded from github (see main README for
more detailed git instructions)
  * https://github.com/benroberts999/PowerSpectrums

e.g.:
  * $ wget http://github.com/benroberts999/PowerSpectrums/archive/master.tar.gz
  * $ tar xf master.tar.gz -C ./ --strip-components=1

This program calls the functions:
simulatedData, whiteData, and simulatedPositions

Uses GSL to perform FFT (in simulatedData)

====== CHANGE LOG ======
170829- Works with class
170912- Can now generate fake data with less than 2880 epochs! {not fully tested
170921- Uses GSL to perform FFT. Now works with <2880 points too.
171004- Added optional routine to "swap" reference clock to a white noise clock.
        the idea is that this will help mimic cross-clock correlations!
171012- Add abilty to simulate base stations (just white for now)
*/
{
  //count the number of clocks (= sum of each type):
  num_satellites=0;
  for(std::size_t iType=0; iType<numClks.size(); iType++)
    num_satellites+=numClks[iType];
  num_receivers=in_receivers;
  num_clocks = num_satellites + num_receivers;
  date="2000-00-00"; //?

  //writes Ref Clock info:
  refprn=sRef;
  refsvn="00";
  refname=sRef;
  refblk="AR"; // Means "reciever", i.e. eaerth-based clock
  refclk="H";  // The reference clocks for JPL files always H-masers (??)

  if(in_epochs>JPLEPOCHS){
    std::cout<<"FAILURE 1594 in JplGpsData::gpsSimulator. Cannot simulate data "
             <<"with "<<in_epochs<<" epochs (at the moment), due to the length "
             <<"of the input PSDs. 2880 is maximum.\n";
    return 2;
  }

  // dif_order is which difference level PSD to use to simulate data.
  int dif_order=1;
  // temp_epochs is the "actual" length of the generated time series.
  // Note: Only store the first num_epochs points of this time series!
  int temp_epochs = JPLEPOCHS - dif_order;

  //set the number of epochs (=2880 by default)
  num_epochs = in_epochs;

  //Now we know the number of clocks/epochs, we can "re-size" the vector arrays:
  reSizeVec(bias,num_clocks,num_epochs);  //clock bias (2d vector array)
  reSizeVec(ferr,num_clocks,num_epochs);  //clock formal error
  reSizeVec(prn,num_clocks); //prn ("G01" etc)
  reSizeVec(svn,num_clocks); //SVN (space vehicle number)
  reSizeVec(clk,num_clocks); //Clock type (Cs, Rb, H)
  reSizeVec(blk,num_clocks); //Block: II, IIA, IIR, IIF, or AR (for recievers)

  //R2. Don't need it, but needs to exist (?)
  reSizeVec(dR2,num_clocks);
  for(int i=0; i<num_clocks; i++){
    dR2[i]=1;
  }

  //minimum and maximum SVNs for each block:
  double IImin=15,IImax=17;
  double IIAmin=22,IIAmax=40;
  double IIRmin=41,IIRmax=61;
  double IIFmin=62,IIFmax=73;//Cs only has 63.. may cause problems?

  //names of each clock/block
  //Note: white clocks always last!
  std::string ClkBlk[8];
    ClkBlk[0]="RbIIF";
    ClkBlk[1]="RbIIA";
    ClkBlk[2]="RbIIR";
    ClkBlk[3]="RbII";
    ClkBlk[4]="CsIIF";
    ClkBlk[5]="CsIIA";
    ClkBlk[6]="CsII";
    ClkBlk[7]="white";
    const int IW=7; //index for the "white" clocks.

  //number of non-white clocks (for progress bar)
  int tempclocks=num_clocks-numClks[IW]; //for progress bar

  //Prepare some things for the GSL
  //Prepare trigonometric lookup tables for an FFT of size n real elements:
  //Note: We only need to form this once (cen be re-used)
  gsl_fft_real_wavetable * real_wt;
  real_wt = gsl_fft_real_wavetable_alloc (temp_epochs);
  //Prepare trigonometric lookup tables for FFT; size n half-complex elements:
  // (this if for the inverse Fourier transform)
  gsl_fft_halfcomplex_wavetable * hc_wt;
  hc_wt = gsl_fft_halfcomplex_wavetable_alloc (temp_epochs);
  //Prepare workspace for GSL Fourier transforms
  gsl_fft_real_workspace * work;
  work = gsl_fft_real_workspace_alloc (temp_epochs);

  int i=0; //clock index! Used in the three loops

  //Generate the white noise base STATION / reciever clocks:
  // Just white for now!
  for(int iSta=0; iSta<num_receivers; iSta++){
    prn[i]="G"+MSC_padIntString(i);
    clk[i]="Wh";
    blk[i]="AR"; //Station!
    svn[i]="00"; //MSC_padIntString(i);
    //generate the white data
    std::vector<double> tempData; //re-sized inside
    whiteData(tempData,sta_sd);
    for(int j=0;j<num_epochs;j++){
      bias[i][j]=tempData[j];
      ferr[i][j]=MFS_randDouble(0.009,0.011); //??
    }
    i++;//increment clock index
  }//loop over number of white STATION clocks

  //loop through each 'type' of non-white clock:
  //"-1" because not white:
  for(std::size_t iType=0; iType<numClks.size()-1; iType++){
    //loop through number of clocks of each type:
    for(int iTsvn=0; iTsvn<numClks[iType]; iTsvn++){
      if(print)MSC_progressBar("Simulating Clock Data ",35,i,tempclocks);

      //check to see if this clock type is OK:
      std::string psdName=pathToPSD+ClkBlk[iType]+"-av-"+sRef+"-"
                          +inLabel+".psd";
      if(!MSC_fexists(psdName)){
        //There were NO clocks of this type (in the PSD files)
        std::cout<<"\n\nFile: "<<psdName
                 <<" DNE. Clocks not used"<<" ("<<numClks[iType]<<")\n"
                 <<" -> "<<numClks[iType]<<" extra white clks used instead\n\n";
        numClks[IW]+=numClks[iType]; //increase number of white clocks
        numClks[iType]=0; //kill these "bad" non-white clocks
        continue;
        return 2;
      }

      //Fill clock name array:
      prn[i]="G"+MSC_padIntString(i);
      clk[i]=ClkBlk[iType].substr(0,2);//take first 2 characters
      blk[i]=ClkBlk[iType].substr(2);//ignore first two characters
      svn[i]=MSC_padIntString(i);
      //svn may be re-written below:

      //Choose SVN (either a random SVN, or the average):
      if(useRandSVN){
        bool foundSVN=false;
        int iTries=0;
        while(!foundSVN){
          //randomly pick SVNs from correct range until a good one is found
          int iSVN=0;
          if(iType==0)iSVN=int(MFS_randDouble(IIFmin,IIFmax+1));
          if(iType==4)iSVN=65; //only CsIIF SVN..
          if(iType==1||iType==5)iSVN=int(MFS_randDouble(IIAmin,IIAmax+1));
          if(iType==2)          iSVN=int(MFS_randDouble(IIRmin,IIRmax+1));
          if(iType==3||iType==6)iSVN=int(MFS_randDouble(IImin,IImax+1));
          psdName=pathToPSD+ClkBlk[iType]+"-"+std::to_string(iSVN)+"-"+sRef+"-"
                  +inLabel+".psd";
          svn[i]=MSC_padIntString(iSVN);
          if(MSC_fexists(psdName))foundSVN=true;
          iTries++;
          if(iTries>250){
            //If it coulen't find a reasonable SVN, just use 'av'
            //Hopefully this never actually happens, just a safe-guard
            psdName=pathToPSD+ClkBlk[iType]+"-av-"+sRef+"-"+inLabel+".psd";
            svn[i]="00"; //MSC_padIntString(i);
            if(MSC_fexists(psdName)){
              foundSVN=true;
            }else{
              //If still couldn't find, use Hmas instead of sRef
              psdName=pathToPSD+ClkBlk[iType]+"-av-Hmas-"+inLabel+".psd";
              if(MSC_fexists(psdName)) foundSVN=true;
            }
          }// if tries > 250
        }//END while
      }//END if randSVN

      //Open the PSD file, read it in to array:
      std::vector<double> power_spectrum; //no size yet, will be "pushed back"
      std::ifstream psdfile;
      psdfile.open(psdName.c_str());
      std::string sLine;
      while(getline(psdfile,sLine)){
        //Look for "comment" lines (marked with '#'):
        std::stringstream ssCheck(sLine);
        std::string checkLine;
        ssCheck>>checkLine;
        if(checkLine=="#")continue; //skip comment lines
        //read in PSD:
        std::stringstream ssin(sLine);
        double sPSD,sPSD0,sPSD1,sPSD2;
        ssin>>sPSD0>>sPSD1>>sPSD2;
        //Choose to use 0th or 1st order differenced data:
        if(dif_order==0) sPSD=sPSD0;
        else sPSD=sPSD1;
        //don't read in the trailing zero's
        if(sPSD!=0) power_spectrum.push_back(sPSD);
      }

      // To store the output time-series:
      std::vector<double>temp_data;

      //Generate simulated GPS time series:
      //Use Alex's method to generate simulated data that matches the
      //given PSD:
      simulatedData(power_spectrum,temp_data,temp_epochs,
                    dif_order,real_wt,hc_wt,work);
      for(int j=0;j<num_epochs;j++){
        bias[i][j]=temp_data[j];
        ferr[i][j]=MFS_randDouble(0.02,0.03); //?? Random formal error
      }

      i++;//increment clock index
    }//loop over clocks of that type (SVNs)
  }//END loop over non-white clock types

  //Clear the memory ascociated with Trig tables and workspace
  gsl_fft_real_workspace_free (work);
  gsl_fft_real_wavetable_free (real_wt);
  gsl_fft_halfcomplex_wavetable_free (hc_wt);

  //loop through and generate the white noise SATtelite clocks:
  for(int iTsvn=0; iTsvn<numClks[IW]; iTsvn++){
    prn[i]="G"+MSC_padIntString(i);
    clk[i]="Wh";
    blk[i]="GAU";
    svn[i]=MSC_padIntString(i);
    //generate the white data
    std::vector<double> tempData; //re-sized inside
    whiteData(tempData,sig);
    for(int j=0;j<num_epochs;j++){
      bias[i][j]=tempData[j];
      ferr[i][j]=MFS_randDouble(0.02,0.03); //??
    }
    i++;//increment clock index
  }//loop over number of white clocks
  if(print&&(i!=numClks[IW])) MSC_progressBar("Simulating Clock Data ",35);

  //re-count number of clocks (check for failure):
  int double_check_clocks=num_receivers;
  for(int iType=0;iType<8;iType++) double_check_clocks+=numClks[iType];
  if(i!=double_check_clocks||i!=num_clocks||num_clocks!=double_check_clocks){
    std::cout<<"\nFAILURE: 1736 in gpsSimulator: "
             <<i<<" "<<num_clocks<<" "<<double_check_clocks
             <<"\n\n"<<std::flush;
    return 2;
  }

  //Simualte cross-clock-correlations, by swapping to a white-noise ref. clock!
  if(swap_white_ref){
    //Routine to swap the reference clock
    std::vector<double> white_ref(num_epochs);
    //Generate un-differenced reference clock data:
    white_ref[0] = MFS_randGausVal(0,ref_sig);
    for(int j=1; j<num_epochs; j++)
      white_ref[j] += MFS_randGausVal(0,ref_sig);
    //Subtract from each actual clock:
    for(int i2=0; i2<num_clocks; i2++){
      //Know the (approximate) original s.d., to "normalise"
      double tempsd=sig; //for white clocks
      if(clk[i2]=="Rb"&&blk[i2]=="II")        tempsd=0.048;
      else if(clk[i2]=="Rb"&&blk[i2]=="IIA")  tempsd=0.040;
      else if(clk[i2]=="Rb"&&blk[i2]=="IIR")  tempsd=0.074;
      else if(clk[i2]=="Rb"&&blk[i2]=="IIF")  tempsd=0.013;
      else if(clk[i2]=="Cs"&&blk[i2]=="II")   tempsd=0.083;
      else if(clk[i2]=="Cs"&&blk[i2]=="IIA")  tempsd=0.089;
      else if(clk[i2]=="Cs"&&blk[i2]=="IIF")  tempsd=0.098;
      else if(blk[i2]=="AR")                  tempsd=sta_sd;
      double norm_const = tempsd/sqrt(pow(tempsd,2) + pow(ref_sig,2) );
      //Subtract the new "reference" clock!
      for(int j=0; j<num_epochs; j++){
        bias[i2][j] = (bias[i2][j] - white_ref[j])*norm_const;
      }
    }
  }

  //re-size the position vectors. Note: entrie vector!
  reSizeVec(pos,num_clocks,num_epochs,iXYZ);
  reSizeVec(refpos,num_epochs,iXYZ);

  //Simulate the satelite and reference clock positions:
  //Semi-random, or completely random sat positions?
  bool bFull_Random=false; //false=use "semi" evenly distributed
  simulatedPositions(bFull_Random);

  return 0;
}


//******************************************************************************
int JplGpsData::simulatedData(std::vector<double> power_spectrum,
                              std::vector<double> &out_z0, int temp_epochs,
                              int dif_order,
                              gsl_fft_real_wavetable * real_wt,
                              gsl_fft_halfcomplex_wavetable * hc_wt,
                              gsl_fft_real_workspace * work
                              )
/*
170323.
Generates simulated bias clock data from a given PSD for a single clock.
Note: it generates data (d^(0)) that does NOT include a 2nd order poly.
Therefore, a polynomial should NOT be subtracted from this data!
Uses method the from Alex Rollings.

    y(t) = random (phase data)
    Py(k) = Power spectrum of y
    Pt(k) = "Target" Power spectrum [of real data]
    z(t) = Simulated time series
    ============================
    Method:
    Generate y(t)
    y*(k) = FT[y(t)]
    z*(k) := y*(k) * Sqrt[ P(k)/Py(k) ]
    z(t) = inverseFT[ z*(k) ]

Takes in a single PSD, inPSD,
outputs a single un-differenced time series, outZ0.
NOTE: At the moment, I use the power spectrum for the first-order differenced
clock solutions. This seems to work more consistently. Therefore the simulated
clock solutions will also be first-order differenced solutions, and have to be
'undifferenced' to bring them into the desired form (so the rest of the program
can run as per the real data).
Note: Uses the GSL Fast Fourier libraries to perorm FFT.
See:
https://www.gnu.org/software/gsl/manual/html_node/Fast-Fourier-Transforms.html#Fast-Fourier-Transforms
https://www.gnu.org/software/gsl/manual/html_node/Mixed_002dradix-FFT-routines-for-real-data.html

*It actually generates a time-series that is temp_epochs long
[Either 2880 for non-differenced, or 2879 if 1-differenced].
Then, just keeps the first num_epochs points.

INPUT:
  power_spectrum  :: input (target) power spectrum
  temp_epochs     :: 2880 or 2879, see above
  dif_order       :: 0 or 1 {always 1 for now}
  real_wt         :: Trig lookup-table for FT [GSL]
  hc_wt           :: "" "" for inverse FT [GSL]
  work            :: FT workspace [GSL]
OUTPUT:
  out_Z0 ::  Output array [num_epochs] of UN-differenced simulated time series


// nb: for 0-order difference, temp_epochs = JPLEPOCHS
// For 1-order, temp_epochs = JPLEPOCHS -1
====== CHANGE LOG ======
170829- Works with class
170912- Fixed the way it works for J<2880. Not 100% efficient, but works(?)
*/
{
  //Length of input power spectrum.
  // When using 1st order diff'd data, this is 1 less than temp_epochs
  int n_psd = power_spectrum.size();

  //generate random "phase" data:
  double y_data[JPLEPOCHS];
  for(int j=0;j<temp_epochs;j++){
    y_data[j]=MFS_randDouble(-1.,1.);
  }

  //Fourier transform random "phase" data:
  int stride = 1;
  gsl_fft_real_transform (y_data, stride, temp_epochs, real_wt, work);

  //Scale by psd
  //Multiply the frequency domain 'phase data' from above with the
  //desired "target" power spectrum, while dividing by the power spectrum of
  //the 'phase data'.
  double z_data[JPLEPOCHS];
  z_data[0] = y_data[0] * sqrt( power_spectrum[0] / pow(y_data[0],2) );
  for (int k=1; k<n_psd; k++){
    double re_y = y_data[2*k-1];
    double im_y = y_data[2*k];
    double Rk = sqrt( power_spectrum[k] / (re_y*re_y + im_y*im_y) );
    z_data[2*k-1] = re_y * Rk;
    z_data[2*k]   = im_y * Rk;
  }

  //Generate the "differenced" time series: z^(1), via inverse FT
  gsl_fft_halfcomplex_inverse (z_data, stride, temp_epochs, hc_wt, work);

  //UN-difference: z^(1) -> z^(0)
  //Note: move from temp_epochs to num_epochs
  reSizeVec(out_z0,num_epochs);
  if(dif_order==1){
    #pragma omp parallel for
    for (int j=0;j<num_epochs;j++){
      out_z0[j]=0; //not needed, but safer(?)
      for(int k=0;k<j;k++){//integrate white noise = random walk
        out_z0[j]+=z_data[k];
      }
    }
  }else if(dif_order==0){
    for (int j=0;j<num_epochs;j++){
      out_z0[j]=z_data[j];
    }
  }

  return 0;
}


//******************************************************************************
int JplGpsData::whiteData(std::vector<double> &out_z0, double sig)
/*
170726.
Generates integrated Gaussian white noise bias clock data for a single clock.
Integrated: turns white-noise into random walk. This is so this data can be
differenced, as per usual.
Note: it generates data (d^(0)) that does NOT include a 2nd order poly.
Therefore, a polynomial should NOT be subtracted from this data!

INPUT:
  sig     ::  Target standard deviation
OUTPUT:
  outZ0   ::  Output array of UN-differenced Gaussian white noise time series

====== CHANGE LOG ======
170829- Works with class

*/
{
  //Generate Gaussian white noise
  std::vector<double> z1(num_epochs);
  for (int j=0;j<num_epochs;j++){
    z1[j]=MFS_randGausVal(0,sig);
  }

  //UN-difference: z^(1) -> z^(0)
  //i.e. take "white" noise -> random walk!
  reSizeVec(out_z0,num_epochs);
  #pragma omp parallel for
  for (int j=0;j<num_epochs;j++){
    out_z0[j]=0; //not needed, but safer(?)
    for(int k=0;k<j+1;k++){//integrate white noise = random walk
      out_z0[j]+=z1[k];
    }
  }

  return 0;
}

//******************************************************************************
int JplGpsData::simulatedPositions(bool bFull_Random)
/*
170726 (written earlier, now moved into 'GPSsimulator').
Generates random or semi-random positions for the satelites and reference clock.
Alternatively, could read in an ECI file, and use actual satellite pos
NOTE: only generates these for a single epoch!!

Uses inverse transform sampling, so that the satellitess don't become "bunched"
up at the poles (i.e. takes solid angle into account).

bFull_Random        :: Semi-randomly distribute sats (=0), or 100% random (=1)
                       --> optional. False by default.

Note: it generates the entire position vector, filled with the same positions.
this is just to make it as compatible as possible with actual data

===== Change Log =====
171012- Added ability to place base-stations on Earth-level
*/
{
  double dp=1./(num_clocks+1); //step-size for semi-evenly distributing sats.
  double th_er=0.1; //s.d. for "uncertainty" in angles (in radians)

  //Generate the positions for the satellite clocks:
  double r,theta,phi,u;
  for(int i=0;i<num_clocks;i++){//loop over each clock

    //generate spherical polar coordinates (r, theta, phi)
    if(i<num_receivers)
      r=MFS_randGausVal(REARTH,0.05*REARTH);  //Earth (station) radius
    else
      r=MFS_randGausVal(RGPS,0.02*RGPS);  //Sat. radius

    //choose angles:
    if(bFull_Random){
      //distribute the satellites completely randomly
      phi  =MFS_randDouble(0,2*PI);
      u    =MFS_randDouble(0,1);
      theta=acos(1.-2.*u);
    }else{
      // semi-evenly distribute the sats (in a spiral)
      phi  =(i*dp)*10.*2*PI + MFS_randGausVal(0.,th_er);
      theta=acos(1.-2.*(i+1)*dp)+ MFS_randGausVal(0.,th_er);
    }

    //calculated cartesian (xyz) positions from sph polar (r,t,p):
    for(int j=0; j<num_epochs; j++){
      //just put same numbers in each epoch. Bit dumb...
      pos[i][j][XX]=r*sin(theta)*cos(phi);  //x position
      pos[i][j][YY]=r*sin(theta)*sin(phi);  //y position
      pos[i][j][ZZ]=r*cos(theta);           //z position
    }

  }//end loop over clocks

  //Generate position of the ref. clock (positioned on the Earth)
  r    =MFS_randGausVal(REARTH,0.05*REARTH);  //Earth radius
  phi  =MFS_randDouble(0,2*PI);
  u    =MFS_randDouble(0,1);
  theta=acos(1.-2.*u);
  for(int j=0; j<num_epochs; j++){
    //just put same numbers in each epoch. Bit dumb...
    refpos[j][XX]=r*sin(theta)*cos(phi);  //x position
    refpos[j][YY]=r*sin(theta)*sin(phi);  //y position
    refpos[j][ZZ]=r*cos(theta);           //z position
  }

  return 0;
}




// //////////////***************************************************************
int JplGpsData::readJplBiasData_1s(std::string path, int in_day)
/*
Code written by Guglielmo.
Reads in the pre-sorted 1s GPS files.
The files are sorted by SVN.
NOTE: Doens't do everything yet, not 100% compatible with rest of code.
e.g., doesn't lookup clock/prn etc.
*/
{
	//Sets GW event name and day
	//std::string event=in_event;
	std::string day=std::to_string(in_day);

	//Store location of the data files
	std::string jpl_file_dir=path;

	//Opens the input file
	std::ifstream inFile(jpl_file_dir);
	if(!inFile){
		std::cout<< "Cannot open file.\n";
	}


	//Count number of data lines in the file
	int iLines=0;
	std::string buffLine;
	while(std::getline(inFile,buffLine)){
		iLines++;
	}
	inFile.close(); //Must re-open

	std::ifstream inFile1(jpl_file_dir);
	//Loads each line of the file into a string array
	std::vector<std::string> fileLine(iLines);
	std::string sLine;
	int iLine=0;

	while(std::getline(inFile1,sLine)){
		if(sLine.size()>0)
		fileLine[iLine]=sLine;
		//std::cout<<fileLine[iLine]<<std::endl;
		iLine++;
	}
	inFile1.close(); //Don't need file anymore

	//Variables that will store file information
	// char rs; //Receiver or Sat clock
  std::string rs;
	std::string clockName;
	int wk, dy, gpsec;
	long double dbias, sigma;
	double nano=pow(10.,9.);

	//Set number of epochs
	std::stringstream ssin(fileLine[0]);
	ssin>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
	int begin_epoch = gpsec;

	std::stringstream ssin1(fileLine[iLines-1]);
	ssin1>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
	int end_epoch = gpsec;

	//int num_epochs=end_epoch-begin_epoch+1;
  // num_epochs is a member of the class
  num_epochs=end_epoch-begin_epoch+1;

	// std::cout<<"The first epoch is: "<<begin_epoch<< std::endl;
	// std::cout<<"The last epoch is: "<<end_epoch<< std::endl;
	// std::cout<<"The total number of epochs is: "<<num_epochs<< std::endl;

	//Count number of clocks
	std::string tempClock1=" ";
	int temp_num_clocks=0;
	for(int i=0; i<iLines; i++){
		std::stringstream tempin(fileLine[i]);
		tempin>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
		if(tempClock1 != clockName){
			tempClock1 = clockName;
			temp_num_clocks++;
		}
	}
  num_clocks = temp_num_clocks; // num_clocks is a class variable
	//std::cout<<"The number of clocks is: "<<num_clocks<<std::endl;

  //XXX Ben:
  reSizeVec(bias,num_clocks,num_epochs);  //clock bias (2d vector array)
  reSizeVec(ferr,num_clocks,num_epochs);  //clock formal error
  reSizeVec(prn,num_clocks); //prn ("G01" etc)
  reSizeVec(svn,num_clocks); //SVN (space vehicle number)
  reSizeVec(clk,num_clocks); //Clock type (Cs, Rb, H)
  reSizeVec(blk,num_clocks); //Block: II, IIA, IIR, IIF, or AR (for recievers)

	//Store clock names and clock types
	// std::string clockNames[num_clocks];
	// std::string clockType[num_clocks];
	std::string tempClock2=" ";
	int count=0;
	for(int i=0; i<iLines; i++){
		std::stringstream tempin1(fileLine[i]);
		tempin1>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
		if(tempClock2 != clockName){
			tempClock2 = clockName;
			// clockNames[count]=tempClock2; //
			// clockType[count]=rs;
      if(rs=="R"){ //XXX I am trying to check if it's a reciever. Right?
        prn[count]=tempClock2;
        svn[count]="00"; //SVN for stations is 0
        blk[count]="AR"; //'AR' is what I call 'block' for stations
        clk[count]="**"; //we don't know yet which is which clock
      }else{
        prn[count]="**";
        svn[count]=tempClock2.substr(2);// substr removes the 'GP' part
        blk[count]="**";
        clk[count]="**"; //we don't know yet which is which clock
      }
			count++;
		}
	}

/*
	std::cout<<"\nClock  TYPE  NAME "<<std::endl;
	for(int i=0; i<num_clocks; i++){
		std::cout<<"Clock["<<i+1<<"]: "<<clockType[i]<<"   "<<clockNames[i]<<std::endl;
	}
*/

	//Size vector arrays for bias and ferr
	//num_epochs=1; //XXX With real file erase this
	// double bias[num_clocks][num_epochs];
	// float ferr[num_clocks][num_epochs];
	//Set values of bias ferr arrays
	for(int i = 0; i<num_clocks; i++){
		for(int j=0; j<num_epochs; j++){
			//std::stringstream tempin2(fileLine[i]);
      std::stringstream tempin2(fileLine[j+i*num_epochs]);
			tempin2>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
      //convert from seconds to nanoseconds
			bias[i][j]=dbias*nano;
			ferr[i][j]=sigma*nano;
		}
	}
	//std::cout<<bias[1][1]<<"   "<<ferr[1][1]<<std::endl;

/*	//print out bias and ferr for each clock
	std::cout<<"\nClock:         BIAS           FERR   "<<std::endl;
	for(int i = 0; i<num_clocks; i++){
		for(int j=0; j<num_epochs; j++){
			std::cout<<"Clock["<<i+1<<"]:  "<<bias[i][j]<<"    "<<ferr[i][j]<<std::endl;
		}
	}
*/

	//reSizeVec(bias,num_clocks,num_epochs);
	//reSizeVec(ferr,num_clocks,num_epochs);

	return 0;
}
