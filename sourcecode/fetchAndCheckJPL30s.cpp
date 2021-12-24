/*
PROGRAM:

Simple program that downloads all the JPL 30s files from
ftp://sideshow.jpl.nasa.gov/pub/jpligsac/
Uses wget, so that part will only work on linux machines.
Also, unzips the files.
For windows, probably best to do this manually using an ftp client.

Any files that are not found on the JPL website are recorded in a log:
"NotFoundJPL.txt"

Before running, it requires a 'y' key-press from the user.

Optionally, it also performs a few checks on the data, and outputs the
results to text files:

  * **PRN_GPS-SVN-faillog.txt**: Any of the SVN mappings in our modified
    PRN_GPS_GPSDM.txt file don't match those in the original PRN_GPS file
  * **ListOfRefs**: Outputs a list of all used reference clocks
  * **NotFoundJPL.txt**: Any day for which there wasn't a corresponding data
    file on the JPL website
  * **ECI-faillog.txt** (Actually produced inside JplGpsDataClass::readEciPos):
    Any day that something went wrong in the ECI position files. E.g.,
    if there were SVNs in the ECI files that weren't in the JPL files (this is
    not actually a problem), if there were SVNs in the JPL files that weren't in
    the ECI files (this is a problem) etc.
  * **PRN_GPS-faillog.txt** (Actually produced inside JplGpsDataClass::mapSVNPRN):
    Any PRN for which the clock/svn wasn't identified


=== Change Log ===
170831- works with class. Merged with 'PRNfail"

=== TO DO ===


*/
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
#include "JplGpsDataClass.h"


//******************************************************************************
//******************************************************************************
int main (void){

  //input parameters:
  std::string path;      //holds file location
  std::string bin_path;  // location for output binary files
  int minweek,maxweek;   //First and last week to consider
  int icheck;            // Check everything OK? or just get files.
  int ibin,ibin_pos;              //write to binary?

  //read the input file:
  std::string junk;
  std::ifstream fInput;
  fInput.open ("fetchAndCheckJPL30s.dat");
    fInput >> path;                        getline(fInput,junk);
    fInput >> minweek >> maxweek;          getline(fInput,junk);
    fInput >> icheck;                      getline(fInput,junk);
    fInput >> ibin >> ibin_pos;            getline(fInput,junk);
    fInput >> bin_path;                    getline(fInput,junk);
  fInput.close();

  //format input to reduce possibility of errors:
  if(minweek==0)minweek=1060; //- first week with 30s data!;
  if(maxweek==0)maxweek=5000;
  int minday=0,maxday=7;

  //Check that everything's ok in the files? (true), or just fetch them (false)
  bool check=false;
  if(icheck==1)check=true;

  bool write_binary=false;
  bool write_binary_pos=false;
  if(ibin==1){
    check=true; //must open the files to convert them!
    write_binary=true;
  }
  if(ibin_pos==1) write_binary_pos=true;

  //Output some info to screen
  std::cout<<"\n*******************************************\n\n";
  std::cout<<"Fetches (and unzips) JPL 30s CLK files from: \n"
       <<"   ftp://sideshow.jpl.nasa.gov/pub/jpligsac/ \n"
       <<"and will place them:\n"
       <<"   "<<path;
  std::cout<<"\n\n*******************************************\n\n";

  if(write_binary){
    std::cout<<"\n Also, will convert files to binary format, placing "
    <<"the outputs here: "<<bin_path<<"\n";
    if(write_binary_pos) std::cout<<"-> Including the positions!\n";
    std::cout<<"\n";
  }

  //To avoid mistakes, requires a "y" key entered by user:
  std::string continue_yn;
  std::cout << "Note: about to download lots of files. "
            <<"Double check the file locations above.\n";
  std::cout<<"\n NOTE: do NOT press cntrl+c to kill this progam, unless you "
    <<"manually check the downloaded files, otherwise only part of the file may"
    <<" be downloaded\n";
  std::cout<<"\n ** Continue? (y/n) \n";
  getline (std::cin, continue_yn);
  if(continue_yn!="y"&&continue_yn!="Y"&&continue_yn!="yes"){
    std::cout<<"Fine. I didn't really want to either. So there.\n";
    return 1;
  }

  //array to store each unique reference clock name!
  std::vector<std::string> listRefs;

  //file to write the "missing JPL files" to
  std::ofstream nojpl;
  std::string nojplName="NotFoundJPL.txt";
  nojpl.open(nojplName.c_str());

  int count=0; //# of non-found files in a row!

  //Geoff only made the ECI files for 1280--1881.
  int min_ECI_week=1280;
  int max_ECI_week=1881;

  //Loop through all possible days:
  for(int iw=minweek;iw<=maxweek;iw++){  //loops through weeks
    for(int jd=minday;jd<maxday;jd++){//loops through each day for given week
      if(iw==1060&&jd<5)continue;//skip non-30s days

      //create empyt data object
      JplGpsData data;

      //print out any messages:
      data.verbose=true;

      //check for (and if needed, download) the files:
      int iok;
      if(check){
        //if we're also 'checking' that everything's ok, need to read files
        iok=data.readJplBiasData(path,iw,jd,true); //true means download!
        std::cout<<"Reading: jpl"<<data.week+data.day<<" ("<<data.date<<") "
                 <<"Ref: "<<data.refname<<std::endl;
      }else{
        data.week = std::to_string(iw);
        data.day = std::to_string(jd);
        iok=data.checkJPLfiles(path,true);
      }

      if(iok==2){
        return 2; //broken. (coudn't create directory)
      }else if(iok==1){
        nojpl<<"jpl"<<iw<<jd<<"\n";
        count++; //add one more "not found in a row"
      }else if(iok==0){
        count=0;
      }

      if((count>15)&&(iw>1970)){
        //once there are 15 non-found JPL files in a row (after week 1970)
        //Program assumes that's because we have fetched all of them
        std::cout<<"Fetched all the available files!\n\n";
        break; //breaks out of 'day' loop only. see below
      }

      if(iok==1) continue;
      if(!check) continue;
      //the rest is only for "checking" all the files/ counting ref clocks etc.

      //find "new" reference clocks:
      bool newRef=true;
      for(size_t k=0; k<listRefs.size(); k++){
        if(data.refprn==listRefs[k]){
          newRef=false;
          break;
        }
      }
      if(newRef){
        listRefs.push_back(data.refprn);
      }

      //Check the SVN assignments:
      JplGpsData dataCopy;
      dataCopy.makeACopyOf(data);
      dataCopy.mapSVNPRN(path,false); //read in the "original" version

      //Compare our updated PRN_GPS file to original one.
      //Clock assignments in our file are more reliable, but
      //SVN's should be ccorrect in original file!
      for(int k=data.num_receivers;k<data.num_clocks;k++){
        if(std::stoi(data.svn[k])!=std::stoi(dataCopy.svn[k])){
          std::cout<<"Error 87: SVN assignments don't match! "<<data.svn[k]<<" "
              <<dataCopy.svn[k]<<"("<<data.prn[k]<<data.clk[k]<<")\n";
          std::ofstream oFile;
          std::string filename="PRN_GPS-SVN-faillog.txt";
          //open file with append
          oFile.open (filename.c_str(),std::ios_base::app);
          oFile<<data.date<<" "<<data.week<<data.day<<" PRN:"<<data.prn[k]
               <<" SVN (our file):"<<data.svn[k]<<" "<<data.clk[k]
               <<". True SVN:"<<dataCopy.svn[k]<<std::endl;
          oFile.close();
        }
      }

      //Check ECI files:
      if((iw>=min_ECI_week)&&(iw<=max_ECI_week)){
      //Geoff only made the ECI files for 1280--1881.
        data.readEciPos(); //write-out happens inside here
      }

      if(write_binary){
        bool write=true;
        if(data.pos.size()==0) write_binary_pos=false;
        data.binaryJpl30s(bin_path,iw,jd,write_binary_pos,write);
      }


    }//END for(int j=minday;j<maxday;j++)
    if((count>15)&&(iw>1970)) break; //break out of week loop
  }//END loop over weeks

  nojpl.close(); //close the "not found on JPL" file

  //write out the reference clocks:
  if(check){
    std::ofstream refFile;
    std::string refFileName="ListOfRefs.txt";
    refFile.open(refFileName.c_str());
    std::cout<<"\n I found "<<listRefs.size()<<" reference clocks:\n";
    for(size_t k=0; k<listRefs.size(); k++){
      std::cout<<listRefs[k]<<std::endl;
      refFile<<listRefs[k]<<std::endl;
    }
    refFile.close();
  }

  return 0;
}
//******************************************************************************
//******************************************************************************
