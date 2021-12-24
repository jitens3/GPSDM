#include "miscFunctions.h"
/*
170827.
Includes "miscelanious" functions.
All functions start with MSC_ designator, 
which means they are found in this file.
*/


//******************************************************************************
bool MSC_fexists(const std::string& filename)
/*
Checks if a file exists, returns true or false.
This might actually be quite slow? (because it has to open the file each time?)
Not really tested, but it's nice to be 'safe' even if it adds a (tiny amount)
of overhead.
*/
{
  std::ifstream infile(filename.c_str());
  return infile.good();
}

//******************************************************************************
bool MSC_direxists(const std::string& DIR)
/*
Checks if a directory exists, returns true or false
Works by trying to creat a 'junk' file. Then, checks if that file exists.
If that file exists, it means the directory exists, return true.
*/
{
 std::string randname=DIR+"/fae9r78q4vynvpw3498qvynw9p3avby.tmp";
 std::ofstream testfile;
 bool bDir;
 if(!MSC_fexists(randname)){
   //don't over-ride on the off-chance such file already exists!
   testfile.open (randname.c_str());
   testfile.close();
   bDir=MSC_fexists(randname);
   MSC_execute("rm -f "+randname);//only linux, but that's ok
 }else{
   bDir=true;
 }
 return bDir;
}

//******************************************************************************
int MSC_execute(std::string exe)
/*
Takes in a string, and executes a (linux) system command!
(Caution: can be dangerous!)
Will only work on Linux machines.
*/
{
  int ret=system(exe.c_str());
  return ret;
}


//******************************************************************************
void MSC_flag(int s1, bool print)
//Flag, for debugging.
// Bool is true by default, so can be called without that option.
//However, that option makes it easy to turn on/off the flags
{
 std::string ss=std::to_string(s1);
 if(print) std::cout<<"#:"+ss<<"\n";
}

//******************************************************************************
void MSC_progressBar(std::string text, int barWidth, 
                    double current, double total, bool OMP)
/*
160913.
Draws a progress bar, with %'ge done to the screen.
works with openmp.
Note: must be called again after the loop with MSC_progressBar(text,barWidth)
to 'finish' the bar (set to 100%, and end the line)
INPUT:
  text      :: Text to display in front of progress bar. Can be none.
  barWidth  :: Integer. How wide to make the progress bar.
  current   :: Integer->double. Current position in loop
  total     :: Integer->double. Total number iterations through the loop
  OMP       :: Optional. If this particular loop is parallisised. (True/False)
               --> optional. False by default

*/
{
  //total=total-1; //??
  int numThreads=1;
  #if defined(_OPENMP)
  if(OMP)numThreads=omp_get_num_threads();
  #endif
  double progress;
  if(total!=0){
    progress=current/total;
  }else{
    progress=0;
  }
  //if(progress==1)MSC_progressBar(text,barWidth);
  if(progress*numThreads<1||progress==0){
    std::cout <<" "<<text<<" [";
    progress*=numThreads;
    //if(progress==0)progress=1; //XXX ???
    if(total==0)progress=1; //new?
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    {
      std::cout.precision(3);
      std::cout << "] " << (progress* 100.0) << " %   \r";
      std::cout.flush();
    }
  }
  if(current==0&&total==0)std::cout<<std::endl;
}




//******************************************************************************
std::string MSC_padIntString ( int in_int, int digits )
/*
Converts integer to padded string, with given number of padded digits.
e.g., if digits=3, then '1' -> '001'
Doesn't shorten the string, so if digits=3, then '1234' -> '1234'
'digits' is an optional parameter, by default, it is 2
*/
{
  
  int abs_int=abs(in_int);
  std::string out_string=std::to_string(abs_int);
  
  //Append the correct number of 0s to front of string:
  for(int i=digits-1; i>0; i--){
    if(pow(10,i)>abs_int){
      out_string="0"+out_string;
    }else{
      break;
    }
  }
  
  //If negative, append the "-":
  if(in_int<0) out_string="-"+out_string;
  
  return out_string;
}








//******************************************************************************
int MSC_ymdtowd(int yy, int mm, int dd)
/*
170320
Converts a date (yyyy,mm,dd) to a JPL week/day wwwwd (string).
 Uses "difftime" to work out the difference in days/weeks between
the given date and 1980-Jan-06 (week 0 day 0 in JPL time)
NB: year indexed from 1900. Month goes 0-11 (?)
http://www.cplusplus.com/reference/ctime/difftime/

Because of overloaded functions, input can be of format:
  int yy, int mm, int dd
  string yyyymmdd
*/
{
 struct tm jpl0,inDate;
 int yb=1980, mb=1, db=6;
 int jplweek=-1,jplday=-1;
 jpl0.tm_hour = 0; jpl0.tm_min = 0; jpl0.tm_sec = 0;
 jpl0.tm_year = yb-1900;  jpl0.tm_mon = mb-1;  jpl0.tm_mday = db;
 inDate.tm_hour = 0; inDate.tm_min = 0; inDate.tm_sec = 0;
 inDate.tm_year = yy-1900;  inDate.tm_mon = mm-1;  inDate.tm_mday = dd; 
 double seconds=difftime(mktime(&inDate),mktime(&jpl0));
 int days=(int)round(seconds/(60*60*24));
 jplweek = int(days/7);
 jplday = days-jplweek*7;
 
 return (10*jplweek+jplday);
}

//--Overloaded---------------------------------
int MSC_ymdtowd(std::string yyyymmdd)
{

  std::stringstream ss(yyyymmdd);
  int yy,mm,dd;

  ss>>yy>>mm>>dd;
  
  if(yy>9999){
    int ymd = yy;
    yy=int(ymd/10000);
    mm=int((ymd-yy*10000)/100);
    dd=int(ymd-yy*10000-mm*100);
  }

  //need the 'abs', since the "-" delimeter is read as a negative
  return MSC_ymdtowd(yy, abs(mm), abs(dd));

}




//******************************************************************************
std::string MSC_wdtoymd(int ww, int dd)
/*
170320
Converts a JPL "wwwwd" format date to "yyyy-mm-dd" format (string).
http://www.cplusplus.com/reference/ctime/
Uses ctime functions to convert JPLs "wwwwd" format to "yyyy-mm-dd" format
Because of overloaded functions, input can be of format:
  int wd
  int w, int d
  string wwwwd
*/
{
 struct tm jpl0;
 int yb=1980, mb=1, db=6;
 jpl0.tm_hour = 0; jpl0.tm_min = 0; jpl0.tm_sec = 0;
 jpl0.tm_year = yb-1900;  jpl0.tm_mon = mb-1;  jpl0.tm_mday = db;
 const time_t days=60*60*24*(7*ww+dd);
 time_t newdate = mktime(&jpl0)+days;
 int oYear=gmtime(&newdate)->tm_year+1900;
 int oMon=gmtime(&newdate)->tm_mon+1;
 int oDay=gmtime(&newdate)->tm_mday;
 
 std::string outDate=std::to_string(oYear)+"-"+MSC_padIntString(oMon)+"-"
                     +MSC_padIntString(oDay);
 return outDate;
}

//--Overloaded---------------------------------
std::string MSC_wdtoymd(int wd)
{
  int ww=int(wd/10);
  int dd=wd-ww*10;
  return MSC_wdtoymd(ww, dd);
}

//--Overloaded---------------------------------
std::string MSC_wdtoymd(std::string wwwwd)
{
  int wd=std::stoi(wwwwd);
  return MSC_wdtoymd(wd);
}
































































