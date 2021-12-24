/*
2017-03-20

Just a short front-end program that converts between JPL week-day date 
to yyyy-mm-dd dates.

Input for the program is given upon running the program.

E.g.,
  * $./dateConvert 19005
will convert the JPL date format week=1900, day=5 and output "2016-06-10"

  * $./dateConvert 2016 6 10
will output "19005"

Note: program is not that smart, and may give rubbish results if input given in
incorrect format. e.g.
  * $./dateConvert 20160610
will interperet this as a wwwwd format, and output "1965-07-13".
This would be simple to fix...but I haven't yet.
  
  
*/
#include "miscFunctions.h"


int main ( int numIn, char* argv[] )
{

  std::string input;
  if(numIn==2){
    std::cout<<"JPL wwwwd to yyyy-mm-dd\n";
    int wd=std::stoi(argv[1]);
    std::cout<<MSC_wdtoymd(wd)<<std::endl;
  }else if(numIn==4){
    std::cout<<"yyyy-mm-dd to JPL wwwwd\n";
    int ya=std::stoi(argv[1]);
    int ma=std::stoi(argv[2]);
    int da=std::stoi(argv[3]);
    std::cout<<MSC_ymdtowd(ya,ma,da)<<std::endl;
  }else{
    std::cout<<"I don't understand the input: ";
    for(int i=1;i<numIn;i++){
      std::cout<<argv[i]<<" ";
    }
    std::cout<<".\n";
    std::cout<<"Either enter a JPL week in the format 'wwwwd',\n"
        <<"or a date in the format 'yyyy mm dd'.\n";
  }
  
  return 0;
}






































