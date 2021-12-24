<<<<<<< HEAD
/*
This function is based of of readJplBiasData with changes that reflect the differences in the 1s data files compared to the 30s data files. Reads in the clk_1s data files, sorts the data and fills the appropriate arrays with data.

It searches the files and tally's the number of newline '\n' characters. After 5, the header is over and data fills the rest of the file. The function then reads the clock data, sorts it and puts it in the relevant arrays. Since clock positions are not necessary for the GW events, the data is sorted by epoch.

int readJplBiasData_1s(std::string path)

INPUT:
	path		:: directory that hold JPL 1s data files
	event		:: name of GW event (ex: GW180817)
	day		:: day related to event

OUTPUT:
	bias, ferr	:: array for clock bias and formal error data
	num_clocks	:: number of total clocks
	clockNames	:: array of clock names in alphabetical order
	clockType	:: array of clock type matching order in clockNames (R-receiver, S-satellite)	




===== CHANGE LOG =====



===== TO DO =====




*/
#include "likelihoods.h"
#include "JplGpsDataClass.h"
#include "miscFunctions.h"
#include "NumericCdfInverseClass.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

int readJplBiasData_1s(std::string path, std::string in_event, int in_day);

int main()
{

	JplGpsData data;

	readJplBiasData_1s("../jpl1s/unr1s/GW170817/unr19625_completeDataSet.txt","",0);
	//readJplBiasData_1s("../jpl1s/unr1s/","GW180817",19621);

	//cout<< "It Worked." <<endl;
	return 0;
}

int readJplBiasData_1s(std::string path, std::string in_event, int in_day)
=======
//class JplGpsData::
#include "JplGpsDataClass.h"

int JplGpsData::readJplBiasData_1s(std::string path, std::string in_event, int in_day)
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
{
	//Sets GW event name and day
	std::string event=in_event;
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
<<<<<<< HEAD
	std::string buffLine;	
=======
	std::string buffLine;
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	while(std::getline(inFile,buffLine)){
		iLines++;
	}
	inFile.close(); //Must re-open
<<<<<<< HEAD
	
=======

>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	std::ifstream inFile1(jpl_file_dir);
	//Loads each line of the file into a string array
	std::vector<std::string> fileLine(iLines);
	std::string sLine;
	int iLine=0;
<<<<<<< HEAD
	
=======

>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	while(std::getline(inFile1,sLine)){
		if(sLine.size()>0)
		fileLine[iLine]=sLine;
		//std::cout<<fileLine[iLine]<<std::endl;
		iLine++;
	}
	inFile1.close(); //Don't need file anymore

	//Variables that will store file information
<<<<<<< HEAD
	char rs; //Receiver or Sat clock
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

	int num_epochs=end_epoch-begin_epoch+1;	
=======
	// char rs; //Receiver or Sat clock
  std::string rs;
	std::string clockName;
	int wk, dy, gpsec;
	long double dbias, sigma;
	//double nano=pow(10.,9.);

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
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6

	std::cout<<"The first epoch is: "<<begin_epoch<< std::endl;
	std::cout<<"The last epoch is: "<<end_epoch<< std::endl;
	std::cout<<"The total number of epochs is: "<<num_epochs<< std::endl;

	//Count number of clocks
	std::string tempClock1=" ";
<<<<<<< HEAD
	int num_clocks=0;
=======
	int temp_num_clocks=0;
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	for(int i=0; i<iLines; i++){
		std::stringstream tempin(fileLine[i]);
		tempin>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
		if(tempClock1 != clockName){
			tempClock1 = clockName;
<<<<<<< HEAD
			num_clocks++;
		}
	}
	std::cout<<"The number of clocks is: "<<num_clocks<<std::endl;

	//Store clock names and clock types
	std::string clockNames[num_clocks];
	std::string clockType[num_clocks];
=======
			temp_num_clocks++;
		}
	}
  num_clocks = temp_num_clocks; // num_clocks is a class variable
	std::cout<<"The number of clocks is: "<<num_clocks<<std::endl;

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
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	std::string tempClock2=" ";
	int count=0;
	for(int i=0; i<iLines; i++){
		std::stringstream tempin1(fileLine[i]);
		tempin1>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;
		if(tempClock2 != clockName){
			tempClock2 = clockName;
<<<<<<< HEAD
			clockNames[count]=tempClock2;
			clockType[count]=rs;
=======
			// clockNames[count]=tempClock2; //
			// clockType[count]=rs;
      if(rs=="R"){ //XXX I am trying to check if it's a reciever. Right?
        prn[count]=tempClock2;
        svn[count]="00"; //SVN for stations is 0
        blk[count]="AR"; //'AR' is what I call 'block' for stations
        clk[count]="**"; //we don't know yet which is which clock
      }else{
        prn[count]="**";
        svn[count]=tempClock2; //This is the SVN right??
        blk[count]="**";
        clk[count]="**"; //we don't know yet which is which clock
      }
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
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
<<<<<<< HEAD
	double bias[num_clocks][num_epochs];
	float ferr[num_clocks][num_epochs];

=======
	// double bias[num_clocks][num_epochs];
	// float ferr[num_clocks][num_epochs];
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
	//Set values of bias ferr arrays
	for(int i = 0; i<num_clocks; i++){
		for(int j=0; j<num_epochs; j++){
			std::stringstream tempin2(fileLine[i]);
<<<<<<< HEAD
			tempin2>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;	
		
			bias[i][j]=dbias;
			ferr[i][j]=sigma;		
		}
	}
	std::cout<<bias[1][1]<<"   "<<ferr[1][1]<<std::endl;
=======
			tempin2>>rs>>clockName>>wk>>dy>>gpsec>>dbias>>sigma;

			bias[i][j]=dbias;
			ferr[i][j]=sigma;
		}
	}
	//std::cout<<bias[1][1]<<"   "<<ferr[1][1]<<std::endl;
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6

/*	//print out bias and ferr for each clock
	std::cout<<"\nClock:         BIAS           FERR   "<<std::endl;
	for(int i = 0; i<num_clocks; i++){
		for(int j=0; j<num_epochs; j++){
<<<<<<< HEAD
			std::cout<<"Clock["<<i+1<<"]:  "<<bias[i][j]<<"    "<<ferr[i][j]<<std::endl;			
=======
			std::cout<<"Clock["<<i+1<<"]:  "<<bias[i][j]<<"    "<<ferr[i][j]<<std::endl;
>>>>>>> 0a172f5bf6cb437420fc25458516b168851321d6
		}
	}
*/

	//reSizeVec(bias,num_clocks,num_epochs);
	//reSizeVec(ferr,num_clocks,num_epochs);

	return 0;
}
