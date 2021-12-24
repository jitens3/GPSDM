ifndef _SEARCHDATA_H
#define _SEARCHDATA_H
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm> //used for "sort"
#include <omp.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "mathematicsFunctions.h"
#include "miscFunctions.h"
#include "JplGpsDataClass.h"


class SearchData{

	friend std::ofstream & operator<<(std::ofstream & of, SearchData & sd);

public:

	SearchData(); //default ctor
	SearchData(std::string path="", int in_jplweek=0, int in_jplday=0,
                       bool detrend=false, bool difference=false); //param ctor
	~SearchData(); //dtor



	void serialize(std::ofstream & of);


	JplGpsData data;
	double snr[2880];
	double odds[2][2880];
	std::string refclk;



private:


};




#endif