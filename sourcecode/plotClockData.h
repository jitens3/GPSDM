#ifndef _PLOTCLK_H
#define _PLOTCLK_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "JplGpsDataClass.h"
#include "miscFunctions.h"
#include "mathematicsFunctions.h"

int plotClocks(JplGpsData data, std::string oneclock, std::string base, 
               std::string filename, double R2min, int jmin, int jmax, 
               std::string title, 
               int iCi, int iCf, int iPlot, int iStyle, double dOfs);

#endif
