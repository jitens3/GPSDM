#ifndef _CONVMOUT_H
#define _CONVMOUT_H
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm> //used for "sort"
#include <omp.h> 
#include "JplGpsDataClass.h"
#include "mathematicsFunctions.h"
#include "miscFunctions.h"

int writeMout(JplGpsData data, int prec, int iFormat, int iNorm);

#endif //here? or at end?
