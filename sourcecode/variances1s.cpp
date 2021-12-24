#include "likelihoods.h"
#include "JplGpsDataClass.h"
#include "miscFunctions.h"
#include "NumericCdfInverseClass.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

int main(void)
{
	//create blank GPS data object
	JplGpsData data;

	data.differenceData(1);
	clkstd = data.calculateStdDev();
	return 0;
}
