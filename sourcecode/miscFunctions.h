#ifndef _MISCFUNS_H
#define _MISCFUNS_H
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>  //?? can I avoid this?
#include <omp.h>

bool MSC_fexists(const std::string& filename);

bool MSC_direxists(const std::string& DIR);

int MSC_execute(std::string exe);

void MSC_flag(int s1, bool print=true);

void MSC_progressBar(std::string text, int barWidth, double current=0,
                    double total=0, bool OMP=false);

std::string MSC_padIntString(int in_int, int digits=2);

//convert between JPL week+day to yyyy-mm-dd format:
int MSC_ymdtowd(int yy, int mm, int dd);
int MSC_ymdtowd(std::string yyyymmdd);
std::string MSC_wdtoymd(int w, int d);
std::string MSC_wdtoymd(int wd);
std::string MSC_wdtoymd(std::string wwwwd);

#endif
