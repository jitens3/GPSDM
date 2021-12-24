# GPS.DM Bayesian dark matter search code

Short descriptions for what each program does.
More detail is given inside the source code, located in /sourcecode directory.

  * For a detailed description of the physics/statistics behind how the code
    and method works, see the "BayesianMethod" pdf in /documentation
    directory (see also the "VectorVelocityDistribution" memo there).

  * For a detailed overview of the noise properties of the GPS clocks, see the
    "NoiseProfilePlots" memo.

  * Full documentation of the code is given inside the source code, at the start
    of each file.

  * Brief documentation for each main program/function is given in the
    "documentation" memo. This document also maps the variable names to the
    notation used in the papers/memos.

********************************************************************************
## Getting the code from github:

All the code is on the GitHub repository. For now, this is a private repository,
which means you need to have a GitHub account to access is (which is free), and
you need to be a "collaborator". Only I can make you a collaborator, so email me
if you want access.

The best way (especially when using a remote server) is to use git.  

**See the full "Using GitHub" instructions for more info.**

If you use the "Desktop" (GUI) version, it's all pretty straight forward.
Click "Clone repository", the name is "benroberts999/GPSDM"

From the command-line:
Initially, you have to 'clone' the repository
  * $ git clone https://github.com/benroberts999/GPSDM.git LOCALDIR

where you should change 'LOCALDIR' to whichever local directory you want the
code to go. This directory shouldn't already exist. You can leave this option
blank, and git will call the directory 'GPSDM'.
Once you have already done this, if you wish to download new/updated version of
the code, you can "pull" any updates down from the remote repository
  * $ git pull

This will only update files that have changed since you last pulled them.
It won't over-write any changes that you've made.
If you do want to overwrite (i.e. discard) all of your changes, use:
  * $ git reset --hard

be careful with this, it will overide all local changes.

This requires the program 'git' be installed.
To install git (on ubuntu), type
  * $ sudo apt-get install git

The git software is installed both on the Oasis and QuickHydrogen servers.

For more info on using git, see the "UsingGit" document I sent around
(also in memos folder on dropbox), or just google around, there are lots of
great resources.


********************************************************************************
## Compilation + Running

Either, run
  * To compile all: _$ make_
  * To compile signle program (program_name): _$ make program_name.x_

All programs are compiled with g++, but you can change this in the make files
(other complilers not tested).
To install g++ (on ubuntu):
  * $ sudo apt-get install g++

The codes make use of a few c++11 features. This means very old compilers won't
be able to compile the code. This is unlikely to be an issue, c++11 has been
standard for many years now. If you have some troubles, just make sure the
compiler is up to date.

All programs make use of the open-source GSL libraries (GNU Scientific
Libraries).
On ubuntu these can be installed with:
  * $ sudo apt-get install libgsl0-dev

You (may) need to also download the lapack and blas libraries.
These can be installed (again, ubuntu) with following command:
  * $ sudo apt-get install libblas-dev libatlas-dev liblapack-dev

Also (optionally) uses gnuplot for some things:
install gnuplot: (on ubuntu)
  * $ sudo apt-get install gnuplot-x11

[if that doesn't work, try just "install gnuplot" (without -x11), or both]


Some programs make use of openMP to parallelise the computations.
This is included in g++, so you don't need to install anything.
However, this can use up all the computers resources, meaning the computer may
run very slowly for other tasks.
If you want to actually use your computer while these codes
are running, it's probably best to limit the number of threads openMP can use.
For some of the resource-heavy programs (in particular, testLikelihoods), you
can do this with an option in the input .dat file.
Or, you can switch of openMP - this can be done in the
"make" files, change "yes" to "no" for the openMP option - should be obvious
once you open the make file.


### Using the code

  * **Run** the program by typing: _$ ./PROGRAM.x_
  * **Input parameters:** All the input parameters are in the file
      "PROGRAM.dat", which is read upon running.

For all programs, detailed descriptions of the input parameters are given
inside the PROGRAM.dat file.

There are many functions located in various files.
Thorough descriptions are available inside the sourcecode for the each of the
given functions.


**NOTE**: the file: "PRN_GPS_GPSDM.txt"
needs to be placed in the same directory as the input jpl 30s files!
See below for instructions.

********************************************************************************
## Some important notes

### PRN to SVN mappings

**NOTE:**
Make sure the updated version of the PRN_GPS file, _PRN_GPS_GPSDM.txt_ is
located in the same directory as the JPL data files.

JPL has a file, PRN_GPS,
(ftp://sideshow.jpl.nasa.gov/pub/gipsy_products/gipsy_params/PRN_GPS.gz)
which has the mapping of each PRN to the specific SVNs.
However, this file does not accurately reflect the true clock assignemtns (i.e.,
Rb or Cs).
So, we have updated the original PRN_GPS file with the accurate clock
assignments from the "operational Advisories" (OAs):
https://www.navcen.uscg.gov/?Do=gpsArchives
The new file is called "PRN_GPS_GPSDM.txt".

This file is generated automatically by a python script, updateSVNPRNmap.py.
Note, this script is located (for now) in a sepperate (public) repository:
https://github.com/benroberts999/updatePRN_SVN

  * See the benroberts999/updatePRN_SVN read-me for more info.

(PRN, Psuedo-Random Noise code, is a number 1-32. Only one satellite at a time
has a given PRN, but the PRN ascociated with any particular satellite may
change. The SVN, Space Vehicle Number, is a unique identifier for each
sattelite. Each SVN may have either a Rb or Cs clock, only one of which is in
operation at any given time.)


### ECI position files

JPL has position files available on their website, however, they are in the
'fixed earth surphace' non-inertial frame.
We need the positions in the Earth-Centred Inertial (ECI) frame.
Geoff has generated these using the GYPSY software for weeks 1280-1881.
Ideally, we would want to automatically generate these files whenever new
data from JPL comes available. It might depend on the GYPSE licence??

-Also: the way the code works is a little dependent on the format of these
files. If it changes, we'd need to update the code.
Also, some of the ECI files seem to be missing some data??


### JPL clock data files, possible future changes

Written for RINEX version 3.00.

We need to be a little wary that if JPL updates the file format, we'll also have
to update that part of the program (only 1 function, so shouldn't matter too much).

Before week 1934, the formal error every 10th epoch was incorrect.
The code dealt with that. After 1934, that error is fixed. The code
knows this too.

Also after week 1934, JPL's naming convention changed [from clk_30s to clk].
Note: New ".clk" files do not include the 'second' midnight; the earlier files
still do. This shouldn't matter.

********************************************************************************
## Issues + still To-Do

  * Automate ECI position file
  * sHs=0 issue
  * Calculate L(x) for each x (?)
  * Large cross-clock correlations properly
  * Include base stations into GPS simulator
  * Read in +/- 1 day of data!? [power of two? nah.] ~100 epochs?
  * More on github.com

********************************************************************************
# Each Program:

********************************************************************************
## plotClockData

Plots the clock solutions for a given day.
Can plot all the clocks, or only a particular set, or a particular PRN/station.

Creates gnuplot files, and places them in a new directory called:
output-plotClockData/
[unless that directory doesn't exist and couldn't be created, in which case
it puts them in the current directory!]
-Can be deleted immediately if you no longer need them.

Either plots "to screen" (using gnuplot),
or to a pdf file (using gnuplot and pdflatex)
If you don't have gnuplot/latex installed, it will still work, it will just
create the .gnu/.tex files without actually plotting them.

Uses 'system' to call gnuplot/latex, so will only fully work on linux machines.
However, you can just call gnuplot/latex sepperately to plot to produced files.

Can plot raw data, polynomial-reduced data, or differenced data.
Plots either all clocks, or a specific clock (by PRN/CLK/Block)



********************************************************************************
## testLiklihoods

Still under development. It's a little messy, but it works.

This program is designed to test the "likelihoods" odds ratio program.
   * It generates random time-series for a user-specified number of clocks.
The noise can be either white (with given standard-deviation), or it can use the
known PSDs from GPS clocks to simulate realistic time series.
It can do this using either the overall averaged PSD, or by randomly assigning
each satellite an SVN, and using the correct PSD for that SVN.
   * It places each satellite in a realistic (though static) position in space
   * Either inject ideal signals into the data (based on user input) to test
     how well program can pick up events, or
   * not inject events, to test for false-positives
   * Then calls 'MClikelihoods' program, which calculates the odds ratio.

Note: you may limit the number of cores you use for parallelisation.
If you are running tests on one of the servers this is a good idea, so you
don't just use up all the servers resources. It is also a good idea if you are
running the code on your own PC to keep at least 1 core free, so you can still
use the PC for other tasks without it being too slow.

The program optionally makes use of power spectrums (for simulating clocks) and
the inverse covariance functions (for covariance matrix). These have been
calculated already, and can be downloaded from github (see main README for
more detailed git instructions)
   * https://github.com/benroberts999/PowerSpectrums

e.g.:
   * $ wget http://github.com/benroberts999/PowerSpectrums/archive/master.tar.gz
   * $ tar xf master.tar.gz -C ./ --strip-components=1

NOTE: At the moment, vmin and vmax are hard-coded in.
J_W (the data "window") is worked out using vmin.
Jw (number of points in Chi^2 sum) is important.
Open question...

********************************************************************************
## skewness

Calculates the skewness (k3) and kurtosis (k4) for GPS by SVN and block.
Also, outputs on daily basis, so look for annual modulation.

********************************************************************************
## impTwWindow

Program searches for "obvious" thin-wall events, using single-differenced data.
Looks for "windows" (in time) in which a certain proprtion of the clocks
experience jumps in both directions of a certain magnitude.
Loops over different time window lengths (J_W), and 'jump' height ranges.
Method described in: http://arxiv.org/abs/1704.06844

Outputs the total number of "events" it found for each window, and jump height.

It requires the 'ACF' files to exist, in order to know the standard deviation
for each clock. If they don't exist, will assume the worst s.d. for each clock.
They can be downloaded from github (benroberts999/PowerSpectrums)

**NOTE** At the moment, uses the averaged standard deviations that are averaged
over SVNs (from the ACF file).
This could easily be changes to either:
 a. use the actual s.d. from this day, or
 b. use the averaged s.d. (from ACF file) for the specific SVN...

This program was used for the results in:
http://arxiv.org/abs/1704.06844

Actually, this is a slightly updated version. However, the results didn't
change.

********************************************************************************
## impTwWindowPattern

Similar to impTwWindow, but this program has only a single window size
(50 epochs), and is designed to output 'pattern' files, for Conner/Wyatt
mathematica "pattern matching" program.
The actual window used by the program is 101 epochs long. But, since we always
give the output with j_r (when the reference clock jumped) in the middle, the
actual maximum window is ~50 epochs, which is still large.

It requires the 'ACF' files to exist, in order to know the standard deviation
for each clock. If they don't exist, will assume the worst s.d. for each clock.
As in "impTwWindow".
They can be downloaded from github (benroberts999/PowerSpectrums)

Format of output file is:
  * date jplFile Clock REF jmin jmax jR
  * S_cut
  * numClocks numEpochs
  * r.n(rproj) rx ry rz sig svn-prn-clk-blk ...data...
  * [one line for each clock; the ref clock is always the last entry]

At the moment, it doesn't seem to work with openMP.
This probably can be fixed easily, however, it doesn't matter much, since
the bottle-neck is reading in the files, and it's fairly quick.


********************************************************************************
## processNoise

Calculates power spectrums, autocorrelation functions, histograms,
Allan variance, and standard deviation (with uncertainty), for each SVN, clock,
and reference clock combination.
Then, averages those functions, to find the above for each block (averaged over)
SVNs, and averaged over all reference clocks.
Also, it can "swap" the reference clock to one of the other satellite clocks.

Note: The list of possible reference clocks is hard-coded in.
It is easy to update.
Probably, this could be automated.
Also hard-coded in, is which SVNs belong to which satellite block.
When block III satellites are launched, this needs to be updated [svnToBlock()].
If not, the newer satellites will be interpreted as block IIF.

All output files are placed into a directory specified by the user.
The output files are names like:
**ClkBlk-svn-REF-outlabel.ext**.
e.g.,
  * RbIIR-61-USN3-test1.acf

holds the autocorrelation function for the Rb IIR clock with svn=61, for days when USN3 reference was used.
Note: usually 'outlabel' will be the weeks used to generate the functions, e.g. '1280-1959'.

Other information about the run is written to the header lines of each file.
The header lines are marked with an '#'.
Generally, the data is given as a table, each row corresponding to an epoch/lag/
frequency component, and each column for 0th, 1st or 2nd order differencing.

The standard deviations are written to the ACF files.
The summary of all standard deviations (including the uncertainty in those
standard deviations) are written to another file:
**sdWunc-outlabel.txt**


********************************************************************************
## stationNoise

Calculates the standard deviations for the station clocks.
(only standard deviations for now, can update for others, ACF/PSD etc later)

 * For now, only does single reference clock as a time (will skip days that
   don't use that reference clock)
 * For each stations clock, calculates s.d. as per usual [will include
 contribution from the station clock and the reference clock]
 * Also calculates std dev. for the reference clock by calculating
 cross-correlations. This should include only a contribution from the reference
 clock.

********************************************************************************
## dateConvert

Just a short front-end program that converts between JPL week-day date
to yyyy-mm-dd dates.

Input for the program is given upon running the program.

E.g.,
  * $./dateConvert 19005

will convert the JPL date format week=1900, day=5 and output "2016-06-10"

  * $./dateConvert 2016 6 10

will output "19005"

**Note:** program is not that smart, and may give rubbish results if input given in
incorrect format. e.g.
  * $ ./dateConvert 20160610

will interpret this as a wwwwd format, and output "1965-07-13".
This would be simple to fix...but I haven't yet.



********************************************************************************
## fetchAndCheckJPL30s

Simple program that downloads all the JPL 30s files from
ftp://sideshow.jpl.nasa.gov/pub/jpligsac/
Uses wget, so that part will only work on linux machines.
Also, unzips the files (again, linux only).
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
    the ECI files (this is a problem) etc. Note: we currently only have ECI
    files for weeks 1280-1881; this range is hard-coded in, and should be updated!
  * **PRN_GPS-faillog.txt**
    (Actually produced inside JplGpsDataClass::mapSVNPRN):
    Any PRN for which the clock/svn wasn't identified (meaning it is missing
    from the PRN_GPS_GPSDM file).

**Note**
Using wget in this way is dangerous!
If the program is quit half-way through downloading a file, everything might
look OK, but actually only half the file downloaded.
I consider this a bug, which I will try to work out how to fix.


********************************************************************************
## convertmout

Converts JPL 30s clock data into a form that can be easily read by mathematica.
Two output formats:
  * a) ".mout" format, suitable for input to Andrei's MMA notebooks
  * b) A simpler ".mout2" format, that is a simple table

Format for ".mout2"/simple version:
1. date (human readable)
2. PRN
3. SVN
4. Block
5. Clock type
6. (6+): Bias data
  * Each column is a different clock, each row an epoch.
  * Does not output the formal error, but can be 'normalised' by it.
  * Good for just checking one clock. The students prefer it this way.

Uses openMP for parallelisation. Parallel by week.
The location of the input JPL files, and all other options are given inside
the input file: "convertmout.dat".
Detailed info on each of the options is also given inside the .dat file.
