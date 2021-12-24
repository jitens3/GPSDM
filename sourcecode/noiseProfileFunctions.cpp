#include "noiseProfileFunctions.h"
/*
170902.
Contains "noise profile" functions:

-Autocorrelation function (ACF)
  Skips frist jmin points, where jmin=diff

-Allan Variance (AVAR)
  ** In units of (ns/epoch)^2 = 1.111*10-21

-Power spectal density (PSD)
  Really, just |FT|^2 - PSD contains factors (J and tau_0) that just cancel
  when we use it anyway.
  NOTE: I only calculate the first half of the PSD (up to Niquist(?) frequncy)
  This is because it is symmetric - NOT including the
  This means the first- and second-order differenced versions have 1 less
  data point (as I skip the zero point).
  ** In units of 1.04 x 10^(-20) s^3

-Histogram
  ** Given in units of: jumps per epoch per bin / bin-width
  i.e., ~probability for given epoch to fall in given bin(?)


XXX NONE of these functions (so far) deal with
  a) outliers
  b) missing data points.
a is probably fine, since that's all part of the noise.
b, however, may be a problem. For now, I'll just ignore entire days
if they have much missing data? Not obvious how to deal with this robustly.

*/

//******************************************************************************
int allanVariance(double d[], int M, double av[], int N, bool para)
/*
Calculates Allan Variance of data, d. ADDS to array av!
Allan variance given in units [ns^2 / epochs^2].
--> to get actual AV, multiply by 1.111x10^-21
NOTE: av array MUST be initialised prior to use, since it is 'added'!
  M:  dimention (length) of array d (number of input data points)
  N:  Length of av[]; how many points in avar to calculate.
      NOTE: must have N<=M/2 !
      If too large N, or if N=0 given, uses max possible (M/2)
Only for 0-differenced data!
*/
{

  //check if N is acceptible. If not, re-size it:
  {
    int tmax=(int)floor(M/2);
    if(N>tmax||N==0)N=tmax;
  }

  av[0]=0;//nb: AVAR of 0 lag undefined (or =0)
  #pragma omp parallel for if (para)
  for(int t=1; t<N; t++){ //double check?
    double s2t=0;
    int jmax=M-2*t;
    for(int j=0; j<jmax; j++){
      s2t+=pow(d[j]-2*d[j+t]+d[j+2*t],2);
    }
    av[t]=s2t/(2*t*t*(jmax));
  }
  return 0;
}



//******************************************************************************
int autoCorrelation(double d[], int M, double acf[], int N, int jmin,
                    bool para)
/*
Calculate auto correlation.
jmin=0 for 0 diff, 1 for 1 diff etc. (skips these 'zero' points)
N<=M, can be much smaller
*/
{

  //calculated standard variance:
  double a=0;
  for(int j=jmin; j<M; j++){
    a+=d[j]*d[j];
  }
  double s2=a/(M-jmin);
  acf[0]=1;

  if(N==M||N==M-1)N=M-jmin;

  #pragma omp parallel for if (para)
  for(int t=1; t<N; t++){
    a=0;
    int jmax=M-t;
    for(int j=jmin; j<jmax; j++){
      a+=d[j]*d[j+t];
    }
    acf[t]=a/(s2*(jmax-jmin));
  }
  for(int t=M-jmin; t<M; t++) acf[t]=0;
  return 0;
}


//******************************************************************************
int powerSpectrumFFT(double d[], int M, double psd[], int N, int dif,
                  gsl_fft_real_workspace * work,
                  gsl_fft_real_wavetable * real_wt
                  )
/*
NOTE: Actually, just square of FT. Doesn't include 'factors'.
To convert to "actual" PSD [in units s^3 (or s^2/Hz), multiply by:
X (t_0/J) x (ns)^2
= (30s / 2880) * [10^-9 s]^2
= 1.04 x 10^(-20) s^3

NOTE: for 1 and 2-order diff, the PSD is actually SHORTER than 'M'!
--> jmin=0 for 0 diff, 1 for 1 diff etc.
Just the first half! NOTE: length is differenent depending on 0,1,2 order
differencing! (it's the LAST points which are excluded)

M is total length of d[]. Length of actual data is M - dif

Uses GSL FFT routines:
https://www.gnu.org/software/gsl/manual/html_node/Fast-Fourier-Transforms.html#Fast-Fourier-Transforms
https://www.gnu.org/software/gsl/manual/html_node/Mixed_002dradix-FFT-routines-for-real-data.html
work is the GSL workspace, and real_wt is the GSL trig look-up table

*/
{
  //first 'half' of PSD!
  int kmax=int((M-dif)/2)+1;

  //Note: gsl_fft_real_transform is destructive, so we need a 'temp' array
  //ALSO: d[] has leading zeros (if differenced), but we need trailing zeros
  double temp_d[2880];
  for(int j=0; j<2880-dif; j++){
    temp_d[j] = d[j+dif];
  }
  for(int j=2880-dif; j<2880; j++)
    temp_d[j] = 0; // Prob not needed..

  //Fourier transform random "phase" data:
  //Output temp_d will be in the GSL "Half Complex" form
  int stride = 1;
  gsl_fft_real_transform (temp_d, stride, M - dif, real_wt, work);

  psd[0] = pow(temp_d[0],2); //DC term
  for(int k=1; k<kmax; k++){
    double re = temp_d[2*k-1];
    double im = temp_d[2*k];
    psd[k] = re*re + im*im; //note: doesn't include factors
  }
  for(int k=kmax; k<N; k++) psd[k]=0; //shouldn't be included in calcs!

  return 0;
}





//******************************************************************************
int histogram(double d[], int M, double hist[], int N, double w, int jmin)
/*
Histogram:
  M=number data points
  N=number of bins - MUST be an odd number!! add a check?
  w=bin width
  Given in units: jumps per epoch per bin / bin-width
*/
{

  int B=1+(N-1)/2; //number of "positive" bins.

  double x=1./(w*(M-jmin)); // 'per epoch'

  for(int j=jmin; j<M; j++){
    double s=d[j];
    int n=0;
    if(s>0)n=N-1; //if not in one of the bins, add to first/last bin
    for(int b=0; b<B; b++){//b denote 'positive' bins. there are N=2B+1 bins.
      if(fabs(s)<(b+0.5)*w){
        if(s>0) n=B-1+b;
        else n=B-1-b;
        break;
      }
    }
    hist[n]+=x;
  }
  return 0;
}
