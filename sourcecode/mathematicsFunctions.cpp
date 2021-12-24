#include "mathematicsFunctions.h"
/*

170829.


Note: abs() only for ints! fabs() for doubles!!! (std!!!!!)

All functions begin with MFS_ designator

Some general "Mathematics functions", including
-wrapper for GSL wieghted least squares
-Random number generators
-Matrix equations
-"Fast" versions of exp, erf, inverse erf

XXX the random number generators need further testing, particularly with
other compilers.

*/


//******************************************************************************
bool MFS_polywfit(int obs,int degree, std::vector<double> dx,
              std::vector<double> dy, std::vector<double> dw,
              std::vector<double> &coef, double &chisq, double &R2)
/*
Polynomal Fit, weighted least squares.
Wrapper function that calls GSL libraries to perform an any-order weighted
least squares polynomial fit.
Note: this function should only really be called by the polynomialDetrend()
function (inside JplGpsData class).

INPUT:
 obs: number of data points (observations)
 degree: degree of polynomial to be fitted
 dx: array of 'x' data points
 dy: array for the 'y' observations (data points)
 dw: array for the weights (NOT sigma!) for weighted fit
OUTPUT:
 coef: array (order degree+1) of the polynomial coef
 chisq: Chi-squared value (output from GSL)
 R2:  R^2 value: calculated manually
.
Chi squared is a little unreliable when the weightings are important??
 Particularly, if some of them are zero!
------
===== Change Log =====
160310- Fixed memory issue: Freed memory of 'w'!!
170702- Changed coef to std::vector
*/
{
  degree=degree+1; //"normal language -> array language"
  gsl_multifit_linear_workspace *work;
  gsl_matrix *cov, *x;
  gsl_vector *y, *w, *c;
  x = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  w = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
  for(int i=0; i < obs; i++) {
    for(int j=0; j < degree; j++) {
      gsl_matrix_set(x, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
    gsl_vector_set(w, i, dw[i]); //Weight vector. (Weight already! not sigma)
  }
  work = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_wlinear (x, w, y, c, cov,&chisq, work);

  /* store result ... */
  for(int i=0; i < degree; i++)
  {
    coef[i] = gsl_vector_get(c, i);
  }

  double ybar=0,Wtot=0;
  for(int i=0;i<obs;i++){
    ybar+=dy[i]*dw[i];
    Wtot+=dw[i];
  }
  ybar=ybar/Wtot; //weighted average value of y (observed, not fit-function)
  double SStot=0,SSres=0;
  for(int i=0;i<obs;i++){
    double funi=0;
    for(int j=0;j<degree;j++){
      funi+=pow(dx[i],j)*coef[j]; //value of the best-fit function
    }
    SStot+=pow((dy[i]-ybar),2)*dw[i];
    SSres+=pow((dy[i]-funi),2)*dw[i];
  }
  R2=1-SSres/SStot; //XXX NOT sure if this is OK (the weighted least squares)

  gsl_multifit_linear_free(work);
  gsl_matrix_free(x);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  return true;
}



//******************************************************************************
double MFS_randGausVal(double x0, double sig)
/*
170713.
Uses inverse transform sampling to create a Gaussian random number.
Returns a random double that is drawn from a Gaussian PDF, with mean x0,
and standard deviation 'sig'.
Makes use of my 'randDouble' function.
It also makes use of my "fastInvErf" function.
MUCH faster than previous GSL version (many orders of magnitude)!!
Also, using 'fastInvErf' is ~3x faster than inverseErf.
Output was tested using 10,000 points. Histogram matched exactly! :)
Works with sig=0 (just returns x0)!
=== Change Log ===
*/
{
 double u=MFS_randDouble(0,1);  //uniform u
 return x0 + 1.41421*sig*MFS_fastInvErf(2*u-1);
}



//thread_local gsl_rng * GSL_RAND_R = gsl_rng_alloc (gsl_rng_mt19937);
//gsl_rng_mt19937
//gsl_rng * GSL_RAND_R = gsl_rng_alloc (gsl_rng_taus);
//unsigned long int GSL_RAND_MAX = gsl_rng_max(GSL_RAND_R);
//******************************************************************************
double MFS_randDouble(double a, double b)
/*
170714. Updated to use c++11
Returns a uniformly distributed real number between a and b.
Uses 'std::mt19937', a c++11 Mersenne Twister pseudo-random generator.
Also uses <thread>, and is safe to use with OpenMP parallelisation
(nb: normal rand() function is not!).
Needs checking with other compliers, works with GCC.
It does not need to be seeded, that is taken care of automatically.

The GSL version is extremely slow on multi-thread machines...???

http://www.cplusplus.com/reference/random/mt19937/

=== Change Log ===
*/
{
  thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> distribution(a, b);
  return distribution(generator);

  //double f = (double)rand() / RAND_MAX;
  //return a + f*(b - a);
  //return a + ((b-a)/GSL_RAND_MAX)*gsl_rng_get(GSL_RAND_R);
  //return a + (b-a)*gsl_rng_uniform(GSL_RAND_R);
}



//******************************************************************************
int MFS_invertMatrix(std::vector< std::vector<double> > inmat,
                        std::vector< std::vector<double> > &outmat)
/*
170316.
Will invert a square matrix of /any/ dimension, n.

https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
https://www.gnu.org/software/gsl/manual/html_node/LU-Decomposition.html
http://www.macapp.net/MyWikiThings/invertmatrix.c
Requires #include <gsl/gsl_linalg.h>

INPUT:
  inmat   :: double 'flat' matix of dimension n*n [call as "(double *)inmat"]
OUTPUT:
  outmat  :: inverted 'flat' matrix [call as "(double *)outmat"]

See also overloaded version, that just takes in 1 matrix, and inverts it!

=== Change Log ===

*/
{
 int iRet=0;

 //size of matrix:
 int n=(int)inmat.size();
 if(inmat[0].size()!=(size_t)n)return 0; //matrix not square!

 //ensure "out" matrix (inverted matrix) is of correct dimension
 outmat.resize(n, std::vector<double>(n));

 // Define all the used matrices (for GSL)
 gsl_matrix * m  = gsl_matrix_alloc (n, n);
 gsl_matrix * inverse = gsl_matrix_alloc (n, n);
 gsl_permutation * perm = gsl_permutation_alloc (n);
 //fill matrix:
 for(int i=0;i<n;i++){
   for(int j=0;j<n;j++){
     //gsl_matrix_set(m,i,j,inmat[i*n+j]);
     gsl_matrix_set(m,i,j,inmat[i][j]);
   }
 }
 //peform LU decomposition (using GSL)
 //and inversion (if non-singular)
 int s;
 gsl_linalg_LU_decomp (m, perm, &s);
 double det=gsl_linalg_LU_det (m, s);
 //if(det!=0)
 gsl_linalg_LU_invert (m, perm, inverse);
 if(det==0)iRet=1;
 //Fill the output matrix:
 for(int i=0;i<n;i++){
   for(int j=0;j<n;j++){
     //outmat[i*n+j]=gsl_matrix_get(inverse,i,j);
     outmat[i][j]=gsl_matrix_get(inverse,i,j);
   }
 }
 //clear memory
 gsl_permutation_free (perm);
 gsl_matrix_free (m);
 gsl_matrix_free (inverse);
 return iRet;
}


//---- Overloaded: -------------------------------------------------------------
int MFS_invertMatrix(std::vector< std::vector<double> > &inmat)
/*
170827.
Overloaded version, that just over-writes the input matrix with the output!
*/
{

  //create "temporary" output (inverted) matrix:
  std::vector< std::vector<double> > outmat; //will be re-sized inside function

  //call function to do the matrix inversion.
  int iRet=MFS_invertMatrix(inmat,outmat);

  int n=(int)outmat.size(); //use size of outmat in case something went wrong

  //over-ride the input matrix with values from output:
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      inmat[i][j]=outmat[i][j];
    }
  }

  return iRet;
}

//---- Overloaded: -------------------------------------------------------------
int MFS_invertMatrix(std::vector< std::vector<float> > inmat,
                        std::vector< std::vector<float> > &outmat)
{

  int n=(int)inmat.size();
  int m=(int)inmat[0].size();
  std::vector< std::vector<double> > dbl_inmat(n,std::vector<double>(m));

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      dbl_inmat[i][j]=inmat[i][j];
    }
  }

  int iret = MFS_invertMatrix(dbl_inmat);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      //dbl_inmat is now the inverted matrix; so make it the output:
      outmat[i][j]=dbl_inmat[i][j];
    }
  }

  return iret;

}

//---- Overloaded: -------------------------------------------------------------
int MFS_invertMatrix(std::vector< std::vector<float> > &inmat)
{

  int n=(int)inmat.size();
  int m=(int)inmat[0].size();
  std::vector< std::vector<double> > dbl_inmat(n,std::vector<double>(m));

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      dbl_inmat[i][j]=inmat[i][j];
    }
  }

  int iret = MFS_invertMatrix(dbl_inmat);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      inmat[i][j]=dbl_inmat[i][j];
    }
  }

  return iret;

}



//******************************************************************************
double MFS_calcDeterminant(std::vector< std::vector<double> > inmat)
/*
170622.
Calculates the determinant of any real square matrix of dimension n.
Uses the GNU 'GSL' libraries:
https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
https://www.gnu.org/software/gsl/manual/html_node/LU-Decomposition.html
Requires #include <gsl/gsl_linalg.h>

INPUT:
  inmat     :: double matix of dimension n*n [from std::vector]
  n         :: integer, dimension of matrices

=== Change Log ===
170702- Uses std::vector input, avoid variable arrays!
*/
{

 //size of array:
 int n=(int)inmat.size();
 if(inmat[0].size()!=(size_t)n)return 0; //matrix not square!

 // Define all the used matrices (for GSL)
 gsl_matrix * m  = gsl_matrix_alloc (n, n);
 gsl_permutation * perm = gsl_permutation_alloc (n);
 //fill matrix:
 for(int i=0;i<n;i++){
   for(int j=0;j<n;j++){
     //gsl_matrix_set(m,i,j,inmat[i*n+j]);
     gsl_matrix_set(m,i,j,inmat[i][j]);
   }
 }
 //peform LU decomposition (using GSL)
 int s;
 gsl_linalg_LU_decomp (m, perm, &s);
 double det=gsl_linalg_LU_det (m, s); //XXX ok as double?
 //clear memory
 gsl_permutation_free (perm);
 gsl_matrix_free (m);
 return det;
}

//---- Overloaded: -------------------------------------------------------------
double MFS_calcDeterminant(std::vector< std::vector<float> > inmat)
{

  int n=(int)inmat.size();
  int m=(int)inmat[0].size();
  std::vector< std::vector<double> > dbl_inmat(n,std::vector<double>(m));

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      dbl_inmat[i][j]=inmat[i][j];
    }
  }

  return MFS_calcDeterminant(dbl_inmat);

}



//******************************************************************************
int MFS_linsolve(std::vector< std::vector<double> > inmat,
                 std::vector<double> invec, std::vector<double> &outvec)
/*
170321.
Solves the linear matrix equation A.x=b for x
Where:
  A=inmat    is an n*n square matrix
  x=outvec   is an n dim vector (the answer/output!)
  b=invec    is an n dim vector (the input)

Uses the GNU 'GSL' libraries:
https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
https://www.gnu.org/software/gsl/manual/html_node/LU-Decomposition.html
Requires #include <gsl/gsl_linalg.h>

INPUT:
  inmat  :: double 'flat' matix of dimension n*n [call as "(double *)inmat"]
  invec  :: double vector of dimension n
  n      :: integer, dimension of matrices
OUTPUT:
  outvec :: solution. output n dimensional vector


*/
{

  int n=(int)inmat.size();

  if(inmat[0].size()!=(size_t)n)return 1; //matrix not square!
  if(invec.size()!=(size_t)n)return 1; // invec incorrect dimension!

  //ensure outvec has correct dimension:
  outvec.resize(n);

  int iRet=0;
  // Define all the used matrices/vectors:
  gsl_matrix *A = gsl_matrix_alloc (n,n);
  gsl_vector *b = gsl_vector_alloc (n);
  gsl_vector *x = gsl_vector_alloc (n);
  gsl_permutation * p = gsl_permutation_alloc (n);
  //fill matrix/vector:
  for(int i=0;i<n;i++){
    gsl_vector_set(b,i,invec[i]);
    for(int j=0;j<n;j++){
      //gsl_matrix_set(A,i,j,inmat[i*n+j]);
      gsl_matrix_set(A,i,j,inmat[i][j]);
    }
  }
  //peform LU decomposition (using GSL)
  //and solve linear equation A.x=b for x (if non-singular)
  int s;
  gsl_linalg_LU_decomp (A, p, &s);
  double det=gsl_linalg_LU_det (A, s);
  if(det!=0)gsl_linalg_LU_solve (A, p, b, x);
  if(det==0)iRet=1;
  //Fill the output vector:
  for(int i=0;i<n;i++){
    outvec[i]=gsl_vector_get(x,i);
  }
  //clear memory
  gsl_permutation_free (p);
  gsl_matrix_free (A);
  gsl_vector_free (x);
  gsl_vector_free (b);
  return iRet;
}









////******************************************************************************
//int MFS_FFT(double data[], int n, int direction)
///*
// NOT FINISHED  (actually, barely started)
//https://www.gnu.org/software/gsl/manual/html_node/Mixed_002dradix-FFT-routines-for-real-data.html
//*/
//{

//  //??? include option for power of 2 data!

//  //Prepare trigonometric lookup tables for an FFT of size n real elements:
//  //Note: technically, we only need to form this once! Make another function?
//  gsl_fft_real_wavetable * real_wt;
//  real_wt = gsl_fft_real_wavetable_alloc (n);
//
//  //Prepare trigonometric lookup tables for FFT; size n half-complex elements:
//  gsl_fft_halfcomplex_wavetable * hc_wt;
//  hc_wt = gsl_fft_halfcomplex_wavetable_alloc (n);

//  //Allocates a workspace for a real transform of length n:
//  //nb: same workspace can be used for forwards/backwards
//  gsl_fft_real_workspace * work;
//  work = gsl_fft_real_workspace_alloc (n);
//
//  //Compute the FFT of data, a real array of length n
//  int stride = 1;
//  if(direction==1){
//    // Exp( +1 * 2Pi i j k / n ) : Inverse transform!
//    gsl_fft_halfcomplex_inverse (data, 1, n, hc, work);
//  }else{
//    // Exp( -1 * 2Pi i j k / n ) : Fourier transform!
//    gsl_fft_real_transform (data, stride, n, real_wt, work);
//  }
//
//  //Free the memory associated with the workspace:
//  gsl_fft_real_workspace_free (work);
//
//  //Free the memory associated with the wavetable and workspace:
//  gsl_fft_real_wavetable_free (real_wt);
//  gsl_fft_halfcomplex_wavetable_free (hc_wt)
//
//
//}


////******************************************************************************
//int MFS_FFT(double data[], int n, int direction)
//{





//}






//******************************************************************************
double MFS_inverseErf(double x)
/*
170525.
Function uses an approximate method to calculate the inverse error function.
There is no cmath inverseErf function.
x must be in interval (-1,1)
It will return '+/-100' if -1 or 1 is given. ok?
It uses a series expansion for small |x|<0.25.
Two regions for the series expansion, so it doesn't have to call large powers
in the case that x is very small.
Then, uses an algorith I found from:
https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
--I couldn't find the original source
Could be made faster by using an approximation for the log!
*/
{
  if(x==0)return 0;
  if(x>=1)return 100.;   //?? allows for error due to floating point errors
  if(x<=-1)return -100.; //??
  int sgn=1;
  if(x<0)sgn=-1;
  double lnx=log(1.-x*x);
  double tt1=4.33+0.5*lnx;
  double tt2=6.803*lnx;
  return sgn*sqrt(sqrt(tt1*tt1-tt2)-tt1);
}


//******************************************************************************
double MFS_fastInvErf(double x)
/*
170622.
Function uses an approximate method to calculate the inverse error function.
There is no cmath inverseErf function.
x must be in interval (-1,1)
It will return '+/-100' if -1 or 1 is given. ok?
It uses a series expansion for small |x|<0.25.
And for 0.25<|x|<0.95, uses series around non-0 points.
Several regions for the series expansion, so it doesn't have to call large
powers in the case that x is very small etc.
Then, uses an algorith I found from:
https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
--I couldn't find the original source
Could be made faster by using an approximation for the log!
This method is ~8x faster than other for |x|<0.95
Accurate to 1e-4 (almost 1e-5)
*/
{
  if(x<0.01&&x>-0.01)return 0.886227*x;
  if(x<0.25&&x>-0.25)return 0.886227*x-0.2320137*pow(x,3)+0.1275562*pow(x,5);
  if(x>=1)return 100.;   //?? allows for error due to floating point errors
  if(x<=-1)return -100.; //??
  double z=fabs(x);
  int sgn=1;
  if(x<0)sgn=-1;
  if(z<=0.95){
    double w;
    if(z<0.55){//order 7 series around 0.4
      double y=z-0.4;
      w=0.3708072+1.016856*y+0.383413*pow(y,2)+0.543234*pow(y,3)
        +0.571544*pow(y,4)+0.777862*pow(y,5)+1.028110*pow(y,6)
        +1.45202*pow(y,7);
    }else if(z<0.85){//order 8 series around 0.7
      double y=z-0.7;
      w=0.7328691+1.516363*y+1.68513*pow(y,2)+3.65912*pow(y,3)+8.6827*pow(y,4)
        +22.4762*pow(y,5)+60.945*pow(y,6)+170.820*pow(y,7)+490.30*pow(y,8);
    }else{//order 8 series around 0.9
      double y=z-0.9;
      w=1.163087+3.42804*y+13.6680*pow(y,2)+86.089*pow(y,3)+621.95*pow(y,4)
        +4846.6*pow(y,5)+39583.*pow(y,6)+3.3382e5*pow(y,7)+2.8817e6*pow(y,8);
    }
    return sgn*w;
  }
  //if z>0.95, just use "normal" approximation.
  double lnx=log(1.-x*x);
  double tt1=4.33+0.5*lnx;
  double tt2=6.803*lnx;
  return sgn*sqrt(sqrt(tt1*tt1-tt2)-tt1);
}


//******************************************************************************
double MFS_fastErf(double x)
/*
170620.
Up to 10x faster!
Uses a series expansion about 0 for values of |x|<0.35.
Two regions for the series expansion, so it doesn't have to call large powers
in the case that x is very small.
For |x|>0.35, reverts to regular cmath erf function.
Accurate to better than 1e-5
See also:Abramowitz and Stegun (equations 7.1.25–28)
  https://en.wikipedia.org/wiki/Error_function
*/
{
 if(x>3.4)return 1.;
 if(x<-3.4)return -1.;
 if(x<0.005&&x>-0.005)return 1.128379*x;
 if(x<0.2&&x>-0.2)return 1.128379*x-0.3761264*pow(x,3)+0.1128379*pow(x,5);
 if(x<0.35&&x>-0.35)return 1.128379*x-0.3761264*pow(x,3)+0.1128379*pow(x,5)
                           -0.02686617*pow(x,7);
 double z=fabs(x);
 double w=0.;
 int sgn=1;
 if(x<0)sgn=-1;
 if(z<0.7){//order 5 series around 0.5
   double y=z-0.5;
   w=0.5204999+0.8787826*y-0.4393913*pow(y,2)-0.1464638*pow(y,3)
            +0.1830797*pow(y,4)+0.007323188*pow(y,5);
 }else if(z<1.2){//order 5 series around 1
   double y=z-1.;
   w=0.84270+0.415107*y-0.415107*pow(y,2)+0.1383692*pow(y,3)
            +0.0691846*pow(y,4)-0.0691846*pow(y,5);
 }else if(z<2.2){//order 8 series around 1.75
   double y=z-1.75;
   w=0.9866717+0.05277500*y-0.09235624*pow(y,2)+0.09015728*pow(y,3)
            -0.04810221*pow(y,4)+0.006624361*pow(y,5)+0.008963045*pow(y,6)
            -0.006058751*pow(y,7)+0.0007300512*pow(y,8);
 }else if(z<=3.4){//order 7 series around 2.75
   double y=z-2.75;
   w=0.9998994+0.0005862772*y-0.001612262*pow(y,2)+0.002760389*pow(y,3)
            -0.003258114*pow(y,4)+0.002755808*pow(y,5)-0.001657327*pow(y,6)
            +0.0006460410*pow(y,7);
 }
 return sgn*w;
 //return erf(x);
}


//******************************************************************************
double MFS_fastExp(double x)
/*
170622.
Function that calculates exponential function very quickly.
It uses series expansion for values of |x|<~10 (three cases,
so that for smaller x values it doesn't have to evaluate large powers).
For ~10<|x|<400, it uses a method based on:
   e^x = lim n->\infty (1+x/n)^n, with n=2^34
For negative values, calls the inverse
For |x|>400, reverts to regular cmath exp function.
Accurate to better than 1e-5
x<1,      10x faster
1<x<10    3x faster
10<x<400  2.3x faster
And about ~30% slower for negative numbers
Another option is to use look-up tables, with extrapolation!?
*/
{
 if(x<0.03&&x>=0)return 1.+x+0.5*x*x;
 if(x<4.1&&x>0)return 1.+x+0.5*x*x+0.1666666667*pow(x,3)
         +0.04166666667*pow(x,4)+0.008333333333*pow(x,5)
         +0.001388888889*pow(x,6)
         +0.0001984126984*pow(x,7)+0.00002480158730*pow(x,8)
         +2.755731922e-6*pow(x,9)+2.755731922e-7*pow(x,10)
         +2.505210839e-8*pow(x,11)+2.087675699e-9*pow(x,12)
         +1.605904384e-10*pow(x,13)+1.147074560e-11*pow(x,14)
         +7.647163732e-13*pow(x,15);
 if(x<10.5&&x>0)return 1.+x+0.5*x*x+0.1666666667*pow(x,3)
         +0.04166666667*pow(x,4)+0.008333333333*pow(x,5)
         +0.001388888889*pow(x,6)
         +0.0001984126984*pow(x,7)+0.00002480158730*pow(x,8)
         +2.755731922e-6*pow(x,9)+2.755731922e-7*pow(x,10)
         +2.505210839e-8*pow(x,11)+2.087675699e-9*pow(x,12)
         +1.605904384e-10*pow(x,13)+1.147074560e-11*pow(x,14)
         +7.647163732e-13*pow(x,15)+4.779477332e-14*pow(x,16)
         +2.811457254e-15*pow(x,17)+1.561920697e-16*pow(x,18)
         +8.220635247e-18*pow(x,19)+4.110317623e-19*pow(x,20)
         +1.957294106e-20*pow(x,21)+8.896791392e-22*pow(x,22)
         +3.868170171e-23*pow(x,23)+1.611737571e-24*pow(x,24)
         +6.446950284e-26*pow(x,25)+2.479596263e-27*pow(x,26)
         +9.183689864e-29*pow(x,27)+3.279889237e-30*pow(x,28)
         +1.130996289e-31*pow(x,29)+3.769987629e-33*pow(x,30);
 if(x<400.&&x>0){
   double myexp=1.+x*5.820766091e-11;
   for(int i=0;i<34;i++)myexp*=myexp;
   return myexp;
 }
 //above approx works better for +ve args:
 if(x<0&&x>-400.)return 1/MFS_fastExp(-x);
 return exp(x);
}
