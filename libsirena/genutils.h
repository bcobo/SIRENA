/***********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      GENUTILS
*
*  File:       genutils.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef GENUTILS_H_
#define GENUTILS_H_

// GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <fftw3.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_errno.h>

// General
#include <fitsio.h>
#include <math.h>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <complex>
#include <getopt.h> // For getopt module
#include <stdarg.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <fstream>
using std::ifstream;
#include <sstream>
#include <ctype.h>
#include <sys/stat.h>
#include "assert.h"

#ifndef EPOK	//Event Processing OK
#define EPOK                 (0)
#endif

#ifndef EPFAIL	//Event Processing Failure
#define EPFAIL                 (1)
#endif

#define EP_EXIT_ERROR(msg,status) (exit_error(__func__, msg, status))

#define EP_PRINT_ERROR(msg,status) (print_error(__func__, msg, status))

using namespace std;

const double pi = 4.0 * atan(1.0);

int polyFit (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b, double *c);
int polyFitLinear (gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b);
int FFT(gsl_vector *invector,gsl_vector_complex *outvector);
int FFTinverse(gsl_vector_complex *invector,gsl_vector *outvector);

//GSL vectors
void gsl_vector_sqrtIFCA(gsl_vector *cvnew,gsl_vector *cv);
void gsl_vector_complex_absIFCA(gsl_vector *cvnew,gsl_vector_complex *cv);
void gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z);
void gsl_vector_complex_argIFCA(gsl_vector *varg,gsl_vector_complex *vin);
int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum);

void print_error( const char* const func, string message, int status);
void writeLog (FILE *fileRef, string type, int verbosity, string message);
void exit_error(const char* const func, string msg,int status);
bool fileExists (const std::string& name);
	
int parabola3Pts (gsl_vector *x, gsl_vector *y, double *a, double *b, double *c);
        
bool isNumber(string s);
    
int hannWindow(gsl_vector **inoutvector);

#endif /*GENUTILS_H_*/
