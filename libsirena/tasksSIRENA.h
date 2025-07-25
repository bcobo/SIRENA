/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      TASKSSIRENA
* 
*  File:       taskssirena.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef TASKSSIRENA_H
#define TASKSSIRENA_H 1

#include "integraSIRENA.h"
#include "pulseprocess.h"
#include <iomanip>      // std::setprecision

#include "versionSIRENA.h"

#include <time.h>

void runDetect(TesRecord* record, 
               int trig_reclength,
               int lastRecord, 
               int nrecord,
               PulsesCollection *pulsesAll, 
               ReconstructInitSIRENA** reconstruct_init, 
               PulsesCollection** pulsesInRecord);

void th_runDetect(TesRecord* record, int trig_reclength,
                  int lastRecord, 
                  int nrecord,
                  PulsesCollection *pulsesAll, 
                  ReconstructInitSIRENA** reconstruct_init, 
                  PulsesCollection** pulsesInRecord);

int createLibrary(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject);
int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, fitsfile **dtcObject);
int filderLibrary(ReconstructInitSIRENA** reconstruct_init, double samprate);
int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double);
int procRecord(ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, double threshold, int windowSize, int offset, PulsesCollection *foundPulses,long num_previousDetectedPulses, int pixid, gsl_vector *phid, int oscillations, int nrecord, double tstartPrevPulse);
int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject);
int writeTestInfo(ReconstructInitSIRENA* reconstruct_init, gsl_vector *recordDERIVATIVE, double threshold, fitsfile *dtcObject);
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, gsl_vector **pulseaverage_B0, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight);
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhisto, gsl_vector **yhisto);
int align(gsl_vector **vector1, gsl_vector ** vector2);
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m);
int weightMatrix (ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, long nonpileupPulses, gsl_vector *nonpileup, gsl_vector *pulseaverage, gsl_matrix **covariance, gsl_matrix **weight);
int weightMatrixReduced (ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, gsl_matrix *TransformedData, gsl_vector *TransformedPulseAverage , gsl_matrix **covariancer, gsl_matrix **weightr);
int eigenVV (gsl_matrix *matrixin, gsl_matrix **eigenvectors, gsl_vector **eigenvalues);
int writeLibrary(ReconstructInitSIRENA **reconstruct_init, double samprate, double estenergy, gsl_vector *pulsetemplate, gsl_vector *pulsetemplate_B0, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject);
int addFirstRow(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int runF0orB0val, gsl_vector *E, gsl_vector *PHEIGHT, gsl_matrix *PULSE, gsl_matrix *PULSEB0, gsl_matrix *MF, gsl_matrix *MFB0, gsl_matrix *COVAR, gsl_matrix *WEIGHT);
int readAddSortParams(ReconstructInitSIRENA *reconstruct_init,fitsfile **inLibObject,double samprate,int eventcntLib, double estenergy, gsl_vector *pulsetemplate, gsl_vector *pulsetemplate_B0, gsl_matrix *covariance, gsl_matrix *weight);
int calculateIntParams(ReconstructInitSIRENA *reconstruct_init, int indexa, int indexb, double samprate, int runF0orB0val, gsl_matrix *modelsaux,gsl_matrix *covarianceaux, gsl_matrix *weightaux,gsl_vector *energycolumn, gsl_matrix **Wabaux,gsl_matrix **TVaux,gsl_vector **tEcolumn,gsl_matrix **XMaux,gsl_matrix **YVaux,gsl_matrix **ZVaux,gsl_vector **rEcolumn, gsl_matrix **Dabaux, gsl_matrix **Sabaux, gsl_matrix **PrecalCOVaux, gsl_matrix **optimalfiltersabFREQaux, gsl_matrix **optimalfiltersabTIMEaux);
int matrix2vector (gsl_matrix *matrixin, gsl_vector **vectorout);
int vector2matrix (gsl_vector *vectorin, gsl_matrix **matrixout);
int convertI2R (char* EnergyMethod,double Ibias, double Imin, double Imax, double ADU_CNV, double ADU_BIAS, double I_BIAS, double Ifit, double V0, double RL, double L, gsl_vector **invector, int real_data);
int filterByWavelets (ReconstructInitSIRENA* reconstruct_init, gsl_vector **invector, int length, int *onlyOnce);
int obtainRiseFallTimes (gsl_vector *recordNOTFILTERED, double samprate, gsl_vector *tstartgsl, gsl_vector *tendgsl, gsl_vector *Bgsl, gsl_vector *Lbgsl, int numPulses, gsl_vector **tauRisegsl, gsl_vector **tauFallgsl);

void runEnergy(TesRecord* record, int lastRecord, int nrecord, int trig_reclength, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, PulsesCollection *pulsesAll);

void th_runEnergy(TesRecord* record, int nrecord, int trig_reclength,
                  ReconstructInitSIRENA** reconstruct_init,
                  PulsesCollection** pulsesInRecord,
                  PulsesCollection *pulsesAll);

int calculus_optimalFilter(int TorF, int intermediate, int opmode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex);
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out);
int find_matchedfilter(int runF0orB0val, double maxDER, gsl_vector *maxDERs, int preBuffer, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, double *Ealpha, double *Ebeta, double margin);
int find_matchedfilterSAB(double maxDER, gsl_vector *maxDERs, int preBuffer, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, gsl_vector **DabFound, double *Ealpha, double *Ebeta, double margin);
int find_optimalfilter(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, double *Ealpha, double *Ebeta, double margin);
int find_optimalfilterSAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, gsl_vector **DabFound, double *Ealpha, double *Ebeta, double margin);
int find_prclcov(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_matrix **PRCLCOVFound, gsl_vector **DabFound,double *Ealpha, double *Ebeta, double margin);
int find_prclofwn(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_matrix **PRCLOFWNFound,double *Ealpha, double *Ebeta, double margin);
int find_Esboundary(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, int *indexEalpha, int *indexEbeta,double *Ealpha, double *Ebeta, double margin);
int pulseGrading (ReconstructInitSIRENA *reconstruct_init, int tstart, long grade1, long grade2, int *pulseGrade, long *OFlength, int nrecord);
int calculateEnergy (gsl_vector *pulse, gsl_vector *filter, gsl_vector_complex *filterFFT, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, gsl_vector *Dab, gsl_matrix *PRCLCOV, gsl_matrix *PRCLOFWN, double *calculatedEnergy, double *tstartNewDev, int *lagsShift, int LowRes, int productSize, int tooshortPulse_NoLags);
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double energy, gsl_vector *optimalfilter, fitsfile **dtcObject);

#endif /* TASKSSIRENA_H */




