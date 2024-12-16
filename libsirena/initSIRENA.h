/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      INITSIRENA
*
*  File:       initSIRENA.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*
***********************************************************************/

#ifndef INITSIRENA_H
#define INITSIRENA_H 1

#include <stdio.h>
#include <stdlib.h>

#include <fitsio.h>

#include "sixt.h"
#include "advdet.h"
#include "testriggerfile.h"
#include "testrigger.h"

#include "integraSIRENA.h"

struct Parameters {
	// File containing the optimal filter
	char OptimalFilterFile[MAXFILENAME];

	// File to reconstruct
	char RecordFile[MAXFILENAME];

	// Ouput event list
	char TesEventFile[MAXFILENAME];

	// File containing the pulse template
	char PulseTemplateFile[MAXFILENAME];

	// 0-padding filter length (only necessary when reconstructing with 0-padding)
	int flength_0pad;
	// 0-padding preBuffer (only necessary when reconstructing with 0-padding)
	int prebuff_0pad;

	// Threshold level
	double Threshold;

	// Calibration factor
	double Calfac;

	// Default size of the event list per record
	int EventListSize;

	// Minimal distance before using OFs after a mireconstruction
	int NormalExclusion;

	// Minimal distance before reconstructing any event after a mireconstruction
	int DerivateExclusion;

	// Saturation level of the ADC curves
	double SaturationValue;

	// Boolean to choose whether to erase an already existing event list
	char clobber;

	// Boolean to choose to save the run parameters in the output file
	char history;

	// File containing the library
	char LibraryFile[MAXFILENAME];

	// Scale Factor for initial filtering
	double scaleFactor;

	// Number of samples for threshold trespassing
	int samplesUp;

    // Number of samples below the threshold to look for other pulse
	int samplesDown;

	// Number of standard deviations in the kappa-clipping process for threshold estimation
	double nSgms;

    // Detect secondary pulses or not
    int detectSP;

	// Calibration run (0) or energy reconstruction run (1)?
	int opmode;

    // DetectionMode: Adjusted Derivative(AD) or Single Threshold Crossing(STC)
	char detectionMode[4];

	// Monochromatic energy for library creation
	double monoenergy;

	// Boolean to add or not pre-calculated values related to COVAR reconstruction method in the library file
	char addCOVAR;

	// Boolean to add or not pre-calculated values related to INT_COVAR reconstruction method in the library file
	char addINTCOVAR;

	// Boolean to add or not pre-calculated values related to Optimal Filtering by using Weight Noise matrix in the library file
	char addOFWN;

	// Running sum length for the RS raw energy estimation, in seconds (only in CALIBRATION)
	double LrsT;

	// Baseline averaging length for the RS raw energy estimation, in seconds (only in CALIBRATION)
	double LbT;

	// Noise filename
	char NoiseFile[MAXFILENAME];

	// Filtering Domain: Time(T) or Frequency(F)
	char FilterDomain[2];

	// Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline) of F0B0 (deleting always the baseline)
	char FilterMethod[5];

	// Energy Method: OPTFILT, 0PAD, INTCOVAR, COVAR, I2R, I2RFITTED or PCA
	char EnergyMethod[10];

    // Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R and I2RFITTED)
    double filtEev;

    // Constant to apply the I2RFITTED conversion
    double Ifit;

	// Noise to use with Optimal Filtering: NSD (Noise Spectral Density) or WEIGHTN (weight matrix)
	char OFNoise[8];

	// LagsOrNot: LAGS == 1 or NOLAGS == 0
	int LagsOrNot;

    // Number of lags (odd number)
	int nLags;

    // Using 3 lags to analytically calculate a parabola (3) or using 5 lags to fit (5)
	int Fitting35;

	// OFIter: Iterate == 1 or NOTIterate == 0
	int OFIter;

	// Boolean to choose whether to use a library with optimal filters or calculate the optimal filter to each pulse
	char OFLib;

	// Optimal Filter by using the Matched Filter (MF) or the DAB as matched filter (MF, DAB)
    char OFInterp[4];

	// Optimal Filter length Strategy: FREE, BASE2, BYGRADE or FIXED
	char OFStrategy[8];

	// Optimal Filter length (taken into account if OFStrategy=FIXED)
	int OFLength;

	// Write intermediate files (Yes:1, No:0)
	int intermediate;

	// File with the output detections
	char detectFile[256];

	// File with the output filter (only in calibration)
	char filterFile[256];

    // Additional error (in samples) added to the detected time"  (Logically, it changes the reconstructed energies)
	int errorT;

    // Sum0Filt: 0-padding: Subtract the sum of the filter (1) or not (0)
	int Sum0Filt;

	// Tstart of the pulses (to be used instead of calculating them if tstartPulse1 =! 0)
	//int tstartPulse1;
    char tstartPulse1[MAXFILENAME]; // Integer number: Sample where the first pulse starts
                                    // or
                                    // nameFile: File where is the tstart (in seconds) of every pulse
	int tstartPulse2;
	int tstartPulse3;

	// Energies for PCA
	double energyPCA1;
	double energyPCA2;

	// XML file with instrument definition
	char XMLFile[MAXFILENAME];
};

typedef struct {
	/** Pointer to the FITS file. */
	fitsfile* fptr;

	/** Number of the current row in the FITS file. The numbering
	starts at 1 for the first line. If row is equal to 0, no row
	has been read or written so far. */
	long row;

	/** Total number of rows */
	long nrows;

	/** Column numbers for time, energy, grade1, grade2, pixID, RA and DEC columns */
	int timeCol,energyCol,avg_4samplesDerivativeCol,E_lowresCol,grade1Col,grade2Col,phiCol,lagsShiftCol,bslnCol,rmsbslnCol,pixIDCol,riseCol,fallCol,phIDCol,raCol,decCol,detxCol,detyCol,gradingCol,srcIDCol,nxtCol,extCol; //SIRENA

} TesEventFileSIRENA;

TesEventListSIRENA* newTesEventListSIRENA(int* const status);
void freeTesEventListSIRENA(TesEventListSIRENA* event_list);

TesEventFileSIRENA* newTesEventFileSIRENA(int* const status);

/** Create and open a new TesEventFile. */
TesEventFileSIRENA* opennewTesEventFileSIRENA(const char* const filename,
				  SixtStdKeywords* keywords,
			      const char* const sirenaVersion,
				  const char clobber,
				  int* const status);

void saveEventListToFileSIRENA(TesEventFileSIRENA* file,TesEventListSIRENA * event_list,
	double start_time,double delta_t,int* const status);

/** Allocates memory for a TesEventList structure for the triggering stage:
 *  only event_index, pulse_height and grades1 */
void allocateTesEventListTriggerSIRENA(TesEventListSIRENA* event_list,int size,int* const status);

void freeTesEventFileSIRENA(TesEventFileSIRENA* file, int* const status);

void MyAssert(int expr, char* msg);

int checkXmls(struct Parameters* const par);

char* subString (const char* input, int offset, int len, char* dest);

int getSamplingrate_trigreclength (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength, int* numfits);
int getSamplingrate_trigreclength_Filei (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength);

int fillReconstructInitSIRENAGrading (struct Parameters par, AdvDet *det, ReconstructInitSIRENA** reconstruct_init_sirena);

int callSIRENA_Filei(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFileSIRENA* outfile);
int callSIRENA(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFileSIRENA* outfile);

int checkpreBuffer(struct Parameters* const par);

#endif /* INITSIRENA_H */
