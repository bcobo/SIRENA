/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      TESRECONS
*
*  File:       tesrecons.c
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#include "tesrecons.h"

/***** SECTION 1 ************************************
* MAIN function: This function is mainly a wrapper to pass a data file to the SIRENA tasks in order to reconstruct the energies.
*
* The user must supply the following input parameters (.par file).
* 
* Parameters:
* 
* - RecordFile: Record FITS file
*   - If RecordFile starts with '@' it provides a file text containing several record input FITS files
* - TesEventFile: Output event list file
* - LibraryFile: Library FITS file to be created
* - XMLFile: XML input FITS file with instrument definition
* - EventListSize: Default size of the event list per record
* - clobber:Overwrite or not output files if exist (1/0)
* - history: write program parameters into output file
* - scaleFactor: Detection scale factor for initial filtering
* - samplesUp: Number of consecutive samples up for threshold trespassing (only used with STC detection mode)
* - samplesDown: Number of consecutive samples below the threshold to look for other pulse (only used with STC detection mode)
* - nSgms: Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm
* - detectionMode: Adjusted Derivative (AD) or Single Threshold Crossing (STC)
* - detectSP: Detect secondary pulses (1) or not (0)
* - LbT: Baseline averaging length (seconds)
* - intermediate: Write or not intermediate files (1/0)
* - detectFile: Intermediate detections file (if intermediate=1)
* - FilterDomain: Filtering Domain: Time (T) or Frequency (F)
******* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline)
* - EnergyMethod: Energy calculation Method: OPTFILT, 0PAD, INTCOVAR, COVAR, I2R or I2RFITTED
* - Ifit: Constant to apply the I2RFITTED conversion
* - OFNoise: Noise to use with Optimal Filtering: NSD or WEIGHTN
* - LagsOrNot: Lags or no lags (1/0)
* - nLags: Number of lags (positive odd number)
* - Fitting35: Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
* - OFIter: Iterate or not iterate (1/0)
* - OFLib: Work or not with a library with optimal filters (yes/no)
* - OFStrategy: Optimal Filter length Strategy: **FREE**, **BYGRADE** or **FIXED**
* - OFLength: Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**)
* - flength_0pad: 0-padding filter length
* - prebuff_0pad: preBuffer when 0-padding
* - int errorT: Additional error (in samples) added to the detected time (Logically, it changes the reconstructed energies)
* - Sum0Filt: If 0-padding, subtract the sum of the filter (1) or not (0)
* - tstartPulse1: Integer number: Sample where the first pulse starts or nameFile: File where the tstart (seconds) of every pulse is
* - tstartPulse2: Tstart (samples) of the second pulse
* - tstartPulse3: Tstart (samples) of the third pulse (if 0 => PAIRS, if not 0 => TRIOS)
* - energyPCA1: First energy (only for PCA)
* - energyPCA2: Second energy (only for PCA)
* 
* Steps:
* 
* - Register HEATOOL
* - Reading all programm parameters by using PIL
* - Read XML info
* - getSamplingrate_trigreclength => Obtain the 'trig_reclength' and the sampling rate
* - Sixt standard keywords structure
* - Open output FITS file
* - Initialize data structures needed for pulse filtering
* - Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
* - Build up TesEventList
* - Call SIRENA to build reconstruct the energies
* - Save GTI extension to event file
* - Free memory
*****************************************************/
int tesrecons_main() {
  printf("Running TESRECONS v%s\n",SIRENA_VERSION);
  time_t ttstart = time(0);
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  par.OFLib = 1;        // Debugger complains about an initialized variable (only the boolean type)
  
  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("tesrecons");
  
  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");

    // Get program parameters
    status=getpar_tesrecons(&par);
    if (status != EXIT_SUCCESS)
    {
    	return(status);
    }
    par.opmode = 1; // Reconstructing the energies

    if ((strcmp(par.EnergyMethod,"I2RFITTED") == 0) && (par.Ifit == 0.0))
    {
        SIXT_ERROR("Ifit value must be provided");
        return(EXIT_FAILURE);
    }
        
    AdvDet *det = newAdvDet(&status);
    CHECK_STATUS_BREAK(status);
    double sampling_rate_XML = -999.; // Sampling rate from XML
    // Read XML info
    det = loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    sampling_rate_XML = det->SampleFreq;
    if ((strcmp(par.EnergyMethod,"0PAD") == 0) && (strcmp(par.OFStrategy,"FIXED") != 0))
    {
        printf("%s","Attention: EnergyMethod=0PAD => OFStrategy set to FIXED\n");
        strcpy(par.OFStrategy,"FIXED");
    }

    // Obtain the 'trig_reclength' and the sampling rate
    double sampling_rate = -999.0;
    int trig_reclength = -999;
    sampling_rate = sampling_rate_XML;
    int numfits = -999;
    status = getSamplingrate_trigreclength (par.RecordFile,par,&sampling_rate,&trig_reclength, &numfits);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Error in 'getSamplingrate_trigreclength' function");
        return(EXIT_FAILURE);
    }

    // Sixt standard keywords structure
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);

    //Open input records file to get the PH_ID dimension
    //(to dimension output PH_ID column in output events file)
    fitsfile* fptr = NULL;
    fits_open_file(&fptr, par.RecordFile, READONLY, &status);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Error opening input RecordFile");
        return(EXIT_FAILURE);
    }
    fits_movnam_hdu(fptr, ANY_HDU, "TESRECORDS", 0, &status);
    if (status != 0) {
        status = 0;
        fits_movnam_hdu(fptr, ANY_HDU, "RECORDS", 0, &status);
        if (status != EXIT_SUCCESS)
        {
            SIXT_ERROR("Error: Neither 'TESRECORDS' nor 'RECORDS' HDUs found in RecordFile");
            return(EXIT_FAILURE);
        }
    }
    int colnum = 0;
    char column_name[] = "PH_ID";
    fits_get_colnum(fptr, CASEINSEN, column_name, &colnum, &status);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Cannot find PH_ID column in RecordFile");
        return(EXIT_FAILURE);
    }
    int naxis = 0;
    long ph_id_column_dim = 0;
    fits_read_tdim(fptr, colnum, 1, &naxis, &ph_id_column_dim, &status);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Error running fits_read_tdim in tesrecons");
        return(EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status != EXIT_SUCCESS)
    {
        SIXT_ERROR("Error closing input records file (RecordFile)");
        return(EXIT_FAILURE);
    }
    
    //Open outfile (events)
    TesEventFileSIRENA* outfile = opennewTesEventFileSIRENA(par.TesEventFile,
                                                 keywords,
                                                 SIRENA_VERSION,
                                                 ph_id_column_dim,
                                                 par.clobber,
                                                 &status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize data structures needed
    ReconstructInitSIRENA* reconstruct_init_sirena = newReconstructInitSIRENA();
    CHECK_STATUS_BREAK(status);
    PulsesCollection* pulsesAll = newPulsesCollection();
    CHECK_STATUS_BREAK(status);  
    OptimalFilterSIRENA* optimalFilter = newOptimalFilterSIRENA();
    CHECK_STATUS_BREAK(status);// define a second structure for calibration
    
    // Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
    status = fillReconstructInitSIRENAGrading (par, det, &reconstruct_init_sirena);
    destroyAdvDet(&det);

    // Build up TesEventList
    /*TesEventListSIRENA* event_list = newTesEventListSIRENA(&status);
    allocateTesEventListTriggerSIRENA(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);*/
            
    // Call SIRENA to reconstruct
    status = callSIRENA(par.RecordFile, keywords, reconstruct_init_sirena, par, sampling_rate, &trig_reclength, pulsesAll, outfile, ph_id_column_dim);

    // Save GTI extension to event file
    GTI* gti=getGTIFromFileOrContinuous("none",keywords->tstart, keywords->tstop,keywords->mjdref, &status);
    saveGTIExt(outfile->fptr, "STDGTI", gti, &status);    
    CHECK_STATUS_BREAK(status);
    
    //Free memory
    if (reconstruct_init_sirena->grading != NULL)
    {
        gsl_vector_free(reconstruct_init_sirena->grading->value); reconstruct_init_sirena->grading->value = 0;
        gsl_matrix_free(reconstruct_init_sirena->grading->gradeData); reconstruct_init_sirena->grading->gradeData = 0;
    }
    free(reconstruct_init_sirena->grading);
    reconstruct_init_sirena->grading = 0;
    freeReconstructInitSIRENA(reconstruct_init_sirena);
    freePulsesCollection(pulsesAll);
    freeOptimalFilterSIRENA(optimalFilter);
    freeTesEventFileSIRENA(outfile,&status);
    //freeTesEventListSIRENA(event_list);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);
 
  } while(0); // END of the error handling loop.
  
  if (EXIT_SUCCESS==status) 
  {
      headas_chat(3, "finished successfully!\n\n");
      time_t ttcurrent = time(0);
      printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
      return(EXIT_SUCCESS);
  } 
  else
  {
      //time_t ttcurrent = time(0);
      //printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
      return(status);
  }
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* getpar_tesrecons function: This function gets the input parameter from the command line or their default values from the tesrecons.par file
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int getpar_tesrecons(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_string("RecordFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the optimal filter file");
    return(status);
  }
  strcpy(par->RecordFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("TesEventFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  }
  strcpy(par->TesEventFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("LibraryFile", &sbuffer);
  strcpy(par->LibraryFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  }
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("EventListSize", &par->EventListSize);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the EventListSize parameter");
	  return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the clobber parameter");
	  return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
  status=ape_trad_query_int("samplesUp", &par->samplesUp);
  status=ape_trad_query_int("samplesDown", &par->samplesDown);
  status=ape_trad_query_double("nSgms", &par->nSgms);

  status=ape_trad_query_string("detectionMode", &sbuffer);
  strcpy(par->detectionMode, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("detectSP", &par->detectSP);
      
  status=ape_trad_query_double("LbT", &par->LbT);

  status=ape_trad_query_int("intermediate", &par->intermediate);
  status=ape_trad_query_string("detectFile", &sbuffer);
  strcpy(par->detectFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("FilterDomain", &sbuffer);
  strcpy(par->FilterDomain, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("FilterMethod", &sbuffer);
  strcpy(par->FilterMethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("EnergyMethod", &sbuffer);
  strcpy(par->EnergyMethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("filtEev", &par->filtEev);
  assert(par->flength_0pad >= 0);

  status=ape_trad_query_double("Ifit", &par->Ifit);

  status=ape_trad_query_string("OFNoise", &sbuffer);
  strcpy(par->OFNoise, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("LagsOrNot", &par->LagsOrNot);
  status=ape_trad_query_int("nLags", &par->nLags);

  status=ape_trad_query_int("Fitting35", &par->Fitting35);
  status=ape_trad_query_int("OFIter", &par->OFIter);

  status=ape_trad_query_bool("OFLib", &par->OFLib);

  strcpy(par->OFInterp, "SAB");

  status=ape_trad_query_string("OFStrategy", &sbuffer);
  strcpy(par->OFStrategy, sbuffer);
  free(sbuffer);

  if ((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(par->EnergyMethod,"I2RDER") == 0))
  {
    if (strcmp(par->OFStrategy,"FREE") == 0)  par->OFLib = 0;
    else if (strcmp(par->OFStrategy,"FIXED") == 0)    par->OFLib = 1;
    else if (strcmp(par->OFStrategy,"BYGRADE") == 0)  par->OFLib = 1;
  }

  status=ape_trad_query_int("OFLength", &par->OFLength);

  status=ape_trad_query_int("flength_0pad", &par->flength_0pad);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the flength_0pad parameter");
    return(status);
  }
  assert(par->flength_0pad > 0);

  status=ape_trad_query_int("prebuff_0pad", &par->prebuff_0pad);
  assert(par->prebuff_0pad >= 0);

  status=ape_trad_query_int("errorT", &par->errorT);
  status=ape_trad_query_int("Sum0Filt", &par->Sum0Filt);
  status=ape_trad_query_string("tstartPulse1", &sbuffer);
  strcpy(par->tstartPulse1, sbuffer);
  free(sbuffer);
  int isNumber = 1;
  for (int i = 0; i < (int)(strlen(par->tstartPulse1)); i++)
  {
      if (isdigit(par->tstartPulse1[i]) == 0)
      {
          isNumber = 0;
          break;
      }
  }
  status=ape_trad_query_int("tstartPulse2", &par->tstartPulse2);
  status=ape_trad_query_int("tstartPulse3", &par->tstartPulse3);
  par->energyPCA1 = -999;
  par->energyPCA1 = -999;
      
  if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading some TESRECONS parameter");
      return(status);
  }

  MyAssert(par->LbT > 0, "LbT must be greater than 0");

  MyAssert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0), "FilterDomain must be T or F");

  //MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0) || (strcmp(par->FilterMethod,"F0B0") == 0),"FilterMethod must be F0 or B0 or F0B0");
  MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");

  //MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"0PAD") == 0) || (strcmp(par->EnergyMethod,"INTCOVAR") == 0) || (strcmp(par->EnergyMethod,"COVAR") == 0) ||
  //(strcmp(par->EnergyMethod,"I2R") == 0) ||	(strcmp(par->EnergyMethod,"I2RFITTED") == 0) ||	(strcmp(par->EnergyMethod,"PCA") == 0)), "EnergyMethod must be OPTFILT, 0PAD, INTCOVAR, COVAR, I2R, I2RFITTED or PCA");
  MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"0PAD") == 0) || (strcmp(par->EnergyMethod,"INTCOVAR") == 0) || (strcmp(par->EnergyMethod,"COVAR") == 0) ||
  (strcmp(par->EnergyMethod,"I2R") == 0) ||	(strcmp(par->EnergyMethod,"I2RFITTED") == 0), "EnergyMethod must be OPTFILT, 0PAD, INTCOVAR, COVAR, I2R or I2RFITTED");

  if ((isNumber == 0) && (strcmp(par->FilterDomain,"F") == 0))    // It is only implemented tstartPulse1 as a file for time domain
  {
      SIXT_ERROR("It is not possible to work in FREQUENCY domain if tstartPulse1 is a file => Change FilterDomain to TIME domain (T) ");
      return(EXIT_FAILURE);
  }

  MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");

  MyAssert((strcmp(par->OFNoise,"NSD") == 0) || (strcmp(par->OFNoise,"WEIGHTN") == 0), "OFNoise must be NSD or WEIGHTN");

  MyAssert((strcmp(par->detectionMode,"AD") == 0) || (strcmp(par->detectionMode,"STC") == 0), "detectionMode must be AD or STC");

  MyAssert((par->LagsOrNot ==0) || (par->LagsOrNot ==1), "LagsOrNot must me 0 or 1");
  if ((par->nLags)%2 == 0)
  {
      SIXT_ERROR("parameter error: nLags must be odd");
      return(EXIT_FAILURE);
  }
  MyAssert((par->Fitting35 ==3) || (par->Fitting35 ==5), "Fitting35 must me 3 or 5");
  if ((par->Fitting35 ==3) && (par->nLags<3))
  {
      SIXT_ERROR("parameter error: nLags must be at least 3");
      return(EXIT_FAILURE);
  }
  if ((par->Fitting35 ==5) && (par->nLags<5))
  {
      SIXT_ERROR("parameter error: nLags must be at least 5");
      return(EXIT_FAILURE);
  }

  MyAssert((par->Sum0Filt ==0) || (par->Sum0Filt ==1), "Sum0Filt must be 0 or 1");

  if ((strcmp(par->EnergyMethod,"INTCOVAR") == 0) && (par->LagsOrNot == 1))
  {
      SIXT_ERROR("parameter error: EnergyMethod=INTCOVAR and Lags not implemented yet");
      return(EXIT_FAILURE);
  }

  MyAssert((par->OFIter ==0) || (par->OFIter ==1), "OFIter must be 0 or 1");

  if ((strcmp(par->EnergyMethod,"0PAD") == 0) && (strcmp(par->FilterDomain,"F") == 0))
  {
      SIXT_ERROR("parameter error: Code is not prepared to run 0-padding in Frequency domain");
      // To run 0-padding in Frequency domain the steps should be:
      //1. Take the 8192-samples-length filter in Time domain
      //2. Cut the 0-padding length first samples (the first 4096 samples, or the first 2048 samples...) => 0-padding filter
      //3. FFT of the 0-padding filter
      //4. FFT of the 0-padding pulse (pulse cut according the 0-padding)
      //5. Scalar product in Frequency domain
      return(EXIT_FAILURE);
  }
  if ((strcmp(par->EnergyMethod,"0PAD") == 0) && (strcmp(par->OFNoise,"WEIGHTN") == 0))
  {
      SIXT_ERROR("parameter error: EnergyMethod=0PAD && OFNoise=WEIGHTN not a valid choice => OFNoise should be NSD");
      return(EXIT_FAILURE);
  }
  if ((strcmp(par->EnergyMethod,"0PAD") == 0) && (par->flength_0pad > par->OFLength))
  {
      SIXT_ERROR("parameter error: EnergyMethod=0PAD => flength_0pad should be less than or equal to the maximum filter length (OFLength)");
      return(EXIT_FAILURE);
  }

  if ((strcmp(par->EnergyMethod,"INTCOVAR") == 0) && (par->OFLib == 0))
  {
      SIXT_ERROR("parameter error: EnergyMethod=INTCOVAR => OFLib should be 'yes'");
      return(EXIT_FAILURE);
  }

  if ((strcmp(par->EnergyMethod,"OPTFILT") == 0) && (strcmp(par->OFNoise,"WEIGHTN") == 0) && (par->OFLib == 0))
  {
      SIXT_ERROR("parameter error: EnergyMethod=OPTFILT && OFNoise=WEIGHTN => OFLib should be 'yes'");
      return(EXIT_FAILURE);
  }

  MyAssert((strcmp(par->OFStrategy,"FREE") == 0) || (strcmp(par->OFStrategy,"BYGRADE") == 0) || (strcmp(par->OFStrategy,"FIXED") == 0),
           "OFStrategy must be FREE, BYGRADE or FIXED");

  MyAssert(par->OFLength > 0, "OFLength must be greater than 0");

  return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
