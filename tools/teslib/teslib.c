/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      TESLIB
*
*  File:       teslib.c
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#include "teslib.h"

/***** SECTION 1 ************************************
* MAIN function: This function is mainly a wrapper to pass a data file to the SIRENA tasks in order to build the library that will be used to reconstruct the energies.
*
* The user must supply the following input parameters (.par file).
* 
* Parameters:
* 
* - RecordFile: Record FITS file
*   - If RecordFile starts with '@' it provides a file text containing several record input FITS files
* - TesEventFile: Output event list file
* - LibraryFile: Library FITS file to be created
* - NoiseFile: Noise FITS file with noise spectrum
* - XMLFile: XML input FITS file with instrument definition
* - EventListSize: Default size of the event list per record
* - clobber:Overwrite or not output files if exist (1/0)
* - history: write program parameters into output file
* - threshold: Threshold to use with the derivative to detect (if -999 it is going to be calculated from noise)
* - nSgms: Number of quiescent-signal standard deviations to establish the threshold (if -999 it is going to be calculated from noise)
* - scaleFactor: Detection scale factor for initial filtering
* - samplesUp: Number of consecutive samples up for threshold trespassing
* - windowSize: Window size used to compute the averaged derivative
* - offset: Window offset
* - LrsT: Running sum length for the RS raw energy estimation (seconds) 
* - LbT: Baseline averaging length (seconds)
* - monoenergy: Monochromatic energy of the pulses in the input FITS file in eV 
* - addCOVAR: Add or not pre-calculated values related to COVAR reconstruction method in the library file (1/0)
* - addINTCOVAR: Add or not pre-calculated values related to INTCOVAR reconstruction method in the library file (1/0)
* - addOFWN: Add or not pre-calculated values related to Optimal Filtering by using Weight Noise matrix in the library file (1/0)
* - FilterDomain: Filtering Domain: Time (T) or Frequency (F)
****** - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline) or F0B0 (deleting always the baseline)
* - FilterMethod: Filtering Method: F0 (deleting the zero frequency bin) or B0 (deleting the baseline)
* - EnergyMethod: Energy calculation Method: OPTFILT, INTCOVAR, COVAR, I2R, I2RALL, I2RNOL, I2RFITTED or PCA
* - Ifit: Constant to apply the I2RFITTED conversion
* - intermediate: Write or not intermediate files (1/0)
* - detectFile: Intermediate detections file (if intermediate=1)
* 
* Steps:
* 
* - Register HEATOOL
* - Reading all programm parameters by using PIL
* - Read XML info
* - getSamplingrate_trigreclength => Obtain the 'trig_reclength' and the sampling rate
* - Sixt standard keywords structure
* - Open output FITS file
* - Initialize data structures
* - Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
* - Build up TesEventList
* - Call SIRENA to build the library
* - Save GTI extension to event file
* - Free memory
*****************************************************/
int teslib_main() {
  printf("Running TESLIB v%s\n",SIRENA_VERSION);
  time_t ttstart = time(0);

  // Containing all programm parameters read by PIL.
  struct Parameters par;
  par.addOFWN = 0;  // Debugger complains about an initialized variable (only the boolean type)
  par.addCOVAR = 0;  // Debugger complains about an initialized variable (only the boolean type)
  par.addINTCOVAR = 0;  // Debugger complains about an initialized variable (only the boolean type)
  par.OFLib = 1;        // Debugger complains about an initialized variable (only the boolean type)
  
  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("teslib");
  //set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");

    // Get program parameters
    status = getpar_teslib(&par);
    //CHECK_STATUS_BREAK(status);
    if (status != EXIT_SUCCESS)
    {
    	return(status);
    }
    par.opmode = 0; // Building the library

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
    TesEventFileSIRENA * outfile = opennewTesEventFileSIRENA(par.TesEventFile,
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
            
    // Call SIRENA to build the library
    status = callSIRENA(par.RecordFile, keywords, reconstruct_init_sirena, par, sampling_rate, &trig_reclength, pulsesAll, outfile, ph_id_column_dim);
    CHECK_STATUS_BREAK(status);

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
* getpar_teslib function: This function gets the input parameter from the command line or their default values from the teslib.par file
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int getpar_teslib(struct Parameters* const par)
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

  status=ape_trad_query_string("NoiseFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the noise file");
    return(status);
  }
  strcpy(par->NoiseFile, sbuffer);
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

  status=ape_trad_query_double("threshold", &par->threshold);
  status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
  status=ape_trad_query_int("samplesUp", &par->samplesUp);
  status=ape_trad_query_double("nSgms", &par->nSgms);

  status=ape_trad_query_int("windowSize", &par->windowSize);
  status=ape_trad_query_int("offset", &par->offset);

  status=ape_trad_query_double("LrsT", &par->LrsT);
  status=ape_trad_query_double("LbT", &par->LbT);

  status=ape_trad_query_double("monoenergy", &par->monoenergy);

  status=ape_trad_query_bool("addCOVAR", &par->addCOVAR);
  status=ape_trad_query_bool("addINTCOVAR", &par->addINTCOVAR);
  status=ape_trad_query_bool("addOFWN", &par->addOFWN);

  status=ape_trad_query_string("EnergyMethod", &sbuffer);
  strcpy(par->EnergyMethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_double("Ifit", &par->Ifit);

  status=ape_trad_query_string("FilterMethod", &sbuffer);
  strcpy(par->FilterMethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("intermediate", &par->intermediate);
  status=ape_trad_query_string("detectFile", &sbuffer);
  strcpy(par->detectFile, sbuffer);
  free(sbuffer);
      
  if (EXIT_SUCCESS!=status) {
      SIXT_ERROR("failed reading some TESLIB parameter");
      return(status);
  }

  MyAssert(par->LbT > 0, "LbT must be greater than 0");

  MyAssert(par->monoenergy > 0, "monoenergy must be greater than 0");

  //MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0) || (strcmp(par->FilterMethod,"F0B0") == 0),"FilterMethod must be F0 or B0 or F0B0");
  MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");

  MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"0PAD") == 0) || (strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RFITTED") == 0),
           "EnergyMethod must be OPTFILT, 0PAD, I2R or I2RFITTED");

  if (((strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RFITTED") == 0)) && (par->addCOVAR || par->addINTCOVAR))
  {
    MyAssert(0,"EnergyMethod=I2R/I2RFITTED incompatible with addCOVAR/addINTCOVAR=yes");
  }

  MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");

  MyAssert(par->windowSize >= 0, "windowSize must be greater or equal than 0");
  if (par->windowSize == 0)
  {
    par->offset=0;
    MyAssert(par->offset >= 0, "windowSize=0 => offset=0");
  }
  else
  {
    MyAssert(par->offset >= 0, "offset must be greater or equal than 0");
  }

  return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
