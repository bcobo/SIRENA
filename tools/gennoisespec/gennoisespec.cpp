/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      GENNOISESPEC
*
*  File:       gennoisespec.cpp
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/
 
 /******************************************************************************
  * DESCRIPTION:
  * 
  * The purpose of this package is the calculation the current noise spectral density.
  * 
  * MAP OF SECTIONS IN THIS FILE:
  * 
  * - 1. gennoisespec_main
  * - 2. inDataIterator
  * - 3. findInterval
  * - 4. findIntervalN
  * - 5. createTPSreprFile
  * - 6. writeTPSreprExten
  * - 7. findPulsesNoise
  * - 8. findTstartNoise
  * - 9. weightMatrixNoise
  * - 10. medianKappaClipping_noiseSigma
  * - 11. getpar_noiseSpec
  * - 12. MyAssert
  * 
  *******************************************************************************/
 
 #include "gennoisespec.hpp"
 #include <tasksSIRENA.h>
 
 #include "versionSIRENA.h"
 
 
 /***** SECTION 1 ************************************
  * MAIN function: This function calculates the current noise spectral density.
  *                If there are pulses in a record, the pulses are rejected and it is going to look for pulse-free intervals of a given size (intervalMinBins).
  *                If there are no pulses in a record, the event is divided into pulse-free intervals of a given size (intervalMinBins).
  *                It is going to look for pulse-free intervals, calculate their FFT(not filtered data) and average them.
  * 		       Another feature is calculate the weight matrix of the noise (in fact, weight matrices of the noise of different lengths).
  * 
  * The output FITS file (_noisespec) contains three columns in two extensions, NOISE and NOISEALL:
  * 	- Frequency
  *  	- Current noise spectral density: Amount of current per unit (density) of frequency (spectral), as a function of the frequency
  * 	- Standard error of the mean (filled out with 0's at the moment)
  * 
  * There is also other extension, WEIGHTMS, where the weight matrices of the noise are stored.
  * 	
  * The user must supply the following input parameters (.par file):
  *
  * - inFile: Name of the input FITS file
  * - outFile: Name of the output FITS file
  * - intervalMinSamples: Minimum length of a pulse-free interval to use (samples) = intervalMinBins
  * - nplPF: Number of pulse lengths after ending the pulse (Tend) to start the pulse-free interval
  * - nintervals: Number of pulse-free intervals to use to calculate the Noise Spectral Density
  * - scaleFactor: Scale factor to apply in order to calculate the LPF box-car length
  * - samplesUp: Consecutive samples over the threshold to locate a pulse
  * - nSgms: Number of Sigmas to establish the threshold
  * - weightMS: Calculate and write the weight matrices if weightMS=yes
  * - EnergyMethod: Transform to resistance space (I2R or I2RFITTED) or not (OPTFILT)
  * - clobber: Re-write output files if clobber=yes
  * - rmNoiseIntervals: Remove some noise intervals before calculating the noise spectrum if rmNoiseIntervals=yes
  * 
  * Steps:
  * 
  * - Reading all programm parameters by using PIL
  * - Open input FITS file
  * - Check if input FITS file have been simulated with TESSIM or XIFUSIM
  * - To calculate 'aducnv' (conversion factor between arbitrary units and A)...
  * -...or read ADU_CNV, I_BIAS and ADU_BIAS
  * - Get structure of input FITS file columns
  * - Read info to transform to resistance space (if it is necessary)
  * - Read and check other input keywords
  * - Read other necessary keywords from ANY HDU
  * - Calculate the sampling rate 
  *     - By using keywords in input FITS file (from DELTAT or TCLOCK+DEC_FAC or NUMROW+P_ROW)
  *     - If necessary read the sampling rate from input FITS file (from the HISTORY in the Primary HDU)
  *     - If not possible, provide an error message to include DELTAT (inverse of sampling rate) in the input FITS file
  * - Initialize variables and transform from seconds to samples
  * - Declare variables
  * - Create structure to run Iteration: inDataIterator
  *   - Read columns (TIME and ADC)
  * - Called iteration function: inDataIterator
  * - Close input FITS file
  * - Generate CSD representation
  *   - Applying medianKappaClipping_noiseSigma in order to remove the noise intervals with a high sigma (if `rmNoiseInterval`= yes)
  *   - FFT calculus (EventSamplesFFT)
  *   - Add to mean FFT samples
  *   - Current noise spectral density
  *   - Extra normalization (further than the FFT normalization factor,1/n) in order to get the apropriate noise level provided by Peille (54 pA/rHz)
  * - Load in noiseIntervals only those intervals with a proper sigma and NumMeanSamples = cnt (in order not to change excesively the code when weightMS) 
  * - Generate WEIGHT representation
  * - Create output FITS File: GENNOISESPEC representation file (*_noisespec.fits)
  * - Write extensions NOISE, NOISEALL and WEIGHTMS (call writeTPSreprExten)
  * - Free allocated GSL vectors
  * - Close output FITS file
  * - Free memory
  * - Finalize the task
  *****************************************************/
 int gennoisespec_main ()
 {
     headas_chat(3, "initialize ...\n");

     message="Running GENNOISESPEC v" + string(SIRENA_VERSION);
     cout<<message<<endl;

     // Setup
     message="=== Setup ===";
     cout<<message<<endl;

     time_t t_start = time(NULL);
     
     int status=EPOK, extver=0;
     
     // Reading all programm parameters by using PIL
     status=getpar_noiseSpec(&par);
     if (status != 0)
     {
         message = "Cannot run getpar_noiseSpec routine to get input parameters";
         EP_EXIT_ERROR(message,status); 
     }
     pulse_length = par.intervalMinSamples;
     
     if ((strcmp(par.EnergyMethod,"I2RFITTED") == 0) && (par.Ifit == 0.0))
     {
         message = "Ifit value must be provided";
         EP_EXIT_ERROR(message,status); 
     }
     
     char str_stat[12];
     char str_stat1[12];
     double cutFreq = 0.;
     int boxLength = 0;
     
     int intervalMinSamples_base2 = pow(2,floor(log2(par.intervalMinSamples)));
     
     // Open input FITS file
     if (fits_open_file(&infileObject, par.inFile,0,&status))
     {
         message = "Cannot open file " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     fits_get_num_hdus(infileObject, &hdunum,&status);
     
     // Check if input FITS file have been simulated with TESSIM or XIFUSIM
     strcpy(extname,"RECORDS");
     fits_movnam_hdu(infileObject, ANY_HDU,extname, extver, &status);
     if (status != 0)
     {
         status = 0;
         strcpy(extname,"TESRECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, extver, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
         else
         {
             tessimOrxifusim = 1;
         }
     }
     else 
     {
         tessimOrxifusim = 0;
     }
     if (tessimOrxifusim == -999)
     {
         message = "Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file";
         EP_EXIT_ERROR(message,status);
     }

     // To calculate 'aducnv'...
     strcpy(extname,"ADCPARAM");
     fits_movnam_hdu(infileObject, ANY_HDU,extname, extver, &status);
     if (status == 0)
     {
         strcpy(keyname,"IMIN");
         if (fits_read_key(infileObject,TDOUBLE,keyname, &Imin,comment,&status))
         {
             status = 0;
         }
         strcpy(keyname,"IMAX");
         if (fits_read_key(infileObject,TDOUBLE,keyname, &Imax,comment,&status))
         {
             status = 0;
         }
     }
     else
     {
         status = 0;
         strcpy(extname,"TESRECORDS");
         fits_movnam_hdu(infileObject,ANY_HDU,extname, 0, &status);
         if (status == 0)
         {
             strcpy(keyname,"IMIN");
             if (fits_read_key(infileObject,TDOUBLE,keyname, &Imin,comment,&status))
             {
                 status = 0;
             }
             strcpy(keyname,"IMAX");
             if (fits_read_key(infileObject,TDOUBLE,keyname, &Imax,comment,&status))
             {
                 status = 0;
             }
         }
         else
         {
             status = 0;
             strcpy(extname,"RECORDS");
             fits_movnam_hdu(infileObject,ANY_HDU,extname, 0, &status);
             strcpy(keyname,"IMIN");
             if (fits_read_key(infileObject,TDOUBLE,keyname, &Imin,comment,&status))
             {
                 status = 0;
             }
             strcpy(keyname,"IMAX");
             if (fits_read_key(infileObject,TDOUBLE,keyname, &Imax,comment,&status))
             {
                 status = 0;
             }
         }
     }

    //...or read ADU_CNV, I_BIAS and ADU_BIAS
    // ADU_CNV(A/ADU)
    strcpy(keyname,"ADU_CNV");
    for (int i=0;i<hdunum;i++)
    {
        fits_movabs_hdu(infileObject, i+1, NULL, &status);
        fits_read_key(infileObject,TDOUBLE,keyname, &adu_cnv,comment,&status);
        if (status == 0)
        {
            adu_cnv_exists = 1;
            break;
        }
        else if ((status != 0) && (i <= hdunum-1))
        {
            status = 0;
        }
    }
    // I_BIAS(A)
    strcpy(keyname,"I_BIAS");
    for (int i=0;i<hdunum;i++)
    {
        fits_movabs_hdu(infileObject, i+1, NULL, &status);
        fits_read_key(infileObject,TDOUBLE,keyname, &i_bias,comment,&status);
        if (status == 0)
        {
            break;
        }
        else if ((status != 0) && (i <= hdunum-1))
        {
            status = 0;
        }
    }
    // ADU_BIAS(ADU)
    strcpy(keyname,"ADU_BIAS");
    for (int i=0;i<hdunum;i++)
    {
        fits_movabs_hdu(infileObject, i+1, NULL, &status);
        fits_read_key(infileObject,TDOUBLE,keyname, &adu_bias,comment,&status);
        if (status == 0)
        {
            break;
        }
        else if ((status != 0) && (i <= hdunum-1))
        {
            status = 0;
        }
    }

     if (strcmp(par.EnergyMethod,"I2RDER") == 0)
     {
        int keyword_exists = 0;

        // V0(V)
        strcpy(keyname,"V0");
        for (int i=0;i<hdunum;i++)
        {
            fits_movabs_hdu(infileObject, i+1, NULL, &status);
            fits_read_key(infileObject,TDOUBLE,keyname, &V0,comment,&status);
            if (status == 0)
            {
                keyword_exists = 1;
                break;
            }
            else if ((status != 0) && (i <= hdunum-1))
            {
                status = 0;
            }
        }
        if (keyword_exists == 0)
        {
            strcpy(extname,"RECORDS");
            fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status);
            if (status != 0)
            {
                status = 0;
                strcpy(extname,"TESRECORDS");
                fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status);
                if (status == 0)
                {
                    IOData objI2RDER;
                    objI2RDER.inObject = infileObject;
                    objI2RDER.nameTable = new char [255];
                    strcpy(objI2RDER.nameTable,"TESPARAM");
                    if (fits_movnam_hdu(infileObject, ANY_HDU,objI2RDER.nameTable, 0, &status))
                    {
                        message = "Cannot move to HDU " + string(objI2RDER.nameTable);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                    }
                    objI2RDER.iniCol = 0;
                    objI2RDER.nameCol = new char [255];
                    objI2RDER.unit = new char [255];
                    objI2RDER.type = TDOUBLE;
                    objI2RDER.iniRow = 1;
                    objI2RDER.endRow = 1;
                    gsl_vector *vectorI2RDER = gsl_vector_alloc(1);
                    strcpy(objI2RDER.nameCol,"V0");
                    if (readFitsSimple (objI2RDER,&vectorI2RDER))
                    {
                        message = "Cannot read V0 keyword in " + string(par.inFile) + " nor " + string(objI2RDER.nameCol) + " column in " + string(extname) + " HDU in " + string(par.inFile);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                    }
                    V0 = gsl_vector_get(vectorI2RDER,0);


                    strcpy(objI2RDER.nameCol,"RL");
                    if (readFitsSimple (objI2RDER,&vectorI2RDER))
                    {
                        message = "Cannot read " + string(objI2RDER.nameCol) + " in " + string(extname) + " HDU in " + string(par.inFile);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                    }
                    RL = gsl_vector_get(vectorI2RDER,0);

                    strcpy(objI2RDER.nameCol,"LIN");
                    if (readFitsSimple (objI2RDER,&vectorI2RDER))
                    {
                        message = "Cannot read " + string(objI2RDER.nameCol) + " in " + string(extname) + " HDU in " + string(par.inFile);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                    }
                    L = gsl_vector_get(vectorI2RDER,0);

                    strcpy(extname,"TESRECORDS");
                    if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
                    {
                        message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
                        EP_EXIT_ERROR(message,status);
                    }

                    gsl_vector_free(vectorI2RDER); vectorI2RDER = 0;
                    delete [] objI2RDER.nameTable; objI2RDER.nameTable = 0;
                    delete [] objI2RDER.nameCol; objI2RDER.nameCol = 0;
                    delete [] objI2RDER.unit; objI2RDER.unit = 0;
                }
            }
        }
        else // keyword_exists(V0) = 1
        {
            keyword_exists = 0;
            // RL(Ohm)
            strcpy(keyname,"RL");
            for (int i=0;i<hdunum;i++)
            {
                fits_movabs_hdu(infileObject, i+1, NULL, &status);
                fits_read_key(infileObject,TDOUBLE,keyname, &RL,comment,&status);
                if (status == 0)
                {
                    keyword_exists = 1;
                    break;
                }
                else if ((status != 0) && (i <= hdunum-1))
                {
                    status = 0;
                }
            }
            if (keyword_exists == 0)
            {
                message = "Cannot read keyword " + string(keyname) + " in input file (to be used in I2RDER)";
                EP_EXIT_ERROR(message,status);
            }
            keyword_exists = 0;

            // L(H)
            strcpy(keyname,"L");
            for (int i=0;i<hdunum;i++)
            {
                fits_movabs_hdu(infileObject, i+1, NULL, &status);
                fits_read_key(infileObject,TDOUBLE,keyname, &L,comment,&status);
                if (status == 0)
                {
                    keyword_exists = 1;
                    break;
                }
                else if ((status != 0) && (i <= hdunum-1))
                {
                    status = 0;
                }
            }
            if (keyword_exists == 0)
            {
                message = "Cannot read keyword " + string(keyname) + " in input file (to be used in I2RDER)";
                EP_EXIT_ERROR(message,status);
            }
        }
     }
     
     if (tessimOrxifusim == 0)
     {
         strcpy(extname,"RECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     else
     {
         strcpy(extname,"TESRECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     
     // Get structure of input FITS file columns
     strcpy(straux,"Time");
     if (fits_get_colnum(infileObject,0,straux,&colnum,&status))
     {
         message = "Cannot get column number for " + string(straux) +" in " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     strcpy(straux,"ADC");
     if (fits_get_colnum(infileObject,0,straux,&colnum,&status))
     {
         message = "Cannot get column number for " + string(straux) +" in " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     
     // Read info to transform to resistance space
     if ((adu_cnv_exists == 0) && ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0)))
     {
         strcpy(extname,"RECORDS");
         fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status);
         if (status != 0)
         {
             status = 0;
             strcpy(extname,"TESRECORDS");
             fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status);
             if (status == 0)
             {                 
                 IOData objI2R;
                 objI2R.inObject = infileObject;
                 objI2R.nameTable = new char [255];
                 strcpy(objI2R.nameTable,"TESPARAM");
                 if (fits_movnam_hdu(infileObject, ANY_HDU,objI2R.nameTable, 0, &status))
                 {
                     message = "Cannot move to HDU " + string(objI2R.nameTable);
                     EP_PRINT_ERROR(message,status); return(EPFAIL);
                 }
                 objI2R.iniCol = 0;
                 objI2R.nameCol = new char [255];
                 objI2R.unit = new char [255];
                 objI2R.type = TDOUBLE;
                 objI2R.iniRow = 1;
                 objI2R.endRow = 1;
                 gsl_vector *vectorI2R = gsl_vector_alloc(1);
                 strcpy(objI2R.nameCol,"I0_START");
                 if (readFitsSimple (objI2R,&vectorI2R))
                 {
                     message = "Cannot read " + string(objI2R.nameCol) + " in " + string(extname) + " HDU in " + string(par.inFile);
                     EP_PRINT_ERROR(message,status); return(EPFAIL);
                 }
                 Ibias = gsl_vector_get(vectorI2R,0);
                 
                 strcpy(extname,"TESRECORDS");
                 if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
                 {
                     message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
                     EP_EXIT_ERROR(message,status);
                 }
                 
                 gsl_vector_free(vectorI2R); vectorI2R = 0;
                 delete [] objI2R.nameTable; objI2R.nameTable = 0;
                 delete [] objI2R.nameCol; objI2R.nameCol = 0;
                 delete [] objI2R.unit; objI2R.unit = 0;
             }
         }
         else
         {
             strcpy(keyname,"I0_START");
             if (fits_read_key(infileObject,TDOUBLE,keyname, &Ibias,NULL,&status))
             {
                 message = "Cannot read neither I0_START in RECORDS HDU in " + string(par.inFile);
                 EP_PRINT_ERROR(message,status); return(EPFAIL);
             }
         }
     }
     
     //Read and check other input keywords
     if (tessimOrxifusim == 0)
     {
         strcpy(extname,"RECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     else
     {
         strcpy(extname,"TESRECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     if (fits_get_num_rows(infileObject,&eventcnt, &status))
     {
         message = "Cannot get number of rows in HDU " + string(extname);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     baseline = gsl_vector_alloc(eventcnt);
     gsl_vector_set_all(baseline,-999.0);
     sigma = gsl_vector_alloc(eventcnt);
     gsl_vector_set_all(sigma,-999.0);
     
     // Read other necessary keywords from ANY HDU
     // Instead of reading TRIGGSZ (xifusim) or RECLEN (tessim), TFORM is used
     snprintf(str_stat,sizeof(str_stat),"%d",colnum);
     message = "TFORM" + string(str_stat);
     strcpy(keyname,message.c_str());
     char readTFORMADC [10];
     fits_read_key(infileObject,TSTRING,keyname,readTFORMADC,comment,&status);
     
     char * pointerTFORM;
     
     pointerTFORM = strstr(readTFORMADC,"(");
     if (pointerTFORM) // There is a parenthesis
     {
         char each_character_after_paren[125];
         char characters_after_paren[125];
     
         pointerTFORM = pointerTFORM + 1; // Pointer to the next character to "(" 
         snprintf(each_character_after_paren,125,"%c",*pointerTFORM);
         snprintf(characters_after_paren,125,"%c",*pointerTFORM);
         while (*pointerTFORM != ')')
         {
             pointerTFORM = pointerTFORM + 1;
             snprintf(each_character_after_paren,125,"%c",*pointerTFORM);
             strcat(characters_after_paren,each_character_after_paren); 
         }
         eventsz = atoi(characters_after_paren);
     }
     else    // There is not a parenthesis
     {       
         string readTFORMADCstring;
         readTFORMADCstring = readTFORMADC;
         
         int index = 0;
         string characters;
         while (isdigit(readTFORMADCstring[index]))
         {
             if (index == 0) characters = readTFORMADCstring[index];
             else            characters = characters + readTFORMADCstring[index];
             
             index++;
         }
         eventsz = stoi(characters);
         
         readTFORMADCstring.clear();
         characters.clear();
     }
     
     // Calculate the sampling rate 
     // By using keywords in input FITS file (from DELTAT or TCLOCK+DEC_FAC or NUMROW+P_ROW)
     int deltat_exists = 0;
     int dec_fac_exists = 0;
     int tclock_exists = 0;
     int numrow_exists = 0;
     int p_row_exists = 0;
     for (int i=0;i<hdunum;i++)
     {
         fits_movabs_hdu(infileObject, i+1, NULL, &status); 
         fits_read_key(infileObject,TDOUBLE,"DELTAT", &deltat,comment,&status);
         if (status == 0)
         {
             deltat_exists = 1;
             samprate = 1/deltat;
             break;
         }
         else if ((status != 0) && (i <= hdunum-1))
         {
             status = 0;
         }
     }
     double tclock;
     if (deltat_exists == 0)
     {
         double dec_fac;
         for (int i=0;i<hdunum;i++)
         {
             fits_movabs_hdu(infileObject, i+1, NULL, &status); 
             fits_read_key(infileObject,TDOUBLE,"DEC_FAC", &dec_fac,comment,&status);
             if (status == 0)
             {
                 dec_fac_exists = 1;
                 break;
             }
             else if ((status != 0) && (i <= hdunum-1))
             {
                 status = 0;
             }
         }
         for (int i=0;i<hdunum;i++)
         {
             fits_movabs_hdu(infileObject, i+1, NULL, &status); 
             fits_read_key(infileObject,TDOUBLE,"TCLOCK", &tclock,comment,&status);
             if (status == 0)
             {
                 tclock_exists = 1;
                 break;
             }
             else if ((status != 0) && (i <= hdunum-1))
             {
                 status = 0;
             }
         }
         if ((dec_fac_exists == 1) && (tclock_exists == 1)) 
         {
             samprate = 1/(tclock*dec_fac);
         }
     }
     if ((deltat_exists == 0) && ((tclock_exists == 0) || (dec_fac_exists == 0)))
     {
         int numrow;
         for (int i=0;i<hdunum;i++)
         {
             fits_movabs_hdu(infileObject, i+1, NULL, &status); 
             fits_read_key(infileObject,TINT,"NUMROW", &numrow,comment,&status);
             if (status == 0)
             {
                 numrow_exists = 1;
                 break;
             }
             else if ((status != 0) && (i <= hdunum-1))
             {
                 status = 0;
             }
         }
         int p_row;
         for (int i=0;i<hdunum;i++)
         {
             fits_movabs_hdu(infileObject, i+1, NULL, &status); 
             fits_read_key(infileObject,TINT,"P_ROW", &p_row,comment,&status);
             if (status == 0)
             {
                 p_row_exists = 1;
                 break;
             }
             else if ((status != 0) && (i <= hdunum-1))
             {
                 status = 0;
             }
         }
         if ((tclock_exists == 1) && (numrow_exists == 1) && (p_row_exists == 1)) 
         {
             samprate =1/(tclock*numrow*p_row);
         }
     }
     
     // If necessary...
     //...read the sampling rate from input FITS file (from the HISTORY in the Primary HDU) 
     //if ((tessimOrxifusim == 1) && (samprate == -999.0))
     if (samprate == -999.0)
     {
         fits_movabs_hdu(infileObject, 1, NULL, &status); // Move to "Primary" HDU
         int numberkeywords;
         char *headerPrimary;
         fits_hdr2str(infileObject, 0, NULL, 0,&headerPrimary, &numberkeywords, &status);   // Reading thee whole "Primary" HDU and store it in 'headerPrimary'
         char * sample_rate_pointer;
         sample_rate_pointer = strstr (headerPrimary,"sample_rate=");    // Pointer to where the text "sample_rate=" is
         if(sample_rate_pointer){
             sample_rate_pointer = sample_rate_pointer + 12; // Pointer to the next character to "sample_rate=" (which has 12 characters)   
             char each_character_after_srate[125];		
             snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
             char characters_after_srate[125];
             snprintf(characters_after_srate,125,"%c",*sample_rate_pointer);
             while (*sample_rate_pointer != ' ')
             {
                 sample_rate_pointer = sample_rate_pointer + 1;
                 snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                 strcat(characters_after_srate,each_character_after_srate); 
             }
             samprate = atof(characters_after_srate);
         }
     }
     
     // If it is not possible neither to calculate nor get the sampling rate from the input FITS file,...
     //...provide an error message to include DELTAT (inverse of sampling rate) in the input FITS file before running again
     if (samprate == -999.0)
     {
         message =  "Cannot read or get the sampling rate from the input file. Please, include the DELTAT keyword (inverse of sampling rate) in the input FITS file before running GENNOISESPEC again";
         
         EP_EXIT_ERROR(message,EPFAIL);
     }
          
     if (tessimOrxifusim == 0)
     {
         strcpy(extname,"RECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     else
     {
         strcpy(extname,"TESRECORDS");
         if (fits_movnam_hdu(infileObject, ANY_HDU,extname, 0, &status))
         {
             message = "Cannot move to HDU " + string(extname) + " in " + string(par.inFile);
             EP_EXIT_ERROR(message,status);
         }
     }
     
     if (eventsz <= 0)
     {
         message = "Legal values for TRIGGSZ (RECORDS) are integer numbers greater than 0";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     ivcal=1.0;
     if ((((Imin != -999.0) && (Imax != -999.0)) && ((Imin != 0) && (Imax != 0))) && (adu_cnv != -999.0) && ((Imax-Imin)/65534 != -1.0*adu_cnv))
     {
         message = "ADU_CNV and (Imax-Imin)/65534 do not match";
         EP_PRINT_ERROR(message,-EPFAIL);
     }
     
     if ((((Imin == -999.0) || (Imax == -999.0)) || ((Imin == 0) || (Imax == 0))) && (adu_cnv == -999.0))
     {
         aducnv = 1.0;
         message = "ADU_CNV not found, or Imin or Imax not found, or both equal to 0 => Conversion factor ('aducnv' to convert adu into A) is fix to 1";
         EP_PRINT_ERROR(message,-999);	// Only a warning
     }
     else if (adu_cnv != -999.0)
     {
         aducnv = -1*adu_cnv;
     }
     else if (((Imin != -999.0) && (Imax != -999.0)) && ((Imin != 0) && (Imax != 0)))
     {
         aducnv = (Imax-Imin)/65534;    // Quantification levels = 65534  // If this calculus changes => Change it also in TASKSSIRENA
     }
     
     asquid = 1.0;
     plspolar = 1.0;
     
     // Declare variables
     intervalMinBins = par.intervalMinSamples;
     if (intervalMinBins > eventsz)
     {
         message = "Illegal value in INTERVALMINSAMPLES parameter. Legal values reals greater than 0 and fewer than record length";
         EP_EXIT_ERROR(message,EPFAIL); 
     }
     
     // Maximum number of noise intervals in all the records eventcnt*floor(eventsz/par.intervalMinSamples)
     noiseIntervals = gsl_matrix_alloc(eventcnt*floor(eventsz/par.intervalMinSamples),intervalMinBins);
     
     EventSamplesFFTMean = gsl_vector_alloc(intervalMinBins);
     gsl_vector_set_zero(EventSamplesFFTMean);
     gsl_vector *mean;
     mean = gsl_vector_alloc(intervalMinBins);
     gsl_vector_set_zero(mean);
     sigmacsdgsl = gsl_vector_alloc(intervalMinBins);
     gsl_vector_set_zero(sigmacsdgsl);
     
     // Create structure to run Iteration
     extern int inDataIterator(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );
     
     long rows_per_loop = 0;		// 0: Use default: Optimum number of rows
     long offset = 0;		// 0: Process all the rows
     iteratorCol cols [2];		// Structure of Iteration
     int n_cols = 2; 		// Number of columns:  Time + ADC
     
     // Read Time column
     strcpy(straux,"TIME");
     status = fits_iter_set_by_name(&cols[0], infileObject, straux, TDOUBLE, InputCol);
     if (status)
     {
         message = "Cannot iterate in column " + string(straux) +" in " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     
     // Read ADC column
     strcpy(straux,"ADC");
     status = fits_iter_set_by_name(&cols[1], infileObject, straux, TDOUBLE, InputCol);
     if (status)
     {
         message = "Cannot iterate in column " + string(straux) +" in " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     
     // Check Boxlength
     cutFreq = 2 * (1/(2*pi*par.scaleFactor));
     boxLength = (int) ((1/cutFreq) * samprate);
     if (boxLength <= 1)
     {
         message = "lpf_boxcar: scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
         EP_PRINT_ERROR(message,-999);
     }

     // Called iteration function
     if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))
     {
         message = "Cannot iterate data using InDataTterator";
         EP_EXIT_ERROR(message,status);
     }
     
     // Close input FITS file
     if (fits_close_file(infileObject,&status))
     {
         message = "Cannot close file " + string(par.inFile);
         EP_EXIT_ERROR(message,status);
     }
     infileObject = 0;
     
     // Applying medianKappaClipping in order to remove the noise intervals with a high sigma
     gsl_vector *interval = gsl_vector_alloc(noiseIntervals->size2);
     double bsln, sgm;
     gsl_vector *sigmaInterval = gsl_vector_alloc(NumMeanSamples);
     int cnt = NumMeanSamples;   // After removing the noise intervals with a too high sigma, it is going to be the number of noise intervals with a proper sigma
     
     double meanThreshold;
     double sgmThreshold;
     for (int i=0;i<NumMeanSamples;i++)
     {
         gsl_matrix_get_row(interval,noiseIntervals,i);
         findMeanSigma(interval, &bsln, &sgm);
         gsl_vector_set(sigmaInterval,i, sgm);
     }
     
     if ((NumMeanSamples > 1) && (par.rmNoiseIntervals == 1))
     {
        if (medianKappaClipping_noiseSigma (sigmaInterval, kappaMKC, stopCriteriaMKC, &meanThreshold, &sgmThreshold))
        {
            message = "Cannot run medianKappaClipping_noiseSigma looking for mean and sigma of the noise sigmas";
            EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
        }
     }

     gsl_vector *intervalsgmOK = gsl_vector_alloc(NumMeanSamples);
     gsl_vector_set_all(intervalsgmOK,1);
     gsl_vector *vector_aux;
     gsl_vector_complex *vector_aux1;
     vector_aux = gsl_vector_alloc(intervalMinBins);
     vector_aux1 = gsl_vector_complex_alloc(intervalMinBins);
     double nSgms_sigmaInterval = 1;
     int NumMeanSamples_afterRm = 0;
     printf("\n"); // New line after ending "Reading the noise file...." (progressbar READING NOISE FILE)
     for (int i=0;i<NumMeanSamples;i++)
     {
         if ((par.rmNoiseIntervals == 1) && (NumMeanSamples > 1) && ((gsl_vector_get(sigmaInterval,i) > meanThreshold+nSgms_sigmaInterval*sgmThreshold)||(gsl_vector_get(sigmaInterval,i) < meanThreshold-nSgms_sigmaInterval*sgmThreshold)))
         {
             // Interval not to be taken account
             cnt --;
             gsl_vector_set(intervalsgmOK,i,0);
         }
         else
         {
             if (NumMeanSamples_afterRm < par.nintervals)
             {
                NumMeanSamples_afterRm++;
                
                gsl_matrix_get_row(interval,noiseIntervals,i);
                
                // Apply a Hanning window to reduce spectral leakage
                /*if (hannWindow(&interval))
                {
                    message = "Cannot run hannWindow routine";
                    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
                }*/
    
                // FFT calculus (EventSamplesFFT)
                if(FFT(interval,vector_aux1))
                {
                    message = "Cannot run FFT routine for vector1";
                    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
                }
                gsl_vector_complex_absIFCA(vector_aux,vector_aux1);
                
                // Add to mean FFT samples
                gsl_vector_mul(vector_aux,vector_aux);
                gsl_vector_add(EventSamplesFFTMean,vector_aux);
             }
         }
         // Progress bar NOISE SPECTRUM
         float progress = (float)i / (NumMeanSamples-1);
         int bar_width = 50;
         int pos = bar_width * progress;
         printf("Generating the noise spectrum |");
         for (int j = 0; j < bar_width; j++) {
             if (j < pos)
                 printf("=");
             else if (j == pos)
                 printf(">");
             else
                 printf(" ");
         }
         printf("| %.2f%%\r", progress * 100);
         fflush(stdout);
         if (i == NumMeanSamples-1)
         {
             printf("\n"); // New line after ending "Generating the noise spectrum...."
         }
     }
     gsl_vector_free(interval); interval = 0;
     gsl_vector_free(sigmaInterval); sigmaInterval = 0;
     gsl_vector_free(vector_aux); vector_aux = 0;
     gsl_vector_complex_free(vector_aux1); vector_aux1 = 0;
     
     // Current noise spectral density
     // sqrt(sum(FFT^2)/NumMeanSamplesCSD) => sqrt(A^2) = A and sqrt(1/NumMeanSamplesCSD)=1/sqrt(Hz)
     gsl_vector_scale(EventSamplesFFTMean,(1.0/(double)NumMeanSamples_afterRm));
     for (int i=0;i<(int)(EventSamplesFFTMean->size);i++)
     {
         if (gsl_vector_get(EventSamplesFFTMean,i)<0)
         {
             message = "EventSamplesFFTMean < 0 for i=" + boost::lexical_cast<std::string>(i);
             EP_EXIT_ERROR(message,EPFAIL); 
         }
     }
     gsl_vector_sqrtIFCA(EventSamplesFFTMean,EventSamplesFFTMean);
     
     // Extra normalization (further than the FFT normalization factor,1/n) in order 
     // to get the apropriate noise level provided by Peille (54 pA/rHz)
     gsl_vector_scale(EventSamplesFFTMean,sqrt(2*intervalMinBins/samprate));
     
     // Load in noiseIntervals only those intervals with a proper sigma and NumMeanSamples = cnt
     // (in order not to change excessively the code when weightMS)
     gsl_matrix *matrixaux = gsl_matrix_alloc(noiseIntervals->size1,noiseIntervals->size2);
     gsl_vector *vectoraux = gsl_vector_alloc(noiseIntervals->size2);
     gsl_matrix_memcpy(matrixaux,noiseIntervals);
     gsl_matrix_free(noiseIntervals); 
     noiseIntervals = gsl_matrix_alloc(cnt,matrixaux->size2);
     int ind = 0;
     for (int i=0;i<(int)(intervalsgmOK->size);i++)
     {
         if (gsl_vector_get(intervalsgmOK,i) == 1)
         {
             gsl_matrix_get_row(vectoraux,matrixaux,i);
             gsl_matrix_set_row(noiseIntervals,ind,vectoraux);
             ind++;
         }
     }
     gsl_matrix_free(matrixaux); matrixaux = 0;
     gsl_vector_free(vectoraux); vectoraux = 0;
     NumMeanSamples = cnt;
     
     // Generate WEIGHT representation
     if ((par.weightMS == 1) && (NumMeanSamples != 0))
     {
         weightpoints = gsl_vector_alloc(floor(log2(intervalMinSamples_base2)));
         for (int i=0;i<(int)(weightpoints->size);i++) 		gsl_vector_set(weightpoints,i,pow(2,floor(log2(intervalMinSamples_base2))-i));
         weightMatrixes = gsl_matrix_alloc(weightpoints->size,intervalMinSamples_base2*intervalMinSamples_base2);
         gsl_matrix_set_all(weightMatrixes,-999.0);
         gsl_matrix_view tempm;
         gsl_matrix *noiseIntervals_weightPoints;
         gsl_matrix *weightMatrix;

         if (NumMeanSamples >= par.nintervals)
         {
             cout<<"Calculating the weight noise matrixes [0"<<"/"<<weightpoints->size<<"][";
             for (int i=0;i<(int)(weightpoints->size);i++) cout<<".";
             cout<<"]"<<flush;
             for (int i=0;i<(int)(weightpoints->size);i++)
             {	
                 weightMatrix = gsl_matrix_alloc(gsl_vector_get(weightpoints,i),gsl_vector_get(weightpoints,i));
                 noiseIntervals_weightPoints = gsl_matrix_alloc(cnt,gsl_vector_get(weightpoints,i));
                 
                 tempm = gsl_matrix_submatrix(noiseIntervals,0,0,cnt,gsl_vector_get(weightpoints,i));
                 gsl_matrix_memcpy(noiseIntervals_weightPoints,&tempm.matrix);
                 
                 weightMatrixNoise(noiseIntervals_weightPoints, &weightMatrix);
                 for (int j=0;j<gsl_vector_get(weightpoints,i);j++)
                 {
                     for (int k=0;k<gsl_vector_get(weightpoints,i);k++)
                     {
                         gsl_matrix_set(weightMatrixes,i,j*gsl_vector_get(weightpoints,i)+k,gsl_matrix_get(weightMatrix,j,k));
                     }
                 }
                 gsl_matrix_free(noiseIntervals_weightPoints);
                 gsl_matrix_free(weightMatrix);
                     
                 // Progress bar WEIGHT NOISE MATRIXES
                 cout<<"\rCalculating the weight noise matrixes [";
                 cout<<i+1;
                 cout<<"/";
                 cout<<weightpoints->size;
                 cout<<"][";
                 for (int j = 0; j < (int)(weightpoints->size); j++)
                 {
                     if (j <= i) {
                         cout << "#";  // Write '#' for each completed matrix
                     } else {
                         cout << ".";  // Write '.' for pending matrixes
                     }
                 }
                 cout << "]" << flush;
                 if (i == (int)(weightpoints->size)-1)
                 {
                     printf("\n"); // New line after ending "Calculating the weight noise matrixes...."
                 }
             } // foreach matrix size
         }
         else
         {
             cout<<"Calculating the weight noise matrixes [0"<<"/"<<weightpoints->size<<"][";
             for (int i=0;i<(int)(weightpoints->size);i++) cout<<".";
             cout<<"]"<<flush;
             for (int i=0;i<(int)(weightpoints->size);i++)
             {	
                 weightMatrix = gsl_matrix_alloc(gsl_vector_get(weightpoints,i),gsl_vector_get(weightpoints,i));
                 noiseIntervals_weightPoints = gsl_matrix_alloc(NumMeanSamples,gsl_vector_get(weightpoints,i));
                 
                 tempm = gsl_matrix_submatrix(noiseIntervals,0,0,NumMeanSamples,gsl_vector_get(weightpoints,i));
                 gsl_matrix_memcpy(noiseIntervals_weightPoints,&tempm.matrix);

                 weightMatrixNoise(noiseIntervals_weightPoints, &weightMatrix);
                 for (int j=0;j<gsl_vector_get(weightpoints,i);j++)
                 {
                     for (int k=0;k<gsl_vector_get(weightpoints,i);k++)
                     {
                         gsl_matrix_set(weightMatrixes,i,j*gsl_vector_get(weightpoints,i)+k,gsl_matrix_get(weightMatrix,j,k));
                     }
                 }
                 gsl_matrix_free(noiseIntervals_weightPoints);
                 gsl_matrix_free(weightMatrix);
                     
                 // Progress bar WEIGHT NOISE MATRIXES
                 cout<<"\rCalculating the weight noise matrixes [";
                 cout<<i+1;
                 cout<<"/";
                 cout<<weightpoints->size;
                 cout<<"][";
                 for (int j = 0; j < (int)(weightpoints->size); j++)
                 {
                     if (j <= i) {
                         cout << "#";  // Write '#' for each completed matrix
                     } else {
                         cout << ".";  // Write '.' for pending matrixes
                     }
                 }
                 cout << "]" << flush;
                 if (i == (int)(weightpoints->size)-1)
                 {
                     printf("\n"); // New line after ending "Calculating the weight noise matrixes...."
                 }
             } // foreach matrix size
         }
     }
     gsl_vector_free(intervalsgmOK); intervalsgmOK = 0;

     if (cnt == 0)
     {
         message = "Pulse-free intervals not found";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     else if	(cnt < par.nintervals)
     {
         snprintf(str_stat,sizeof(str_stat),"%d",cnt);
         snprintf(str_stat1,sizeof(str_stat1),"%d",intervalMinBins);
         if (par.weightMS == 1)
            message = "Not enough pulse-free intervals for calculus. CSD and W" + string(str_stat1) + " matrix calculated with " + string(str_stat);
         else
            message = "Not enough pulse-free intervals for calculus. CSD calculated with " + string(str_stat);
         cout<<message<<endl;
     }
     else if	(cnt >= par.nintervals)
     {
         snprintf(str_stat,sizeof(str_stat),"%d",par.nintervals);
         if (par.weightMS == 1)
            message = "CSD and all Wx matrices calculated with " + string(str_stat) + " intervals";
         else
            message = "CSD calculated with " + string(str_stat) + " intervals";
         cout<<message<<endl;
     }
     
     // Create output FITS File: GENNOISESPEC representation file (*_noisespec.fits)
     if(createTPSreprFile())
     {
         message = "Cannot create file " +  string(par.outFile);
         EP_EXIT_ERROR(message,EPFAIL);
     }	
     
     if (fits_open_file(&gnoiseObject,par.outFile,1,&status))
     {
         message = "Cannot open file " +  string(par.outFile);
         EP_EXIT_ERROR(message,status);
     }

     //HDpar_stamp(gnoiseObject, 0, &status);  // Write the whole list of input parameters in HISTORY______ Does not work in C++

     // Write extensions NOISE and NOISEALL (call writeTPSreprExten)
     if(writeTPSreprExten ())
     {
         message = "Cannot write extensions in " +  string(par.outFile);
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     // Free allocated GSL vectors
     gsl_vector_free(EventSamplesFFTMean); EventSamplesFFTMean = 0;
     gsl_vector_free(mean); mean = 0;
     gsl_vector_free(sigmacsdgsl); sigmacsdgsl = 0;
     gsl_vector_free(startIntervalgsl); startIntervalgsl = 0;
     
     gsl_vector_free(baseline); baseline = 0;
     gsl_vector_free(sigma); sigma = 0;     
     
     gsl_matrix_free(noiseIntervals); noiseIntervals = 0;
     if (weightMS == 1)
     {
         gsl_vector_free(weightpoints); weightpoints = 0;
         gsl_matrix_free(weightMatrixes); weightMatrixes = 0;
     }
     
     // Close output FITS file
     if (fits_close_file(gnoiseObject,&status))
     {
         message = "Cannot close file " + string(par.outFile);
         EP_EXIT_ERROR(message,status);
     }	
     gnoiseObject = 0;
     
     // Free memory
     delete [] obj.nameTable; obj.nameTable = 0;
     delete [] obj.nameCol; obj.nameCol = 0;
     delete [] obj.unit; obj.unit = 0;
     
     // Finalize the task
     time_t t_end = time(NULL);
     sprintf(straux,"%f",(double) (t_end - t_start));
     message = "Elapsed time:" + string(straux);
     cout<<message<<endl;
     
     message = "Gennoisespec Module OK";
     cout<<message<<endl;
     
     message.clear();
     
     return EPOK;
 }
 /*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 2 ************************************************************
  * inDataIterator: This function takes the optimum number of rows to read the input FITS file
  *                 and works iteratively
  *
  * - Declare variables
  * - Allocate input GSL vectors
  * - Read iterator
  * - Processing each record
  *     - Information has been read by blocks (with nrows per block)
  *     - Just in case the last record has been filled out with 0's => Last record discarded
  * 	- Convert to the resistance space if necessary
  *     - To avoid taking into account the pulse tails at the beginning of a record as part of a pulse-free interval
  * 	- Low-pass filtering
  *   	- Differentiate 
  *   	- Finding the pulses: Pulses tstarts are found (call findPulsesNoise)
  *       - Finding the pulse-free intervals in each record
  *  	    - If there are pulses => Call findInterval
  *	    - No pulses => The whole event is going to be used (DIVIDING into intervals of intervalMinBins size) => Call findIntervalN
  *       - Calculating the mean and sigma of the intervals without pulses together => BSLN0 & NOISESTD
  *       - Preparing the CSD calculus (not filtered data)
  * - Free allocated GSL vectors
  ****************************************************************************/
 int inDataIterator(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
 {
     int status = EPOK;
     
     // Declare variables
     // To read from the input FITS file
     double *time_block, *timein;	// Vector of TIME column
     double *iout_block, *ioutin;	// Vector of ADC column
     gsl_vector *timegsl_block;
     gsl_matrix *ioutgsl_block;
     
     // To differentiate
     gsl_vector *derSGN = gsl_vector_alloc(eventsz);
     
     gsl_vector *tstartgsl = gsl_vector_alloc(eventsz);
     gsl_vector *qualitygsl = gsl_vector_alloc(eventsz);
     int nPulses;
     double threshold;

     gsl_vector *tstartDERgsl = gsl_vector_alloc(eventsz);
     gsl_vector *tmaxDERgsl = gsl_vector_alloc(eventsz);
     gsl_vector *maxDERgsl = gsl_vector_alloc(eventsz);
     gsl_vector *tendDERgsl = gsl_vector_alloc(eventsz);
     
     int pulseFound = 0;	// 0->The function findTstart has not found any pulse
     // 1->The function findTstart has found at least one pulse
     
     int tail_duration;
     
     double baselineIntervalFreeOfPulses;
     double sigmaIntervalFreeOfPulses;
     
     // Auxiliary variables
     gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)
     
     // To calculate the FFT
     gsl_vector *EventSamples = gsl_vector_alloc(intervalMinBins);
     
     // Allocate input GSL vectors
     timegsl_block = gsl_vector_alloc(nrows);
     timegsl = gsl_vector_alloc(nrows); 			// GSL of input TIME column
     ioutgsl_block = gsl_matrix_alloc(nrows,eventsz);
     ioutgsl = gsl_vector_alloc(eventsz); 			// GSL of input ADC column
     gsl_vector *ioutgsl_aux = gsl_vector_alloc(eventsz);	//In case of working with ioutgsl without filtering
     gsl_vector *ioutgslNOTFIL = gsl_vector_alloc(eventsz);
     gsl_vector *ioutgslFIL = gsl_vector_alloc(eventsz);
     
     // Read iterator
     // NOTE: fits_iter_get_array because in this fits function the 1st element
     //  of the output array is the null pixel value!
     timein = (double *) fits_iter_get_array(&cols[0]);
     time_block = &timein[1];
     if (toGslVector(((void **)&time_block), &timegsl_block, nrows, 0, TDOUBLE))
     {
         message = "Cannot run routine toGslVector for TIME vector";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     ioutin = (double *) fits_iter_get_array(&cols[1]);
     iout_block = &ioutin[1];
     if(toGslMatrix(((void **)&iout_block), &ioutgsl_block, eventsz, nrows, TDOUBLE,	0))
     {
         message = "Cannot run routine toGslMatrix for IOUT ";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     
     // Processing each record
     for (int i=0; i< nrows; i++)
     {      
         // Information has been read by blocks (with nrows per block)
         // Now, information is going to be used by rows
         gsl_vector_memcpy(timegsl,timegsl_block);
         gsl_matrix_get_row(ioutgsl,ioutgsl_block,i);
         gsl_vector_scale(ioutgsl,ivcal);		//IVCAL to change arbitrary units of voltage to non-arbitrary units of current (Amps)

         // Just in case the last record has been filled out with 0's => Last record discarded
         if ((gsl_vector_get(ioutgsl,ioutgsl->size-1) == 0) && (gsl_vector_get(ioutgsl,ioutgsl->size-2) == 0))		break;
         
         // Convert to the resistance space if necessary
         if (strcmp(par.EnergyMethod,"OPTFILT") != 0)
         {
             if ((tessimOrxifusim = 1) && (hdunum == 2))
             {
                 real_data = 1;
             }
             if (convertI2R(par.EnergyMethod,Ibias,Imin,Imax,adu_cnv,adu_bias,i_bias,par.Ifit,V0,RL,L,&ioutgsl,real_data))
             {
                 message = "Cannot run routine convertI2R";
                 EP_EXIT_ERROR(message,EPFAIL);
             }
         }
         else
         {

            gsl_vector_scale(ioutgsl,aducnv);
         }
         
         // Assigning positive polarity (by using ASQUID and PLSPOLAR)
         gsl_vector_memcpy(ioutgsl_aux,ioutgsl);
         if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))	gsl_vector_scale(ioutgsl_aux,-1.0);
         gsl_vector_memcpy(ioutgslNOTFIL,ioutgsl_aux);
         
         // To avoid taking into account the pulse tails at the beginning of a record as part of a pulse-free interval
         tail_duration = 0;
         
         // Low-pass filtering
         status = lpf_boxcar(&ioutgsl_aux,ioutgsl_aux->size,par.scaleFactor,samprate);
         if (status == EPFAIL)
         {
             message = "Cannot run routine lpf_boxcar for low pass filtering";
             EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
         }
         if (status == 3)
         {
             status = EPOK;
         }
         if (status == 4)
         {
             message = "lpf_boxcar: scaleFactor too high => Cut-off frequency too low";
             EP_PRINT_ERROR(message,EPFAIL); 
         }
         gsl_vector_memcpy(ioutgslFIL,ioutgsl_aux);
         
         // Differentiate
         if (differentiate (&ioutgsl_aux,ioutgsl_aux->size))
         {
             message = "Cannot run routine differentiate for differentiating after filtering";
             EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
         }
         
         //Finding the pulses: Pulses tstarts are found
         if (findPulsesNoise (ioutgsl_aux, &tstartgsl, &qualitygsl,
             &nPulses, &threshold,
             par.scaleFactor, samprate,
             par.samplesUp, par.nSgms,
             stopCriteriaMKC,
             kappaMKC))
         {
             message = "Cannot run routine findPulsesNoise";
             EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
         }
         
         if (nPulses != 0)	pulseFound = 1;
         
         nPulses = 0;
         
         if ((pulseFound == 1) || (tail_duration != 0))
         {
             // Finding the pulse-free intervals in each record 
             if (findInterval(tail_duration, ioutgsl_aux, tstartgsl, nPulses, pulse_length, par.nplPF, intervalMinBins, &nIntervals, &startIntervalgsl))
             {
                 message = "Cannot run routine findInterval";
                 EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
             }
         }
         else if (pulseFound == 0)
         {
             // The whole event is going to be used: It must be divided into intervals of intervalMinBins size
             if (findIntervalN(ioutgsl_aux,intervalMinBins,&nIntervals,&startIntervalgsl))
             {
                 message = "Cannot run routine findIntervalN";
                 EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
             }
         }
         
         if (par.scaleFactor != 0)
         {
             gsl_vector_memcpy(ioutgsl,ioutgslFIL);
         }
         
         // Calculating the mean and sigma of the intervals without pulses together => BSLN0 & NOISESTD
         gsl_vector *intervalsWithoutPulsesTogether = gsl_vector_alloc(nIntervals*intervalMinBins);
         for (int k=0; k<nIntervals; k++)
         {
             for (int j=0; j<intervalMinBins; j++)
             {
                 gsl_vector_set(intervalsWithoutPulsesTogether,intervalMinBins*k+j,gsl_vector_get(ioutgsl,j+gsl_vector_get(startIntervalgsl,k)));
             }
         }
         
         if (nIntervals != 0)
         {
             findMeanSigma (intervalsWithoutPulsesTogether, &baselineIntervalFreeOfPulses, &sigmaIntervalFreeOfPulses);
             gsl_vector_set(baseline,indexBaseline,baselineIntervalFreeOfPulses);
             gsl_vector_set(sigma,indexBaseline,sigmaIntervalFreeOfPulses);
             indexBaseline++;
         }

         gsl_vector_free(intervalsWithoutPulsesTogether); intervalsWithoutPulsesTogether = 0;

         // Preparing the CSD calculus (not filtered data)
         for (int k=0; k<nIntervals;k++)
         {
             temp = gsl_vector_subvector(ioutgsl,gsl_vector_get(startIntervalgsl,k), intervalMinBins);
             gsl_vector_memcpy(EventSamples,&temp.vector);
          
             gsl_matrix_set_row(noiseIntervals,NumMeanSamples,EventSamples);
                 
             NumMeanSamples = NumMeanSamples + 1;
         }

         if (ntotalrows != totalrows)
         {
             // Progress bar READING NOISE FILE
             float progress = (float)ntotalrows / (totalrows-1);
             int bar_width = 50;
             int pos = bar_width * progress;
             printf("Reading the noise file        |");
             for (int j = 0; j < bar_width; j++) {
                 if (j < pos)
                     printf("=");
                 else if (j == pos)
                     printf(">");
                 else
                     printf(" ");
             }
             printf("| %.2f%%\r", progress * 100);
             fflush(stdout);
         }
         
         ntotalrows++;
     }

     // Free allocated GSL vectors
     gsl_vector_free(timegsl); timegsl = 0;
     gsl_vector_free(ioutgsl); ioutgsl = 0;
     gsl_vector_free(ioutgsl_aux); ioutgsl_aux = 0;
     gsl_vector_free(ioutgslFIL); ioutgslFIL = 0;
     gsl_vector_free(timegsl_block); timegsl_block = 0;
     gsl_matrix_free(ioutgsl_block); ioutgsl_block = 0;
     gsl_vector_free(derSGN); derSGN = 0;
     gsl_vector_free(tstartgsl); tstartgsl = 0;
     gsl_vector_free(tstartDERgsl); tstartDERgsl = 0;
     gsl_vector_free(tmaxDERgsl); tmaxDERgsl = 0;
     gsl_vector_free(maxDERgsl); maxDERgsl = 0;
     gsl_vector_free(tendDERgsl); tendDERgsl = 0;
     
     gsl_vector_free(EventSamples); EventSamples = 0;
     
     return (EPOK);
 }
 /*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 3 ************************************************************
  * findInterval: This function finds the pulse-free intervals when the input vector has pulses.
  *               The pulse-free intervals must have a minimum length (intervalMinBins).
  *               The interval with pulse is Tstart,Tend+nPF*pulse_length (being Tend=pulse_length).
  *
  * Parameters:
  * - tail_duration: Length of the tail of a previous pulse
  * - invector: Input vector WITH pulses
  * - startpulse: Vector with the Tstart of all the pulses of the input vector (samples)
  * - npin: Number of pulses in the input vector
  * - pulse_length: Pulse length (samples)
  * - nPF: Number of pulse lengths after ending the pulse to start the pulse-free interval
  * - interval: Minimum length of the interval (samples)
  * - ni: Number of pulse-free intervals in the input vector
  * - startinterval: Vector with the starting time of each pulse-free interval (samples).
  *
  * - Declare variables
  * - Processing if the input vector has more pulses
  * 	- It looks for pulse-free intervals between pulses
  * - Processing if there are no more pulses in the input vector
  * 	-It looks for pulse-free intervals at the end of the event and the search for more pulse-free intervals is finished
  *****************************************************************************/
 int findInterval(int tail_duration, gsl_vector *invector, gsl_vector *startpulse, int npin, int local_pulse_length, int nPF, int interval,
                  int *ni, gsl_vector **startinterval)
 {
     // Declare variables
     long start;		// Auxiliary variable, possible starting time of a pulse-free interval
     int niinc;		// Increase the number of pulse-free intervals
     *ni = 0;
     
     // startinterval is allocated with the maximum possible pulse-free intervals
     // Not necessary handling division by zero because of interval being always greater than zero
     *startinterval = gsl_vector_alloc(invector->size/interval+1);
     
     if ((tail_duration != 0) && (tail_duration >= gsl_vector_get(startpulse,0)))	start = gsl_vector_get(startpulse,0);
     else										start = tail_duration;
     
     for (int i=0; i<npin; i++)
     {
         //If the input vector has pulses, it looks for pulse-free intervals between pulses
         if (gsl_vector_get(startpulse,i)-start >= 0)
         {
             niinc = (gsl_vector_get(startpulse,i)-start)/interval;
             
             for (int j=*ni; j<*ni+niinc; j++)
             {
                 gsl_vector_set(*startinterval,j,start+(j-*ni)*interval);
             }
             start = gsl_vector_get(startpulse,i)+local_pulse_length+nPF*local_pulse_length;
             *ni = *ni+niinc;
         }
         else
             //The tend+nPF*pulse_length of a pulse is aliased to the start of the next pulse
         {
             if (gsl_vector_get(startpulse,i)-gsl_vector_get(startpulse,i-1)<0)
             {
                 start = gsl_vector_get(startpulse,i-1)+local_pulse_length+nPF*local_pulse_length;
             }
             else
             {
                 start = gsl_vector_get(startpulse,i)+local_pulse_length+nPF*local_pulse_length;
             }
         }
     }
     
     if ((int) (invector->size-start) >=interval)
     {
         //If there are no more pulses in the event, it looks for pulse-free intervals at the end of the event
         //and the search for more pulse-free intervals is finished
         
         niinc = (invector->size-start)/interval;
         
         for (int j=*ni; j<*ni+niinc; j++)
         {
             gsl_vector_set(*startinterval,j,start+(j-*ni)*interval);
         }
         *ni = *ni+niinc;
     }
     
     return (EPOK);
 }
 /*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 4 ************************************************************
  * findIntervalN: This function finds the pulse-free intervals when the input vector has No pulses.
  *                The pulse-free intervals must have a minimum length (intervalMinBins).
  *
  * Parameters:
  * - invector: Input vector WITHOUT pulses
  * - interval: Minimum length of the interval (samples)
  * - ni: Number of pulse-free intervals in the input vector
  * - startinterval: Vector with the starting time of each pulse-free interval (samples)
  *****************************************************************************/
 int findIntervalN (gsl_vector *invector, int interval, int *ni, gsl_vector **startinterval)
 {
     *ni = invector->size/interval;	// Due to both numerator and denominator are integer numbers =>
     // Result is the integer part of the division
     
     // startinterval is allocated with the number of pulse-free intervals
     *startinterval = gsl_vector_alloc(*ni);
     
     for (int i=0;i<*ni;i++)
     {
         gsl_vector_set(*startinterval,i,i*interval);
     }
     
     return (EPOK);
 }
 /*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 5 ************************************************************
  * createTPSreprFile: This function creates the gennoisespec output FITS file (_noisespec.fits).
  *
  * - Create the NOISE representation file (if it does not exist already)
  * - Create the extensions NOISE, NOISEALL and WEIGHTMS
  * - Write keywords
  ****************************************/
 int createTPSreprFile ()
 {
     int status = EPOK, extver=0;
     
     // Create output FITS files: If it does not exist already
     // Create NOISE representation File
     
     if ( fileExists(string(par.outFile)) && par.clobber==1)
     {
         if (remove(par.outFile))
         {
             message = "Output FITS file ("+ string(par.outFile)+") already exits & cannot be deleted ("+string(strerror(errno))+")";
             EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
         }
     }
     else if(fileExists(string(par.outFile)) && par.clobber==0)
     {
         message = "Output FITS file ("+ string(par.outFile)+") already exits: must not be overwritten";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     if(!fileExists(string(par.outFile)))
     {
         if(fits_create_file(&gnoiseObject, par.outFile, &status))
         {
             message = "Cannot create output gennoisespec file " + string(par.outFile);
             EP_PRINT_ERROR(message,status); return(EPFAIL);
         }
     }
     
     message = "Created gennoisespec FITS File (_noisespec): " + string(par.outFile);
     cout<<message<<endl;
     
     // Create extensions: NOISE
     strcpy(extname,"NOISE");
     if (fits_create_tbl(gnoiseObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
     {
         message = "Cannot create table " + string(extname);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     // Create extensions: NOISEALL
     strcpy(extname,"NOISEALL");
     if (fits_create_tbl(gnoiseObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
     {
         message = "Cannot create table " + string(extname);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     if (weightMS == 1)
     {
         // Create extensions: WEIGHTMS
         strcpy(extname,"WEIGHTMS");
         if (fits_create_tbl(gnoiseObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
         {
             message = "Cannot create table " + string(extname);
             EP_PRINT_ERROR(message,status); return(EPFAIL);
         }
     }
     
     // Create keywords in NOISE HDU
     strcpy(extname,"NOISE");
     if (fits_movnam_hdu(gnoiseObject, ANY_HDU,extname, extver, &status))
     {
         message = "Cannot move to HDU " + string(extname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     strcpy(keyname,"EVENTCNT");
     keyvalint = intervalMinBins/2;
     if (fits_write_key(gnoiseObject,TINT,keyname,&keyvalint,comment,&status))
     {
         message = "Cannot write keyword " + string(keyname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     // Create keywords in NOISEall HDU
     strcpy(extname,"NOISEALL");
     if (fits_movnam_hdu(gnoiseObject, ANY_HDU,extname, extver, &status))
     {
         message = "Cannot move to HDU " + string(extname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     strcpy(keyname,"EVENTCNT");
     if (fits_write_key(gnoiseObject,TINT,keyname,&intervalMinBins,comment,&status))
     {
         message = "Cannot write keyword " + string(keyname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     // Primary HDU
     strcpy(extname,"Primary");
     int *hdutype = nullptr;
     if (fits_movabs_hdu(gnoiseObject, 1, hdutype, &status))
     {
         message = "Cannot move to HDU " + string(extname) + " in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     const char *history_text = " ";
     if (fits_write_history(gnoiseObject, history_text, &status)) {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
        EP_PRINT_ERROR(message,status); return(EPFAIL);
     }

     string strhistory = "START PARAMETER list for gennoisespec_" + string(SIRENA_VERSION);
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status)) {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
        EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     int index = 1;
     strhistory = "P" + to_string(index) + " inFile = " + string(par.inFile);
     int num_pieces = strhistory.length()/65+1;    // 65 is a bit less than the length line allowed to write in HISTORY
     string piece_i;
     for (int i=0; i<num_pieces; i++)
     {
         if (i == 0) piece_i = strhistory.substr(0+64*i,64+64*i);
         else        piece_i = strhistory.substr(0+64*i+1,64+64*i);
         
         strcpy(keyvalstr,piece_i.c_str());
         if (fits_write_history(gnoiseObject, piece_i.c_str(), &status))
         {
              message = "Cannot write HISTORY in noise file " + string(par.outFile);
              EP_PRINT_ERROR(message,status); return(EPFAIL);
         }
     }
     index = index +1;
     
     strhistory = "P" + to_string(index) + string(" outFile = ") + string(par.outFile);
     num_pieces = strhistory.length()/65+1;    // 65 is a bit less than the length line allowed to write in HISTORY
     for (int i=0; i<num_pieces; i++)
     {
         if (i == 0) piece_i = strhistory.substr(0+64*i,64+64*i);
         else        piece_i = strhistory.substr(0+64*i+1,64+64*i);
         
         strcpy(keyvalstr,piece_i.c_str());
         if (fits_write_history(gnoiseObject, piece_i.c_str(), &status))
         {
             message = "Cannot write HISTORY in noise file " + string(par.outFile);
             EP_PRINT_ERROR(message,status); return(EPFAIL);
         }
     }
     index = index +1;
     
     char str_intervalMinSamples[125];		sprintf(str_intervalMinSamples,"%d",par.intervalMinSamples);
     strhistory="P" + to_string(index) + string(" intervalMinSamples = ") + string(str_intervalMinSamples);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_nplPF[125];			sprintf(str_nplPF,"%d",par.nplPF);
     strhistory="P" + to_string(index) + string(" nplPF = ") + string(str_nplPF);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_nintervals[125];		sprintf(str_nintervals,"%d",par.nintervals);
     strhistory="P" + to_string(index) + string(" nintervals = ") + string(str_nintervals);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_scaleFactor[125];		sprintf(str_scaleFactor,"%f",par.scaleFactor);
     strhistory="P" + to_string(index) + string(" scaleFactor = ") + string(str_scaleFactor);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_samplesUp[125];		sprintf(str_samplesUp,"%d",par.samplesUp);
     strhistory="P" + to_string(index) + string(" samplesUp = ") + string(str_samplesUp);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_nSgms[125];	    		sprintf(str_nSgms,"%f",par.nSgms);
     strhistory="P" + to_string(index) + string(" nSgms = ") + string(str_nSgms);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_weightMS[125];      sprintf(str_weightMS,"%d",par.weightMS);
     strhistory="P" + to_string(index) + string(" weightMS = ") + string(str_weightMS);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     strhistory="P" + to_string(index) + string(" EnergyMethod = ") + string(par.EnergyMethod);
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_Ifit[125];			sprintf(str_Ifit,"%f",par.Ifit);
     strhistory="P" + to_string(index) + string(" Ifit = ") + string(str_Ifit);
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;

     char str_rmNoiseIntervals[125];      sprintf(str_rmNoiseIntervals,"%d",par.rmNoiseIntervals);
     strhistory="P" + to_string(index) + string(" rmNoiseIntervals = ") + string(str_rmNoiseIntervals);
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     char str_clobber[125];      sprintf(str_clobber,"%d",par.clobber);
     strhistory="P" + to_string(index) + string(" clobber = ") + string(str_clobber);
     strcpy(keyvalstr,strhistory.c_str());
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status))
     {
         message = "Cannot write HISTORY in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     index = index +1;
     
     strhistory = "END PARAMETER list for gennoisespec_" + string(SIRENA_VERSION);
     if (fits_write_history(gnoiseObject, strhistory.c_str(), &status)) {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
        EP_PRINT_ERROR(message,status); return(EPFAIL);
     }

     if (fits_write_history(gnoiseObject, history_text, &status)) {
        message = "Cannot write HISTORY in noise file " + string(par.outFile);
        EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     strcpy(keyname,"CREADATE");
     time_t rawtime;
     struct tm * timeinfo;
     time ( &rawtime );
     timeinfo = localtime ( &rawtime );
     const char * chardate = asctime (timeinfo);  
     strcpy(keyvalstr,chardate);
     if (fits_write_key(gnoiseObject,TSTRING,keyname,keyvalstr,comment,&status))
     {
         message = "Cannot write keyword " + string(keyname) + " in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     strcpy(keyvalstr,SIRENA_VERSION);
     if (fits_write_key(gnoiseObject,TSTRING,"SIRENAV",keyvalstr,comment,&status))
     {
         message = "Cannot write keyword SIRENAV in noise file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     return EPOK;
 }
 /*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 6 ************************************************************
  * writeTPSreprExten: This function writes the noisespec output FITS file (_noisespec.fits).
  *
  * - Allocate GSL vectors
  * - Write the data in the output FITS file (print only half of FFT to prevent aliasing)
  * - NOISE HDU only contains positive frequencies (=>Multiply by 2 the amplitude)
  * - NOISEALL HDU contains negative and positive frequencies => It is the HDU read to build the optimal filters
  * - WEIGHTMS HDU
  *****************************************************************************/
 int writeTPSreprExten ()
 {	
     int status = EPOK;
     double SelectedTimeDuration;
     
     // Allocate GSL vectors
     gsl_vector *local_freqgsl = gsl_vector_alloc(intervalMinBins/2);
     gsl_vector *local_csdgsl = gsl_vector_alloc(intervalMinBins/2);
     gsl_vector *sigmacsdgslnew = gsl_vector_alloc(intervalMinBins/2);
     gsl_vector *freqALLgsl = gsl_vector_alloc(intervalMinBins);
     gsl_vector *csdALLgsl = gsl_vector_alloc(intervalMinBins);
     
     SelectedTimeDuration=intervalMinBins/((double)samprate);
     
     // Print only half of FFT to prevent aliasing
     for (int i=0; i< (intervalMinBins/2); i++)
     {
         gsl_vector_set(local_freqgsl,i,i/SelectedTimeDuration);
         gsl_vector_set(local_csdgsl,i,2*gsl_vector_get(EventSamplesFFTMean,i)); // *2 if only positive frequencies
         gsl_vector_set(sigmacsdgslnew,i,gsl_vector_get(sigmacsdgsl,i));
         gsl_vector_set(freqALLgsl,i,i/SelectedTimeDuration);
         gsl_vector_set(csdALLgsl,i,gsl_vector_get(EventSamplesFFTMean,i));
     }
     gsl_vector_set(freqALLgsl,intervalMinBins/2,(intervalMinBins/2)/SelectedTimeDuration);
     gsl_vector_set(csdALLgsl,intervalMinBins/2,gsl_vector_get(EventSamplesFFTMean,intervalMinBins/2));
     for (double i=1; i<(intervalMinBins/2); i++)
     {
         gsl_vector_set(freqALLgsl,i+intervalMinBins/2,-(i+intervalMinBins/2-i*2)/SelectedTimeDuration);
         gsl_vector_set(csdALLgsl,i+intervalMinBins/2,gsl_vector_get(EventSamplesFFTMean,i+intervalMinBins/2));
     }
     
     obj.inObject = gnoiseObject;
     obj.nameTable = new char [255];
     strcpy(obj.nameTable,"NOISE");
     obj.iniRow = 1;
     obj.endRow = intervalMinBins/2;
     obj.iniCol = 0;
     obj.nameCol = new char [255];
     strcpy(obj.nameCol,"FREQ");
     obj.type = TDOUBLE;
     obj.unit = new char [255];
     strcpy(obj.unit,"Hz");
     if (writeFitsSimple(obj, local_freqgsl))
     {
         message = "Cannot run routine writeFitsSimple for freqgsl";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     gsl_vector_free(local_freqgsl); local_freqgsl = 0;
     
     obj.inObject = gnoiseObject;
     obj.nameTable = new char [255];
     strcpy(obj.nameTable,"NOISE");
     obj.iniRow = 1;
     obj.endRow = intervalMinBins/2;
     obj.iniCol = 0;
     obj.nameCol = new char [255];
     strcpy(obj.nameCol,"CSD");
     obj.type = TDOUBLE;
     obj.unit = new char [255];
     if (strcmp(par.EnergyMethod,"OPTFILT") == 0)
     {
         strcpy(obj.unit,"A/sqrt(Hz)");
     }
     else
     {
         strcpy(obj.unit,"adu/sqrt(Hz)");
     }
     //strcpy(obj.unit,"A/sqrt(Hz)");
     if (writeFitsSimple(obj, local_csdgsl))
     {	
         message = "Cannot run routine writeFitsSimple for csdgsl";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     gsl_vector_free(local_csdgsl); local_csdgsl = 0;
     
     obj.inObject = gnoiseObject;	
     obj.nameTable = new char [255];
     strcpy(obj.nameTable,"NOISE");
     obj.iniRow = 1;
     obj.endRow = intervalMinBins/2;
     obj.iniCol = 0;
     obj.nameCol = new char [255];
     strcpy(obj.nameCol,"SIGMACSD");
     obj.type = TDOUBLE;	
     obj.unit = new char [255];
     if (strcmp(par.EnergyMethod,"OPTFILT") == 0)
     {
         strcpy(obj.unit,"A/sqrt(Hz)");
     }
     else
     {
         strcpy(obj.unit,"adu/sqrt(Hz)");
     }
     if (writeFitsSimple(obj, sigmacsdgslnew))
     {
         message = "Cannot run routine writeFitsSimple for sigmacsdgsl";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     gsl_vector_free(sigmacsdgslnew); sigmacsdgslnew = 0;
     
     strcpy(keyname,"NOISEBSL");       // Real calculated baseline
     double sumBaseline;
     gsl_vector_Sumsubvector(baseline, 0, indexBaseline, &sumBaseline);
     double keyvaldouble;
     if (strcmp(par.EnergyMethod,"OPTFILT") == 0)
     {
        keyvaldouble = (sumBaseline/indexBaseline)/aducnv;
     }
     else
     {
         keyvaldouble = (sumBaseline/indexBaseline);
     }
     if (fits_write_key(gnoiseObject,TDOUBLE,keyname,&keyvaldouble,comment,&status))
     {
         message = "Cannot write keyword " + string(keyname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     strcpy(keyname,"NOISESTD");
     double sumSigma;
     gsl_vector_Sumsubvector(sigma, 0, indexBaseline, &sumSigma);
     keyvaldouble = (sumSigma/indexBaseline)/aducnv;
     if (fits_write_key(gnoiseObject,TDOUBLE,keyname,&keyvaldouble,comment,&status))
     {
         message = "Cannot write keyword " + string(keyname) + " in file " + string(par.outFile);
         EP_PRINT_ERROR(message,status); return(EPFAIL);
     }
     
     obj.inObject = gnoiseObject;		
     obj.nameTable = new char [255];
     strcpy(obj.nameTable,"NOISEALL");
     obj.iniRow = 1;
     obj.endRow = intervalMinBins;
     obj.iniCol = 0;
     obj.nameCol = new char [255];
     strcpy(obj.nameCol,"FREQ");
     obj.type = TDOUBLE;
     obj.unit = new char [255];
     strcpy(obj.unit,"Hz");
     if (writeFitsSimple(obj, freqALLgsl))
     {
         message = "Cannot run routine writeFitsSimple for freqALLgsl";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     gsl_vector_free(freqALLgsl); freqALLgsl = 0;
     
     strcpy(obj.nameCol,"CSD");
     if (strcmp(par.EnergyMethod,"OPTFILT") == 0)
     {
        strcpy(obj.unit,"A/sqrt(Hz)");
     }
     else
     {
         strcpy(obj.unit,"adu/sqrt(Hz)");
     }
     if(writeFitsSimple(obj, csdALLgsl))
     {
         message = "Cannot run routine writeFitsSimple for csdALLgsl";
         EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
     }
     gsl_vector_free(csdALLgsl); csdALLgsl = 0;
     
     if (weightMS == 1)
     {
         obj.inObject = gnoiseObject;		
         obj.nameTable = new char [255];
         strcpy(obj.nameTable,"WEIGHTMS");
         obj.iniRow = 1;
         obj.endRow = 1;
         obj.iniCol = 0;
         obj.nameCol = new char [255];
         char str_length[125];
         gsl_vector *weightMatrixesrow= gsl_vector_alloc(weightMatrixes->size2);
         gsl_vector_view temp;
         
         // All the Wx (being x 'intervalMinSamples') columns appear in the WEIGHTMS HDU.
         for (int i=0; i<(int)(weightpoints->size);i++)
         {
             snprintf(str_length,125,"%d",(int) gsl_vector_get(weightpoints,i));
             strcpy(obj.nameCol,(string("W")+string(str_length)).c_str());
             obj.type = TDOUBLE;
             obj.unit = new char [255];
             strcpy(obj.unit," ");
             gsl_matrix_get_row(weightMatrixesrow,weightMatrixes,i);
             gsl_matrix *weightMatrixes_matrix = gsl_matrix_alloc(1,gsl_vector_get(weightpoints,i)*gsl_vector_get(weightpoints,i));
             temp = gsl_vector_subvector(weightMatrixesrow,0,gsl_vector_get(weightpoints,i)*gsl_vector_get(weightpoints,i));
             gsl_matrix_set_row(weightMatrixes_matrix,0,&temp.vector);
             if (writeFitsComplex(obj,weightMatrixes_matrix))
             {
                 message = "Cannot run routine writeFitsSimple for weightMatrixes_matrix";
                 EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
             }
             gsl_matrix_free(weightMatrixes_matrix); 
         }
         gsl_vector_free(weightMatrixesrow);
     }
     
     return (EPOK);
 }
 /*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 7 ************************************************************
  * findPulsesNoise function: This function is going to find the pulses in a record by using the function findTstartNoise 
  * 
  * - Declare variables
  * - Establish the threshold (call medianKappaClipping)
  * - Find pulses (call findTstartNoise) 
  * - Free allocated GSL vectors
  *
  * Parameters:
  * - vectorinDER: Derivative of the low-pass filtered 'vectorin'
  * - tstart: Starting time of the found pulses into the record (in samples)
  * - quality: Quality of the found pulses into the record
  * - nPulses: Number of found pulses
  * - threshold: Threshold used to find the pulses (output parameter because it is necessary out of the function)
  * - scalefactor: Scale factor to calculate the LPF box-car length
  * - samplingRate: Sampling rate
  * - samplesup: Number of consecutive samples over the threshold to locate a pulse ('samplesUp')
  * - nsgms: Number of Sigmas to establish the threshold
  * - stopCriteriamkc: Used in medianKappaClipping_noiseSigma (%)
  * - kappamkc: Used in medianKappaClipping_noiseSigma
  * 
  ****************************************************************************/
 int findPulsesNoise
 (
  gsl_vector *vectorinDER,
  gsl_vector **tstart,
  gsl_vector **quality,
  
  int *nPulses,
  double *threshold,
  
  double scalefactor,
  double samplingRate,
  
  int samplesup,
  double nsgms,
  
  double local_stopcriteriamkc,
  double kappamkc)
 {
     // Declare variables
     double thresholdmediankappa;	// Threshold to look for pulses in the first derivative
     
     gsl_vector *maxDERgsl = gsl_vector_alloc(vectorinDER->size);
      
     gsl_vector_set_zero(*quality);
     
     // First step to look for single pulses: Establish the threshold (call medianKappaClipping)
     if (medianKappaClipping (vectorinDER, kappamkc, local_stopcriteriamkc, nsgms, (int)(pi*samplingRate*scalefactor), &thresholdmediankappa))
     {
         message = "Cannot run medianKappaClipping looking for single pulses";
         EP_PRINT_ERROR(message,EPFAIL);
     }
     *threshold = thresholdmediankappa;
     
     if (findTstartNoise (10000,vectorinDER, thresholdmediankappa, samplesup, nPulses, tstart, quality, &maxDERgsl))
     {
         message = "Cannot run findTstartNoise with two rows in models";
         EP_PRINT_ERROR(message,EPFAIL);
     }
       
     // Free allocated GSL vectors
     if (maxDERgsl != NULL)      {gsl_vector_free(maxDERgsl); maxDERgsl = 0;}
     
     return(EPOK);
 }
 /*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 8 ************************************************************
  * findTstartNoise function: This function finds the pulses tstarts in the input vector (first derivative of the filtered record)
  *
  * This function scans all values the derivative of the (low-pass filtered) record until it finds nSamplesUp consecutive values (due to the noise more than 1 value is
  * required) over the threshold. To look for more pulses, it finds nSamplesUp consecutive values
  * (due to the noise) under the threshold and then, it starts to scan again.
  *
  * - Declare variables
  * - Allocate GSL vectors
  * - Obtain tstart of each pulse in the derivative:
  * 	- If der_i>threshold and foundPulse=false, it looks for nSamplesUp consecutive samples over the threshold
  *   		- If not, it looks again for a pulse crossing over the threshold
  *	        - If yes, a pulse is found (truncated if it is at the beginning)
  *	- If der_i>threshold and foundPulse=true, it looks for a sample under the threshold
  *   		- If not, it looks again for a sample under the threshold
  *       	- If yes, it looks for nSamplesUp consecutive samples under the threshold and again it starts to look for a pulse
  *
  * Parameters:
  * - maxPulsesPerRecord: Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
  * - der: First derivative of the (low-pass filtered) record
  * - adaptativethreshold: Threshold
  * - nSamplesUp: Number of consecutive samples over the threshold to 'find' a pulse
  * - numberPulses: Number of found pulses
  * - tstartgsl: Pulses tstart (in samples)
  * - flagTruncated: Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
  * - maxDERgsl: Maximum of the first derivative of the (low-pass filtered) record inside each found pulse
  ******************************************************************************/
 int findTstartNoise
 (
    int maxPulsesPerRecord,
  
    gsl_vector *der,
    double adaptativethreshold,
    int nSamplesUp,

    int *numberPulses,

    gsl_vector **tstartgsl,
    gsl_vector **flagTruncated,
    gsl_vector **maxDERgsl)
 {
     char valERROR[256];
     
     // Declare variables
     int szRw = der->size;	 // Size of segment to process
     *numberPulses = 0;
     bool foundPulse = false;
     int i = 0;		// To go through the elements of a vector
     
     // Allocate GSL vectors
     // It is not necessary to check the allocation because 'maxPulsesPerRecord'='EventListSize'(input parameter) must already be > 0
     *tstartgsl = gsl_vector_alloc(maxPulsesPerRecord);	
     *flagTruncated = gsl_vector_alloc(maxPulsesPerRecord);
     gsl_vector_set_zero(*flagTruncated);
     *maxDERgsl = gsl_vector_alloc(maxPulsesPerRecord);	// Maximum of the first derivative
     gsl_vector_set_all(*maxDERgsl,-1E3);
     
     int cntUp = 0;
     int cntDown = 0;
     double possibleTstart = 0.0;
     double possiblemaxDER = 0.0;
     
     
     // It looks for a pulse
     // If a pulse is found (foundPulse==true) => It looks for another pulse
     do
     {
         foundPulse = false;
         
         // It looks for a pulse since the beginning (or the previous pulse) to the end of the record
         while (i < szRw-1)
         {
             if (foundPulse == false)
             {
                 // The first condition to detect a pulse is that the adjustedDerivative was over the threshold
                 if (gsl_vector_get(der,i) > adaptativethreshold)
                 {
                     cntUp++;
                     if (cntUp == 1)
                     {
                         possibleTstart = i;
                         possiblemaxDER = gsl_vector_get(der,i);
                     }
                     else if (cntUp == nSamplesUp)
                     {
                         if (*numberPulses == maxPulsesPerRecord)
                         {
                             sprintf(valERROR,"%d",__LINE__+6);
                             string str(valERROR);
                             message = "Found pulses in record>'EventListSize'(input parameter) => Change EventListSize or check if the threshold is too low => Setting i-th element of vector out of range in line " + str + " (" + __FILE__ + ")"; 
                             str.clear();
                             EP_PRINT_ERROR(message,EPFAIL);
                         }
                         gsl_vector_set(*maxDERgsl,*numberPulses,possiblemaxDER);
                         gsl_vector_set(*tstartgsl,*numberPulses,possibleTstart);
                         if (possibleTstart == 0)	gsl_vector_set(*flagTruncated,*numberPulses,1);
                         *numberPulses = *numberPulses +1;
                         foundPulse = true;
                         cntUp = 0;
                     }
                     cntDown = 0;
                 }
                 else cntUp = 0;
             }
             else
             {
                 if (gsl_vector_get(der,i) > gsl_vector_get(*maxDERgsl,*numberPulses-1))
                 {
                     gsl_vector_set(*maxDERgsl,*numberPulses-1,gsl_vector_get(der,i));
                 }
                 else if (gsl_vector_get(der,i) < adaptativethreshold)
                 {
                     cntUp = 0;
                     cntDown++;
                     if (cntDown == nSamplesUp) // nSamplesUp samples under the threshold in order to look for another pulse
                     {
                         foundPulse = false; 
                         cntDown = 0;
                     }
                     
                 }
                 else if (gsl_vector_get(der,i) > adaptativethreshold)
                 {
                     cntDown = 0;
                 }
             }
             i++;
         }
     } while (foundPulse == true);
     
     return (EPOK);
 }
 /*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 9 ************************************************************
  * weightMatrixNoise function: This function calculates the weight matrix of the noise.
  *
  * Di: Pulse free interval
  * V: Covariance matrix
  * 
  *  Vij = E[DiDj]-E[Di]E[Dj] 
  * 
  * Di^p: Value of the pth-sample of the pulse-free interval i
  * N: Number of samples
  *
  *       (i)
  * <DiDj> = E[DiDj] - E[Di]E[Dj] =
  *       (ii)
  *        =(1/N)sum(p=1,N){(Di^p)(Dj^p)} - (1/N)sum(p=1,N){(Di^p) * (1/N)sum(p=1,N){(Dj^p)
  *
  * (i) Var(X) = <X> = E[(X-E[X])^2] = E[X^2] - (E[X])^2
  * (ii) E[X] = (1/N)sum(i=1,N){(xi}
  * 
  * 		        |<D1D1> <D1D2>...<D1Dn>|
  *  Vij = <DiDj> = |<D2D1> <D2D2>...<D2Dn>|	where n is the pulse-free interval length     V => (nxn)
  *                 |...                   |
  *                 |<DnD1> <DnD2>...<DnDn>|
  *  W = 1/V
  *
  * - Calculate the elements of the diagonal of the covariance matrix
  * - Calculate the elements out of the diagonal of the covariance matrix
  * - Calculate the weight matrix
  *
  * Parameters:
  * - intervalMatrix: GSL matrix containing pulse-free intervals whose baseline is 0 (baseline previously subtracted) [nintervals x intervalMinSamples]
  * - weight: GSL matrix with weight matrix
  ******************************************************************************/
 int weightMatrixNoise (gsl_matrix *intervalMatrix, gsl_matrix **weight)
 {
     char valERROR[256];
     
     double elementValue;
     double elementValue1 = 0.0;
     double elementValue2 = 0.0;
     double elementValue3 = 0.0;
     gsl_permutation *perm = gsl_permutation_alloc((*weight)->size1);
     int s=0;
     
     gsl_matrix *covariance = gsl_matrix_alloc((*weight)->size1,(*weight)->size2);
     
     clock_t t;
     t=clock();
     // Elements of the diagonal of the covariance matrix
     for (int i=0;i<(int)(intervalMatrix->size2);i++)
     {
         for (int p=0;p<(int)(intervalMatrix->size1);p++)
         {
             elementValue1 = elementValue1 + pow(gsl_matrix_get(intervalMatrix,p,i),2.0);
             elementValue2 = elementValue2 + gsl_matrix_get(intervalMatrix,p,i);
         }
         
         elementValue1 = elementValue1/intervalMatrix->size1;
         elementValue2 = elementValue2/intervalMatrix->size1;
         elementValue = elementValue1-elementValue2*elementValue2;
         
         gsl_matrix_set(covariance,i,i,elementValue);
         
         elementValue = 0.0;
         elementValue1 = 0.0;
         elementValue2 = 0.0;
     }
     //cout<<"Matrix diagonal ended "<<covariance->size1<<"x"<<covariance->size2<<endl;
     t = clock() - t;
     //cout<<"Consumed "<<((float)t)/CLOCKS_PER_SEC<<" sec"<<endl;
     
     t = clock();
     // Other elements
     for (int i=0;i<(int)(intervalMatrix->size2);i++)
     {
         for (int j=i+1;j<(int)(intervalMatrix->size2);j++)
         {
             for (int p=0;p<(int)(intervalMatrix->size1);p++)
             {
                 elementValue1 = elementValue1 + (gsl_matrix_get(intervalMatrix,p,i)*gsl_matrix_get(intervalMatrix,p,j));		
                 elementValue2 = elementValue2 + gsl_matrix_get(intervalMatrix,p,i);
                 elementValue3 = elementValue3 + gsl_matrix_get(intervalMatrix,p,j);
             }
             
             elementValue1 = elementValue1/intervalMatrix->size1;
             elementValue2 = elementValue2/intervalMatrix->size1;
             elementValue3 = elementValue3/intervalMatrix->size1;
             elementValue = elementValue1-elementValue2*elementValue3;
             
             gsl_matrix_set(covariance,i,j,elementValue);
             gsl_matrix_set(covariance,j,i,elementValue);
             
             elementValue = 0.0;
             elementValue1 = 0.0;
             elementValue2 = 0.0;
             elementValue3 = 0.0;
         }
     }
     
     //cout<<"Elements out of the matrix diagonal ended "<<covariance->size1<<"x"<<covariance->size2<<endl;
     t = clock() - t;
     //cout<<"Consumed "<<((float)t)/CLOCKS_PER_SEC<<" sec"<<endl;
     
     t = clock();
     // Calculate the weight matrix
     // It is not necessary to check the allocation because 'covarianze' size must already be > 0
     gsl_matrix *covarianceaux = gsl_matrix_alloc(covariance->size1,covariance->size2);
     gsl_matrix_memcpy(covarianceaux,covariance);
     //cout<<"Preparation to the inversion ended "<<covariance->size1<<"x"<<covariance->size2<<endl;
     t = clock() - t;
     //cout<<"Consumed "<<((float)t)/CLOCKS_PER_SEC<<" sec"<<endl;
     t = clock();
     gsl_linalg_LU_decomp(covarianceaux, perm, &s);
     if (gsl_linalg_LU_invert(covariance, perm, *weight) != 0) 
     {
         sprintf(valERROR,"%d",__LINE__-2);
         string str(valERROR);
         message = "Singular matrix in line " + str + " (" + __FILE__ + ")";
         str.clear();
         EP_PRINT_ERROR(message,EPFAIL);	return(EPFAIL);
     }
     
     gsl_matrix_free(covarianceaux); covarianceaux=0;
     gsl_matrix_free(covariance); covariance=0;
     gsl_permutation_free(perm); perm=0;
     
     return (EPOK);
 }
 /*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 10 ************************************************************
  * medianKappaClipping_noiseSigma function: This function provides the mean and the sigma of an input vector (with noise sigmas) 
  *                               by using a Kappa-clipping method (replacing points beyond mean+-Kappa*sigma with the
  *                               median).
  *
  * First, mean and sigma are calculated and invector values out of (mean+Kappa*sigma,mean-Kappa*sigma) are replaced
  * with the median (it is trying to look for the baseline). And this process is iteratively repeated until there are
  * no points beyond mean+-Kappa *sigma. Finally, the mean and sigma of the resulting vector are provided.
  *
  * - Declare variables
  * - Calculate the median
  * - Iterate until there are no points out of the maximum excursion (kappa*sigma)
  * - Calculate mean and sigma
  *
  * Parameters:
  * - invector: First derivative of the (filtered) record
  * - kappa: To establish the range around of the mean
  * - stopCriteria: It is given in %
  * - mean: Mean value of the invector (no points beyond mean+-Kappa *sigma)
  * - local_sigma: Sigma value of the invector (no points beyond mean+-Kappa *sigma)
  ******************************************************************************/
 int medianKappaClipping_noiseSigma (gsl_vector *invector, double kappa, double stopCriteria, double *mean, double *local_sigma)
 {
     char valERROR[256];
     
     // Declare variables
     int size = invector->size; // Size of the input vector
     double mean1, sg1;
     double mean2, sg2;
     gsl_vector_view temp;
     // Variables to remove input vector elements higher than the maximum excursion (kappa*sg)
     int i;							// To go through the elements of a vector
     int cnt;						// Number of points inside the excursion (mean+-excursion)
     // It is not necessary to check the allocation because 'invector' size must already be > 0
     gsl_vector *invectorNew = gsl_vector_alloc(size);	// Auxiliary vector
     
     double median;
     int boxLPF = 0;
     
     gsl_vector_memcpy(invectorNew,invector);
     gsl_sort_vector(invectorNew);
     if (size%2 == 0)	//Even
     {
         median = (gsl_vector_get(invectorNew,(int) (size/2)-1)+gsl_vector_get(invectorNew,(int) (size/2)))/2;
     }
     else                    //Odd
     {
         median = gsl_vector_get(invectorNew,(int) (size/2));
     }
     
     gsl_vector_memcpy(invectorNew,invector);
     
     // Iterate until no points out of the maximum excursion (kappa*sigma)
     do
     {
         if ((size-boxLPF-1 < 1) || (size-boxLPF-1 >(int) invectorNew->size))
         {
             sprintf(valERROR,"%d",__LINE__+5);
             string str(valERROR);
             message = "View goes out of scope the original vector in line " + str + " (" + __FILE__ + ")";  str.clear();
             EP_PRINT_ERROR(message,EPFAIL);
         }
         temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
         if (findMeanSigma (&temp.vector, &mean1, &sg1))
         {
             message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
             EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
         }
         i = 0;
         cnt = 0;
         
         while (i<((int)invectorNew->size))
         {
             if ((gsl_vector_get(invectorNew,i) >= mean1 + kappa*sg1) || (gsl_vector_get(invectorNew,i) <= mean1 - kappa*sg1))
                 // HARDPOINT!!!!!!!!!!!!!!!!!!! (kappa)
             {
                 if ((i < 0) || (i >((int)(invectorNew)->size)-1))
                 {
                     sprintf(valERROR,"%d",__LINE__+5);
                     string str(valERROR);
                     message = "Setting with <= 0 size in line " + str + " (" + __FILE__ + ")"; str.clear();
                     EP_PRINT_ERROR(message,EPFAIL);
                 }
                 gsl_vector_set(invectorNew,i,median);
                 cnt++;
             }
             i++;
         }
         
         if (cnt != 0)
             // Some points of the invector have been replaced with the median
         {
             temp = gsl_vector_subvector(invectorNew,0,size-boxLPF-1);
             if (findMeanSigma (&temp.vector, &mean2, &sg2))
             {
                 message = "Cannot run findMeanSigma routine for kappa-sigma iteration after replacement with the median";
                 EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
             }
         }
         else
             // No points of the invector have been replaced with the median
         {
             mean2 =mean1;
             sg2 = sg1;
         }
         
     } while (fabs((mean2-mean1)/mean1)>(stopCriteria/100.0));	// HARDPOINT!!!!!!!!!!!!!!!!!!! (stopCriteria)
     
     *mean = mean2;
     *local_sigma =sg2;
     
     gsl_vector_free(invectorNew); invectorNew= 0;
     
     return EPOK;
 }
 /*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 11 ************************************************************
  * getpar_noiseSpec function: This function gets the input parameter from the command line or their default values from the gennoisespec.par file
  *
  * Parameters:
  * - local_par: Structure containing the input parameters
  ******************************************************************************/
 int getpar_noiseSpec(struct Parameters* const local_par)
 {
     // String input buffer.
     char* sbuffer=NULL;
     
     // Error status.
     int status=EXIT_SUCCESS;
     
     status=ape_trad_query_string("inFile", &sbuffer);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the input FITS file name";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     strcpy(local_par->inFile, sbuffer);
     free(sbuffer);
     
     status=ape_trad_query_string("outFile", &sbuffer);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the output FITS file name";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     strcpy(local_par->outFile, sbuffer);
     free(sbuffer);
     
     status=ape_trad_query_int("intervalMinSamples", &local_par->intervalMinSamples);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the intervalMinSamples parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     if ((local_par->intervalMinSamples < 2) || ((local_par->intervalMinSamples)%2 != 0))
     {
         message = "intervalMinSamples must be even and greater than 0";
         return(EXIT_FAILURE);
     }
     
     status=ape_trad_query_int("nplPF", &local_par->nplPF);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the nplPF parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     status=ape_trad_query_int("nintervals", &local_par->nintervals);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the nintervals parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     if (local_par->nintervals == 0) MyAssert(local_par->nintervals > 0, (char *)("nintervals must be greater than 0"));
     
     status=ape_trad_query_double("scaleFactor", &local_par->scaleFactor);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the scaleFactor parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     status=ape_trad_query_int("samplesUp", &local_par->samplesUp);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the samplesUp parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     if (local_par->samplesUp == 0) MyAssert(local_par->samplesUp > 0, (char *)("samplesUp must be greater than 0"));
     
     status=ape_trad_query_double("nSgms", &local_par->nSgms);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the nSgms parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     if (local_par->nSgms == 0) MyAssert(local_par->nSgms >= 1, (char *)("nSgms must be greater than 1"));
     
     status=ape_trad_query_bool("weightMS", &local_par->weightMS);
     if (local_par->weightMS == 1)    weightMS = 1;
     
     status=ape_trad_query_string("EnergyMethod", &sbuffer);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the EnergyMethod parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     strcpy(local_par->EnergyMethod, sbuffer);
     free(sbuffer);
     
     MyAssert((strcmp(local_par->EnergyMethod,"OPTFILT") == 0) || (strcmp(local_par->EnergyMethod,"I2R") == 0) ||	(strcmp(local_par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(local_par->EnergyMethod,"I2RDER") == 0),
              (char *)("EnergyMethod must be OPTFILT, I2R, I2RFITTED, I2RDER or PCA"));
     
     status=ape_trad_query_double("Ifit", &local_par->Ifit);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the Ifit parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     status=ape_trad_query_bool("clobber", &local_par->clobber);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the clobber parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
    status=ape_trad_query_bool("rmNoiseIntervals", &local_par->rmNoiseIntervals);
     if (EXIT_SUCCESS!=status) {
         message = "failed reading the rmNoiseIntervals parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     if (EXIT_SUCCESS!=status) 
     {
         message = "Failed reading some SIRENA parameter";
         EP_EXIT_ERROR(message,EPFAIL);
     }
     
     return(status);
 }
 /*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
 
 /***** SECTION 12 ************************************************************
  * MyAssert function: This function displays an error message if the condition in 'expr' is true
  *
  * Parameters:
  * - expr: Condition to be true in order to display the error message
  * - msg: Message to be displayed
  ******************************************************************************/
 void MyAssert(int expr, char* msg)
 {
     if (expr == 0)
     {
         printf("%s %s %s"," Assertion failure: ",msg,"\n");
         abort();
     }
 }
 /*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 
