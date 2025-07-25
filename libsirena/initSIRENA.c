/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      INITSIRENA
*
*  File:       initSIRENA.c
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*
***********************************************************************/

#include "initSIRENA.h"
#include <stdlib.h>

/***** SECTION 1 ************************************************************
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
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* checkxmls function: Check if the XML file used to build the library is the same to be used to recconstruct (by checking the checksums)
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int checkXmls(struct Parameters* const par)
{
    // Error status.
    int status=EXIT_SUCCESS;

    // Read XMLCKSM from library file or calculate it
    fitsfile* libptr = NULL;
    // Move to "Primary" HDU of the library file
    fits_open_file(&libptr, par->LibraryFile, READONLY, &status);
    if (status != 0)
    {
        SIXT_ERROR("File given in LibraryFile does not exist");
        return(EXIT_FAILURE);
    }
    if (fits_movabs_hdu(libptr, 1, NULL, &status))
    {
        return(EXIT_FAILURE);
    }
    // Read XMLCKSM from library
    int keywordvalue;
    unsigned checksum_libXMLfile;
    char libXMLfile[1024] = "x";
    fits_read_key(libptr, TINT, "XMLCKSM", &keywordvalue, NULL, &status);
    if (status == EXIT_SUCCESS)
    {
        checksum_libXMLfile = keywordvalue;
    }
    else // or read full Primary HDU and store it in 'libheaderPrimary' and calculate XML checksum
    {
        status = EXIT_SUCCESS;
        int numberkeywords;
        char *libheaderPrimary = NULL;
        if (fits_hdr2str(libptr, 0, NULL, 0,&libheaderPrimary, &numberkeywords, &status))
        {
            free(libheaderPrimary);
            return(EXIT_FAILURE);
        }

        char *xml_pointer = NULL;
        char *XMLFile_pointer = NULL;
        char wholelibXMLfile[1024] = "x";
        xml_pointer = strstr(libheaderPrimary,".xml");
        if(!xml_pointer)
        {
            SIXT_ERROR("XML file info not included in Primary HDU in library file");
            return(EXIT_FAILURE);
        }
        XMLFile_pointer = strstr(libheaderPrimary,"XMLFile = ");
        if(!XMLFile_pointer)
        {
            SIXT_ERROR("XML file info not included in Primary HDU in library file");
            return(EXIT_FAILURE);
        }
        strncpy(wholelibXMLfile,XMLFile_pointer+10,xml_pointer-XMLFile_pointer-10+4);  // 10 -> "XMLFile = ", 4 -> ".xml"

        int HISTORYnum = 0;
        int k;
        int lengthTOTAL = (int)strlen(wholelibXMLfile);
        int spaces[lengthTOTAL];
        k=0;
        do
        {
            if(wholelibXMLfile[k]==' ')
            {
                if (HISTORYnum % 2 == 0)
                {  // Even
                    spaces[HISTORYnum] = k-7;
                }
                else
                {
                    spaces[HISTORYnum] = k;
                }

                HISTORYnum++;
            }
            k++;
        }while(k<=lengthTOTAL);

        if (HISTORYnum != 0) HISTORYnum = HISTORYnum/2;

        char libXMLfile2[1024] = "x";

        if (HISTORYnum == 0)
        {
            strcpy(libXMLfile,wholelibXMLfile);
        }
        else
        {
            for (int i=0;i<=HISTORYnum;i++)
            {
                if (i == 0) // First
                {
                    subString (wholelibXMLfile, 0, spaces[i], libXMLfile);
                }
                else
                {
                    if (i == HISTORYnum) // Last
                    {
                        subString (wholelibXMLfile, spaces[2*i-1]+1, lengthTOTAL-spaces[2*i-1]-1, libXMLfile2);
                        strcat(libXMLfile,libXMLfile2);
                        memset(libXMLfile2,0,1024);
                    }
                    else
                    {
                        subString (wholelibXMLfile, spaces[2*i-1]+1, spaces[2*i]-spaces[2*i-1]-1, libXMLfile2);
                        strcat(libXMLfile,libXMLfile2);
                        memset(libXMLfile2,0,1024);
                    }
                }
            }
        }

        if (libheaderPrimary != NULL)
        {
            free(libheaderPrimary);
        }

        FILE *fp_libXMLfile;
        size_t len_libXMLfile;
        char buf_libXMLfile[4096];
        if (NULL == (fp_libXMLfile = fopen(libXMLfile, "rb")))
        {
            printf("Unable to open XML from library, %s, for reading it and calculate its checksum.\n", libXMLfile);
            return -1;
        }
        len_libXMLfile = fread(buf_libXMLfile, sizeof(char), sizeof(buf_libXMLfile), fp_libXMLfile);
        checksum_libXMLfile = checksum(buf_libXMLfile, len_libXMLfile, 0);
    }
    fits_close_file(libptr,&status);

    // Calculate checksum from the XML input FITS file
    FILE *fp_reconsXMLfile;
    size_t len_reconsXMLfile;
    unsigned checksum_reconsXMLfile;
    char reconsXMLfile[1024] = "x";
    char buf_reconsXMLfile[4096];
    strcpy(reconsXMLfile,par->XMLFile);
    if (NULL == (fp_reconsXMLfile = fopen(reconsXMLfile, "rb")))
    {
          printf("Unable to open provided XML, %s, for reading it and calculate its checksum\n", reconsXMLfile);
          return -1;
    }
    len_reconsXMLfile = fread(buf_reconsXMLfile, sizeof(char), sizeof(buf_reconsXMLfile), fp_reconsXMLfile);
    checksum_reconsXMLfile = checksum(buf_reconsXMLfile, len_reconsXMLfile, 0);

    // Compare checksums
    if (checksum_libXMLfile == checksum_reconsXMLfile) status = 0;
    else status = 1;

    if (status != 0)
    {
        SIXT_ERROR("XML file from library FITS file and from input parameter do not match (checksum)");
        printf("The checksum of XML from library, %s, is %#x\n", libXMLfile, checksum_libXMLfile);
        printf("The checksum of provided XML, %s, is %#x\n", reconsXMLfile, checksum_reconsXMLfile);
        return(EXIT_FAILURE);
    }

    return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* subString: This function extracts some elements from an array of characters
*
* Parameters:
* - input: Array of characters from which some elements are extracted
* - offset: Offset
* - len: Length (number of elements to extract)
* - dest: Array of characters into which the extracted characters are written
******************************************************************************/
char* subString (const char* input, int offset, int len, char* dest)
{
  int input_len = strlen (input);

  if (offset + len > input_len)
  {
     return NULL;
  }

  strncpy (dest, input + offset, len);
  return dest;
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* getSamplingrate_trigreclength_Filei: This function gets the sampling rate and the trig_reclength
*
* Steps:
* - Open FITS file
* - Check if input FITS file have been simulated with TESSIM or XIFUSIM
* - Check if input XML file and XMl file to build the library to be used to
*   reconstruct are the same
* - Get the sampling rate from the HISTORY keyword from the input FITS file
*   and check with sampling rate from XML file
* - If xifusim file => Get 'trig_reclength' from the HISTORY keyword from the
*   input FITS file 'trig_reclength' is necessary if SIRENA is going to run in
*   THREADING mode
* - Close FITS file
* - Free memory
*
* Parameters:
* - inputFile: Input file name
* - par: Input parameters
* - samplingrate: (In) Sampling rate from XML file => (Out) Sampling rate
* - trigreclength: Necessary if SIRENA is going to run in THREADING mode
******************************************************************************/
int getSamplingrate_trigreclength_Filei (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength)
{
    double samplingrate_XML = *samplingrate; // To work with a more meaningful name

    // Error status.
    int status=EXIT_SUCCESS;

    // Open file
    fitsfile* fptr = NULL;
    fits_open_file(&fptr, inputFile, READONLY, &status);
    if (status != 0)
    {
        SIXT_ERROR("File given in RecordFile does not exist");
        return(EXIT_FAILURE);
    }

    // Check if input FITS file have been simulated with TESSIM or XIFUSIM
    char extname[20];
    int extver = 0;
    int tessimOrxifusim = -999;     // 0: tessim, 1: xifusim
    strcpy(extname,"RECORDS");
    fits_movnam_hdu(fptr, ANY_HDU, extname, extver, &status);
    if (status != 0)
    {
        status = 0;
        strcpy(extname,"TESRECORDS");
        fits_movnam_hdu(fptr, ANY_HDU, extname, extver, &status);
        if (status != 0)
        {
            SIXT_ERROR("Cannot move to TESRECORDS HDU in input FITS file");
            return(EXIT_FAILURE);
        }
        else    // TESRECORDS => xifusim
        {
            tessimOrxifusim = 1;
        }
    }
    else // RECORDS => tessim
    {
        tessimOrxifusim = 0;
    }
    if (tessimOrxifusim == -999)
    {
        SIXT_ERROR("Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file");
        return(EXIT_FAILURE);
    }

    // Check if input XML file and XML file to build the library to be used to reconstruct are the same
    if (par.opmode == 1)
    {
        status = checkXmls(&par);
    }

    // Get the sampling rate from the HISTORY keyword from the input FITS file
    // and check with sampling rate from XML file
    // Move to "Primary" HDU
    fits_movabs_hdu(fptr, 1, NULL, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    // Read full Primary HDU and store it in 'headerPrimary'
    char *headerPrimary = NULL;
    int numberkeywords;
    if (fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status))
    {
        free(headerPrimary);
        return(EXIT_FAILURE);
    }
    // Pointer to where the text "sample_rate=" is in HISTORY block
    char *sample_rate_pointer = NULL;
    sample_rate_pointer = strstr (headerPrimary,"sample_rate=");
    // If no 'sample_rate' in HISTORY => sampling frequency from XML is used
    if (sample_rate_pointer)
    {
        // Pointer to the next character to "sample_rate=" (12 characters)
        sample_rate_pointer = sample_rate_pointer + 12;
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
        double samplingrate_HISTORY = atof(characters_after_srate);
        if (((samplingrate_XML != -999.0) && (samplingrate_HISTORY != -999.0)) && (samplingrate_XML != samplingrate_HISTORY))
        {
            SIXT_ERROR("Sampling rate from input FITS file and from XML file do not match");
            return(EXIT_FAILURE);
        }
        *samplingrate = samplingrate_HISTORY;
    }

    // If xifusim file => Get 'trig_reclength' from the HISTORY keyword from the input FITS file
    // 'trig_reclength' is necessary if SIRENA is going to run in THREADING mode
    if (tessimOrxifusim == 1) //xifusim simulated file (with TESRECORDS)
    {
        // Pointer to where the text "trig_reclength=" is in HISTORY block
        char *trig_reclength_pointer = NULL;
        trig_reclength_pointer = strstr (headerPrimary,"trig_reclength=");
        if(!trig_reclength_pointer)
        {
            printf("%s","Attention: 'trig_reclength' is not in the input FITS file (necessary if SIRENA is going to run in THREADING mode)\n");
        }
        else
        {
            // Pointer to the next character to "trig_reclength=" (15 characters)
            trig_reclength_pointer = trig_reclength_pointer + 15;
            char each_character_after_treclength[125];
            snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
            char characters_after_treclength[125];
            snprintf(characters_after_treclength,125,"%c",*trig_reclength_pointer);
            while (*trig_reclength_pointer != ' ')
            {
                trig_reclength_pointer = trig_reclength_pointer + 1;
                snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                strcat(characters_after_treclength,each_character_after_treclength);
            }
            *trigreclength = atoi(characters_after_treclength);
        }
        // Free memory
        if (trig_reclength_pointer != NULL)
        {
            trig_reclength_pointer = NULL;
            free(trig_reclength_pointer);
        }
    }

    // Close file
    fits_close_file(fptr,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Free memory
    if (sample_rate_pointer != NULL)
    {
        sample_rate_pointer = NULL;
        free(sample_rate_pointer);
    }
    if (headerPrimary != NULL)
    {
        free(headerPrimary);
    }

    return(status);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* getSamplingrate_trigreclength: This function gets the sampling rate and the trig_reclength no
*                                matter if only a FITS file or more (inputFile can start with '@'
*                                or not, input file or files can have been simulated with TESSIM
*                                or XIFUSIM)
*
* Parameters:
* - inputFile: Input file name
* - par: Input parameters
* - samplingrate: (In) Sampling rate from XML file => (Out) Sampling rate
* - trigreclength: Necessary if SIRENA is going to run in THREADING mode
* - numfits: Number of FITS files to work with
******************************************************************************/
int getSamplingrate_trigreclength (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength, int* numfits)
{
    // Error status.
    int status=EXIT_SUCCESS;
    char firstcharAux = inputFile[0];
    char firstchar[2] = {firstcharAux , '\0'};

    if (strcmp(firstchar,"@") != 0) // Only A input FITS file
    {
        *numfits = 1;
        status = getSamplingrate_trigreclength_Filei(inputFile, par, samplingrate, trigreclength);
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    }
    else // More than A input FITS file
    {
        // Open @ file to know the number of FITS files to work with
        int numfitsAux = 0;
        FILE *filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");
        if (filetxt == NULL)
        {
            printf("%s","File given in RecordFile does not exist\n");
            status = 104;
        }
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        char filefits[256];
        while(fscanf(filetxt,"%s",filefits)!=EOF)
        {
            numfitsAux++;
        }
        *numfits = numfitsAux;
        fclose(filetxt);

        filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");

        for (int j=0;j<(*numfits);j++)   // For every FITS file
        {
            if (fgets(filefits, 256, filetxt) == NULL)
            {
                fprintf(stderr, "Error reading file line in TXT file with file names\n");
            }
            strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)

            status = getSamplingrate_trigreclength_Filei(inputFile, par, samplingrate, trigreclength);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        }

        fclose(filetxt);
    }

    return(status);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* fillReconstructInitSIRENAGrading: This function reads the grading data from the XML file
*                            and store it in 'reconstruct_init_sirena->grading'
*
*  It also checks if prebuff_0pad input parameter value (preBuffer when 0-padding) is possible
*  depending on the prebuffer values in the XML file
*
* `reconstruct_init_sirena->grading` number of rows = Number of grades in the XML file
* `reconstruct_init_sirena->grading` number of columns = 3 (0->pre, 1->filter length inlcuding prebuffer, 2->prebuffer values)
*
* SIRENA's values (grading=>pre,post and pB) or format XML file (grading=>pre,post and filtlen)
*  post=8192, pB=1000                          pre=494, post=7192, filtlen=8192
*                                                  preBuffer=filtlen-post
*
* Parameters:
* - par: Input parameters
* - det: Pixel detector
* - reconstruct_init_sirena: Parameters to run SIRENA
******************************************************************************/
int fillReconstructInitSIRENAGrading (struct Parameters par, AdvDet *det, ReconstructInitSIRENA** reconstruct_init_sirena)
{
    // Error status.
    int status=EXIT_SUCCESS;

    (*reconstruct_init_sirena)->grading = NULL;
    (*reconstruct_init_sirena)->grading = (Grading*)malloc(sizeof(Grading));
    (*reconstruct_init_sirena)->grading->ngrades = 0;
    (*reconstruct_init_sirena)->grading->value  = NULL;
    (*reconstruct_init_sirena)->grading->gradeData = NULL;

    if ((det->nrecons == 0) && (det->npix == 0))
    {
        SIXT_ERROR("The provided XMLFile does not have the grading info");
        return(EXIT_FAILURE);
    }
    else if ((det->nrecons == 0) && (det->npix != 0))
    {
        if (det->pix->grades == NULL)
        {
            SIXT_ERROR("The provided XMLFile does not have the grading info");
            return(EXIT_FAILURE);
        }
        (*reconstruct_init_sirena)->grading->ngrades=det->pix->ngrades;
        (*reconstruct_init_sirena)->grading->gradeData = gsl_matrix_alloc(det->pix->ngrades,3);
        for (int i=0;i<det->pix->ngrades;i++)
        {
            if ((int) (det->pix->grades[i].gradelim_post) != 8)
            {
                if ((int) (det->pix->grades[i].grade_preBuffer) >= (int) (det->pix->grades[i].gradelim_post))
                {
                    printf("%s %d %s %d %s","preBuffer=",(int) (det->pix->grades[i].grade_preBuffer)," for filter length ",(int) (det->pix->grades[i].gradelim_post),"\n");
                    SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                    return(EXIT_FAILURE);
                }
            }
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,0,(int) (det->pix->grades[i].gradelim_pre));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,1,(int) (det->pix->grades[i].gradelim_post));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,(int) (det->pix->grades[i].grade_preBuffer));
            /*if (((int) (det->pix->grades[i].gradelim_post) == 8) && ((int) (det->pix->grades[i].grade_preBuffer) != 0))
            {
                printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,0);
            }*/
        }
    }
    else if(((det->nrecons != 0) && (det->npix == 0)) || ((det->nrecons == 1) && (det->npix == 1)))
    {
        if (det->recons->grades == NULL)
        {
            SIXT_ERROR("The provided XMLFile does not have the grading info");
            return(EXIT_FAILURE);
        }
        (*reconstruct_init_sirena)->grading->ngrades=det->recons->ngrades;
        (*reconstruct_init_sirena)->grading->gradeData = gsl_matrix_alloc(det->recons->ngrades,3);
        for (int i=0;i<det->recons->ngrades;i++)
        {
            if ((int) (det->recons->grades[i].gradelim_post) != 8)
            {
                if ((int) (det->recons->grades[i].grade_preBuffer) >= (int) (det->recons->grades[i].gradelim_post))
                {
                    printf("%s %d %s %d %s","preBuffer=",(int) (det->recons->grades[i].grade_preBuffer)," for filter length ",(int) (det->recons->grades[i].gradelim_post),"\n");
                    SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                    return(EXIT_FAILURE);
                }
            }
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,0,(int) (det->recons->grades[i].gradelim_pre));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,1,(int) (det->recons->grades[i].gradelim_post));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,(int) (det->recons->grades[i].grade_preBuffer));
            /*if (((int) (det->recons->grades[i].gradelim_post) == 8) && ((int) (det->recons->grades[i].grade_preBuffer) != 0))
            {
                printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,0);
            }*/
        }
    }

    //Check prebuff_0pad input parameter (preBuffer when 0-padding)
    if ((par.opmode == 1) && (strcmp(par.EnergyMethod,"0PAD") == 0))
    {
        gsl_vector *pBsXML = gsl_vector_alloc((*reconstruct_init_sirena)->grading->gradeData->size1);
        gsl_matrix_get_col(pBsXML,(*reconstruct_init_sirena)->grading->gradeData,2);
        if (par.prebuff_0pad > gsl_vector_max(pBsXML))
        {
            SIXT_ERROR("Prebuffer to be used with 0-padding (prebuff_0pad) is bigger than the preBuffers of the optimal filters in the library");
            return(EXIT_FAILURE);
        }
    }

    // Loading in the reconstruct_init structure values related to grading and preBuffer values from the XML file
    gsl_vector *pBi = gsl_vector_alloc((*reconstruct_init_sirena)->grading->ngrades);   // preBuffer values
    gsl_matrix_get_col(pBi,(*reconstruct_init_sirena)->grading->gradeData,2);
    (*reconstruct_init_sirena)->preBuffer_max_value = gsl_vector_max(pBi);
    (*reconstruct_init_sirena)->preBuffer_min_value = gsl_vector_min(pBi);
    gsl_vector *posti = gsl_vector_alloc(1); // Filter length (including preBuffer)
    // post in (grading=>pre,post and pB)
    // filtlen in (grading=>pre,post and filtlen)
    posti = gsl_vector_alloc((*reconstruct_init_sirena)->grading->ngrades);
    gsl_matrix_get_col(posti,(*reconstruct_init_sirena)->grading->gradeData,1);
    (*reconstruct_init_sirena)->post_max_value = gsl_vector_max(posti);
    (*reconstruct_init_sirena)->post_min_value = gsl_vector_min(posti);

    if (pBi != NULL) {gsl_vector_free(pBi); pBi = 0;}
    if (posti != NULL) {gsl_vector_free(posti); posti = 0;}

    return(status);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* callSIRENA_Filei: This function calls SIRENA to build a library or reconstruct energies
*
* Steps:
* - Open record file
* - initializeReconstructionSIRENA
* - Build up TesRecord to read the file
* - Iterate of records and run SIRENA:
*   - reconstructRecordSIRENA
*   - Save events to the event_list
*   - Copy trigger keywords to event file
*   - Close file
*
* Parameters:
* - inputFile: Input file name
* - keywords: Sixt standard keywords structure
* - reconstruct_init_sirena: Parameters to run SIRENA
* - par: Input parameters
* - sampling_rate: Sampling rate
* - trig_reclength: Necessary if SIRENA is going to run in THREADING mode
* - pulsesAll: Structure containing the detected pulses
* - outfile: Output events FITS file
* - ph_id_column_dim: Dimension of PH_ID column
******************************************************************************/
int callSIRENA_Filei(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFileSIRENA* outfile, long ph_id_column_dim)
{
    // Error status.
    int status = EXIT_SUCCESS;

    // Open record file
    TesTriggerFile* record_file;
    record_file = openexistingTesTriggerFile(inputFile,keywords,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr,
        par.LibraryFile, par.TesEventFile, par.flength_0pad, par.prebuff_0pad,
        par.threshold, par.windowSize, par.offset, par.scaleFactor, par.samplesUp, par.samplesDown, par.nSgms,
        par.detectSP, par.opmode, par.detectionMode, par.LrsT,
        par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod,
        par.filtEev, par.Ifit, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter,
        par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.monoenergy,
        par.addCOVAR, par.addINTCOVAR, par.addOFWN, par.intermediate, par.detectFile,
        par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1,
        par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Build up TesRecord to read the file
    TesRecord* record;
    record = newTesRecord(&status);
    if ((record_file->delta_t == -999) && (sampling_rate == -999))
    {
        SIXT_ERROR("Cannot read or get the sampling rate neither from the input FITS file nor the XML file. Please, include the DELTAT keyword (inverse of sampling rate) in the input FITS file before running the task again");
        return(EXIT_FAILURE);
    }
    if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
    if ((*trig_reclength) == -999) *trig_reclength = record_file->trigger_size;
    allocateTesRecord(record,*trig_reclength,record_file->delta_t,0,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Iteration of records and running SIRENA
    // Read numrecords (NAXIS2) to simulate a progress bar
    fitsfile *fptr;
    if (fits_open_file(&fptr, inputFile, READONLY, &status)) {
        SIXT_ERROR("Cannot open the input FITS file to read the number of records (NAXIS2)(in 'callSIRENA_Filei' function).");
        return status;
    }
    fits_movnam_hdu(fptr, ANY_HDU,"TESRECORDS", 0, &status);
    if (status != EXIT_SUCCESS)
    {
        status = EXIT_SUCCESS;
        fits_movnam_hdu(fptr, ANY_HDU,"RECORDS", 0, &status);
        if (status != EXIT_SUCCESS)
        {
            SIXT_ERROR("Cannot find the TESRECORDS (or RECORDS) HDU in the input FITS file to read the number of records (NAXIS2)(in 'callSIRENA_Filei' function).");
            return status;
        }
    }
    int numrecords;
    fits_read_key(fptr, TINT, "NAXIS2", &numrecords, NULL, &status);
    fits_close_file(fptr, &status);

    fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    int lastRecord = 0, nrecord = 0, startRecordGroup = 0;    //last record required for SIRENA library creation
    // Build up TesEventList to recover the results of SIRENA
    TesEventListSIRENA* event_list = newTesEventListSIRENA(&status);
    allocateTesEventListTriggerSIRENA(event_list,par.EventListSize, ph_id_column_dim, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    while(getNextRecordSIRENA(record_file,record,&lastRecord,&startRecordGroup,ph_id_column_dim,&status))
    {
        // Progress bar
        /*float progress;
        if (numrecords == 1)
            progress = 1.0;
        else
            progress = (float)nrecord / (numrecords-1);
        int bar_width = 50;
        int pos = bar_width * progress;
        if (par.opmode == 0)        printf("Building the library |");
        else if (par.opmode == 1)   printf("Reconstructing |");
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
        if (nrecord == numrecords-1)
        {
            printf("\n"); // New line after ending "Building the library..../Reconstructing..."
        }*/

        nrecord = nrecord + 1;
        /*if(nrecord < 5887)
        {
            continue;
        }
        else if(nrecord > 5920)
        {
            status=1;
            CHECK_STATUS_BREAK(status);
        }*/
        /*if(nrecord > 1)
        {
          	status=1;
            CHECK_STATUS_BREAK(status);
        }*/
        if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0)
            || (strcmp(par.EnergyMethod,"I2RDER") == 0))
        {
            strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
        }

        //printf("%s %d %s","** nrecord = ",nrecord,"\n");
        reconstructRecordSIRENA(record,*trig_reclength, event_list,reconstruct_init_sirena,
            lastRecord, startRecordGroup, &pulsesAll, &status);
        CHECK_STATUS_BREAK(status);

        if (lastRecord == 1)
        {
            fits_write_key(outfile->fptr, TINT, "EVENTS", &(pulsesAll->ndetpulses), "Number of detected pulses (including fakes)", &status);
            CHECK_STATUS_BREAK(status);
            fits_write_key(outfile->fptr, TINT, "EVENTSFK", &(pulsesAll->nfakepulses), "Number of invented pulses", &status);
            CHECK_STATUS_BREAK(status);
        }
        if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
        {
            // In THREADING mode, saveEventListToFileSIRENA is not called until finishing with calculus
            // (ordering is necessary previously)
            if(!is_threading()){
                saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,&status);
                CHECK_STATUS_BREAK(status);
                //Reinitialize event list
                event_list->index=0;
            }
            IntegrafreeTesEventListSIRENA(event_list);
        }
    }

    if(is_threading())
    {
        //printf("%s","**Threading...waiting \n");
        th_end(reconstruct_init_sirena, &pulsesAll);
        int i = 1;
        int aux = 1;
        while((aux = th_get_event_list(&event_list, &record)) == 1)
        {
            saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,&status);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
            ++i;
        }
        IntegrafreeTesEventListSIRENA(event_list);
    }

    if (pulsesAll->ndetpulses == 0)
        printf("%s","WARNING: no pulses have been detected\n");
    if (pulsesAll->nfakepulses != 0)
        printf("%s","WARNING: some pulses have been invented\n");

    // Copy trigger keywords to event file
    //copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
    //if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Messages providing info of some columns
    char keywordvalue[9];
    char comment[MAXMSG];

    fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    fits_read_key(outfile->fptr, TSTRING, "TTYPE1", &keywordvalue, NULL, &status);
    strcpy(comment, "Starting time");
    fits_update_key(outfile->fptr, TSTRING, "TTYPE1", keywordvalue, comment, &status);

    fits_read_key(outfile->fptr, TSTRING, "TTYPE2", &keywordvalue, NULL, &status);
    strcpy(comment, "Reconstructed-uncalibrated energy");
    fits_update_key(outfile->fptr, TSTRING, "TTYPE2", keywordvalue, comment, &status);

    fits_read_key(outfile->fptr, TSTRING, "TTYPE3", &keywordvalue, NULL, &status);
    strcpy(comment, "Average first 4 samples (derivative)");
    fits_update_key(outfile->fptr, TSTRING, "TTYPE3", keywordvalue, comment, &status);

    fits_read_key(outfile->fptr, TSTRING, "TTYPE4", &keywordvalue, NULL, &status);
    strcpy(comment, "Optimal filter length");
    fits_update_key(outfile->fptr, TSTRING, "TTYPE4", keywordvalue, comment, &status);

    fits_read_key(outfile->fptr, TSTRING, "TTYPE5", &keywordvalue, NULL, &status);
    strcpy(comment, "Starting time-starting time previous event");
    fits_update_key(outfile->fptr, TSTRING, "TTYPE5", keywordvalue, comment, &status);

    freeTesTriggerFile(&record_file,&status);
    freeTesRecord(&record);

    return (status);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* callSIRENA: This function calls SIRENA to build a library or reconstruct energies no
*             matter if only a FITS file or more (inputFile can start with '@' or not)
*
* Parameters:
* - inputFile: Input file name
* - keywords: Sixt standard keywords structure
* - reconstruct_init_sirena: Parameters to run SIRENA
* - par: Input parameters
* - sampling_rate: Sampling rate
* - trig_reclength: Necessary if SIRENA is going to run in THREADING mode
* - pulsesAll: Structure containing the detected pulses
* - outfile: Output events FITS file
* - ph_id_column_dim: Dimension of PH_ID column
******************************************************************************/
int callSIRENA(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFileSIRENA* outfile, long ph_id_column_dim)
{
    // Error status.
    int status=EXIT_SUCCESS;

    char firstcharAux = inputFile[0];
    char firstchar[2] = {firstcharAux , '\0'};

    if (strcmp(firstchar,"@") != 0) // Only A input FITS file
    {
        status = callSIRENA_Filei(inputFile, keywords, reconstruct_init_sirena, par, sampling_rate, trig_reclength, pulsesAll, outfile, ph_id_column_dim);
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    }
    else // More than A input FITS file
    {
        // Open @ file to know the number of FITS files to work with
        int numfits = 0;
        FILE *filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");
        if (filetxt == NULL)
        {
            printf("%s","File given in RecordFile does not exist\n");
            status = 104;
        }
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        char filefits[256];
        while(fscanf(filetxt,"%s",filefits)!=EOF)
        {
            numfits++;
        }
        fclose(filetxt);

        filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");

        for (int j=0;j<numfits;j++)   // For every FITS file
        {
            if (fgets(filefits, 256, filetxt) == NULL)
            {
                fprintf(stderr, "Error reading file line in TXT file with file names\n");
            }
            strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)

            status = callSIRENA_Filei(inputFile, keywords, reconstruct_init_sirena, par, sampling_rate, trig_reclength, pulsesAll, outfile, ph_id_column_dim);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        }

        fclose(filetxt);
    }

    return (status);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* checkpreBuffer function: Check if the preBuffer value (HISTORY) in the library is the same to preBuffer input parameter
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int checkpreBuffer(struct Parameters* const par)
{
    // Error status.
    int status = EXIT_SUCCESS;

    fitsfile* libptr = NULL;

    // If the library exists => Open it
    fits_open_file(&libptr, par->LibraryFile, READONLY, &status);
    if (status == EXIT_SUCCESS)
    {
        // Move to "Primary" HDU of the library file
        if (fits_movabs_hdu(libptr, 1, NULL, &status))
        {
            SIXT_ERROR("Error moving to Primary HDU in the LibraryFile");
            return(EXIT_FAILURE);
        }

        char *libheaderPrimary = NULL;
        int numberkeywords;
        if (fits_hdr2str(libptr, 0, NULL, 0,&libheaderPrimary, &numberkeywords, &status))
        {
            free(libheaderPrimary);
            SIXT_ERROR("Error reading Primary HDU in the LibraryFile");
            return(EXIT_FAILURE);
        }

        fits_close_file(libptr,&status);
    }
    else status = EXIT_SUCCESS;

    return(status);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


TesEventListSIRENA* newTesEventListSIRENA(int* const status){
	TesEventListSIRENA* event_list= malloc(sizeof*event_list);
	if (NULL==event_list) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventList failed");
		return(event_list);
	}
	// Initialize pointers with NULL
	event_list->event_indexes=NULL;
	event_list->pulse_heights=NULL;
	event_list->avgs_4samplesDerivative=NULL;
	event_list->Es_lowres=NULL;
	event_list->grading=NULL;
	event_list->phis=NULL;
	event_list->lagsShifts=NULL;
	event_list->bsln=NULL;
	event_list->rmsbsln=NULL;
	event_list->energies=NULL;
	event_list->grades1=NULL;
	event_list->grades2=NULL;
	event_list->ph_ids_array = NULL;
    event_list->ph_ids_array_size1 = 0;
    event_list->ph_ids_array_size2 = 0;
    event_list->pix_ids=NULL;
	event_list->risetimes=NULL;
	event_list->falltimes=NULL;

	// Initialize values
	event_list->size=0;
	event_list->size_energy=0;
	event_list->index=0;

	return(event_list);
}


void freeTesEventListSIRENA(TesEventListSIRENA* event_list)
{
    if (event_list == NULL) return; // Avoid dereferencing a NULL pointer

    // Free dynamically allocated arrays
    if (event_list->event_indexes) { free(event_list->event_indexes); event_list->event_indexes = NULL; }
    if (event_list->grades1) { free(event_list->grades1); event_list->grades1 = NULL; }
    if (event_list->pulse_heights) { free(event_list->pulse_heights); event_list->pulse_heights = NULL; }
    if (event_list->energies) { free(event_list->energies); event_list->energies = NULL; }
    if (event_list->avgs_4samplesDerivative) { free(event_list->avgs_4samplesDerivative); event_list->avgs_4samplesDerivative = NULL; }
    if (event_list->Es_lowres) { free(event_list->Es_lowres); event_list->Es_lowres = NULL; }
    if (event_list->phis) { free(event_list->phis); event_list->phis = NULL; }
    if (event_list->lagsShifts) { free(event_list->lagsShifts); event_list->lagsShifts = NULL; }
    if (event_list->bsln) { free(event_list->bsln); event_list->bsln = NULL; }
    if (event_list->rmsbsln) { free(event_list->rmsbsln); event_list->rmsbsln = NULL; }
    if (event_list->grading) { free(event_list->grading); event_list->grading = NULL; }
    if (event_list->grades2) { free(event_list->grades2); event_list->grades2 = NULL; }
    if (event_list->pix_ids) { free(event_list->pix_ids); event_list->pix_ids = NULL; }
    if (event_list->tstarts) { free(event_list->tstarts); event_list->tstarts = NULL; }
    if (event_list->tends) { free(event_list->tends); event_list->tends = NULL; }
    if (event_list->risetimes) { free(event_list->risetimes); event_list->risetimes = NULL; }
    if (event_list->falltimes) { free(event_list->falltimes); event_list->falltimes = NULL; }

    // Free the 2D array ph_ids_array
    if (event_list->ph_ids_array) {
        for (int i = 0; i < event_list->ph_ids_array_size2; i++) {
            if (event_list->ph_ids_array[i]) {
                free(event_list->ph_ids_array[i]);
                event_list->ph_ids_array[i] = NULL;
            }
        }
        free(event_list->ph_ids_array);
        event_list->ph_ids_array = NULL;
    }

    // Free the struct itself
    free(event_list);
    event_list = NULL;
}

TesEventFileSIRENA* newTesEventFileSIRENA(int* const status){
	TesEventFileSIRENA* file=malloc(sizeof(*file));
	if (NULL==file) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventFile failed");
		return(file);
	}

	// Initialize pointers with NULL.
	file->fptr    =NULL;

	// Initialize values.
	file->row  	    =1;
	file->nrows     =0;
	file->timeCol   =1;
	file->energyCol =2;
	file->avg_4samplesDerivativeCol =3;
	file->E_lowresCol =4;
	file->grade1Col =5;
	file->grade2Col =6;
	file->phiCol =7;
	file->lagsShiftCol =8;
	file->bslnCol =9;
	file->rmsbslnCol =10;
	file->pixIDCol  =11;
	file->phIDCol   =12;
	file->riseCol   =13;
	file->fallCol   =14;
	file->raCol     =15;
	file->decCol    =16;
	file->detxCol   =17;
	file->detyCol   =18;
	file->gradingCol=19;
	file->srcIDCol  =20;
	file->nxtCol    =21;
	file->extCol    =22;

	return(file);
}

/** TesEventFileSIRENA Destructor. */
void freeTesEventFileSIRENA(TesEventFileSIRENA* file, int* const status){
	if (NULL!=file) {
		if (NULL!=file->fptr) {
			fits_close_file(file->fptr, status);
			CHECK_STATUS_VOID(*status);
			headas_chat(5, "closed TesEventFile list file\n");
		}
		free(file);
		file=NULL;
	}
}

TesEventFileSIRENA* opennewTesEventFileSIRENA(const char* const filename,
				  SixtStdKeywords* keywords,
				  const char* const sirenaVersion,
                  int ph_id_size,
				  const char clobber,
				  int* const status){
	TesEventFileSIRENA* file = newTesEventFileSIRENA(status);
	CHECK_STATUS_RET(*status, file);

	int exists;
	char buffer[MAXFILENAME];
	sprintf(buffer,"%s",filename);
	fits_file_exists(buffer, &exists, status);
	CHECK_STATUS_RET(*status,file);
	if (0!=exists) {
		if (0!=clobber) {
			// Delete the file.
			remove(buffer);
		} else {
			// Throw an error.
			char msg[MAXMSG];
            snprintf(msg, sizeof(msg), "file '%.*s' already exists",
                (int)(sizeof(msg) - strlen("file '' already exists") - 1), buffer);
			SIXT_ERROR(msg);
			*status=EXIT_FAILURE;
			CHECK_STATUS_RET(*status,file);
		}
	}
	fits_create_file(&file->fptr,buffer, status);
	CHECK_STATUS_RET(*status,file);
	int logic=(int)'T';
	int bitpix=8;
	int naxis=0;
	fits_update_key(file->fptr, TLOGICAL, "SIMPLE", &(logic), NULL, status);
	fits_update_key(file->fptr, TINT, "BITPIX", &(bitpix), NULL, status);
	fits_update_key(file->fptr, TINT, "NAXIS", &(naxis), NULL, status);
	sixt_add_fits_stdkeywords(file->fptr,1,keywords,status);
	CHECK_STATUS_RET(*status,file);

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    const char * chardate = asctime (timeinfo);
    char keyvalstr[1000];
    strcpy(keyvalstr,chardate);
    fits_write_key(file->fptr,TSTRING,"CREADATE",keyvalstr,NULL,status);
    CHECK_STATUS_RET(*status,file);

	fits_write_key(file->fptr,TSTRING,"SIRENAV",(char *)sirenaVersion,NULL,status);
	CHECK_STATUS_RET(*status,file);

	// Create table
    char temp[10];
    sprintf(temp, "%dJ", ph_id_size);
	char   *ttype[]={"TIME","SIGNAL","AVG4SD","ELOWRES","GRADE1","GRADE2","PHI","LAGS","BSLN","RMSBSLN","PIXID","PH_ID","RISETIME","FALLTIME","RA","DEC","DETX","DETY","GRADING","SRC_ID","N_XT","E_XT"};
	char *tform[]={"1D", "1D", "1D", "1D", "1K", "1K", "1D", "1J", "1D", "1D", "1J", temp, "1D", "1D", "1D", "1D", "1E", "1E", "1I", "1J", "1I", "1D"};
	char *tunit[]={"s", "keV", "", "keV", "", "", "", "", "", "", "", "", "s", "s", "deg", "deg", "m", "m", "", "", "", "keV"};

	fits_create_tbl(file->fptr, BINARY_TBL, 0, 22, ttype, tform, tunit,"EVENTS", status);
	sixt_add_fits_stdkeywords(file->fptr,2,keywords,status);
	CHECK_STATUS_RET(*status,file);

	int firstpix=0,lastpix=0,numberpix=0;
	float monoen=-1.;
	long nes_tot=0,net_tot=0;
	fits_update_key(file->fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in record file", status);
	fits_update_key(file->fptr, TINT, "LASTPIX", &lastpix, "Last pixel in record file", status);
	fits_update_key(file->fptr, TINT, "NPIX", &numberpix, "Number of pixels in record file", status);
	fits_update_key(file->fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);
	fits_update_key(file->fptr, TLONG, "NESTOT", &nes_tot, "Total number of events simulated", status);
	fits_update_key(file->fptr, TLONG, "NETTOT", &net_tot, "Total number of events actually triggered", status);
	CHECK_STATUS_RET(*status,file);

	return(file);

}


void saveEventListToFileSIRENA(TesEventFileSIRENA* file,TesEventListSIRENA * event_list,
double start_time,double delta_t,int* const status){
	//Save time, PIXID and dummy grading column
	double time;
	int dummy_grading = 0;
	for (int i = 0 ; i<event_list->index ; i++){
		time = start_time + delta_t*event_list->event_indexes[i];
		fits_write_col(file->fptr, TDOUBLE, file->timeCol,
					   file->row, 1, 1, &time, status);
		fits_write_col(file->fptr, TINT, file->pixIDCol,
					file->row, 1, 1, &event_list->pix_ids[i], status);
		fits_write_col(file->fptr, TINT, file->gradingCol,
				file->row, 1, 1, &dummy_grading, status);
		CHECK_STATUS_VOID(*status);
		file->row++;
	}
	file->row = file->row - event_list->index;
	CHECK_STATUS_VOID(*status);

	//Save energy column
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
					file->row, 1, event_list->index, event_list->energies, status);
	CHECK_STATUS_VOID(*status);

	//Save avgs_4samplesDerivative (AVG4SD) column
 	fits_write_col(file->fptr, TDOUBLE, file->avg_4samplesDerivativeCol,
					file->row, 1, event_list->index, event_list->avgs_4samplesDerivative, status);
	CHECK_STATUS_VOID(*status);

    //Save Es_lowres (ELOWRES) column
 	fits_write_col(file->fptr, TDOUBLE, file->E_lowresCol,
					file->row, 1, event_list->index, event_list->Es_lowres, status);
	CHECK_STATUS_VOID(*status);

	//Save phis (PHI) column
 	fits_write_col(file->fptr, TDOUBLE, file->phiCol,
					file->row, 1, event_list->index, event_list->phis, status);
	CHECK_STATUS_VOID(*status);

	//Save lagsShifts (LAGS) column
 	fits_write_col(file->fptr, TINT, file->lagsShiftCol,
					file->row, 1, event_list->index, event_list->lagsShifts, status);
	CHECK_STATUS_VOID(*status);

	//Save bsln (BSLN) column
 	fits_write_col(file->fptr, TDOUBLE, file->bslnCol,
					file->row, 1, event_list->index, event_list->bsln, status);
	CHECK_STATUS_VOID(*status);

	//Save rmsbsln (RMSBSLN) column
 	fits_write_col(file->fptr, TDOUBLE, file->rmsbslnCol,
					file->row, 1, event_list->index, event_list->rmsbsln, status);
	CHECK_STATUS_VOID(*status);

	//Save grading (GRADING) column
 	fits_write_col(file->fptr, TINT, file->gradingCol,
					file->row, 1, event_list->index, event_list->grading, status);
	CHECK_STATUS_VOID(*status);

	//Save grade1 column
	fits_write_col(file->fptr, TLONGLONG, file->grade1Col,
					file->row, 1, event_list->index, event_list->grades1, status);
	CHECK_STATUS_VOID(*status);

	//Save grade2 column
	fits_write_col(file->fptr, TLONGLONG, file->grade2Col,
					file->row, 1, event_list->index, event_list->grades2, status);
	CHECK_STATUS_VOID(*status);

    //Save RISETIME column
 	fits_write_col(file->fptr, TDOUBLE, file->riseCol,
					file->row, 1, event_list->index, event_list->risetimes, status);
	CHECK_STATUS_VOID(*status);

    //Save FALLTIME column
 	fits_write_col(file->fptr, TDOUBLE, file->fallCol,
					file->row, 1, event_list->index, event_list->falltimes, status);
	CHECK_STATUS_VOID(*status);

    //Save PH_ID column

    if (event_list->ph_ids_array_size1 != 0)
    {
        // Bidimensional matrix => Unidimensional matrix
        long *flat_ph_ids_array = malloc(event_list->ph_ids_array_size1 * event_list->ph_ids_array_size2 * sizeof(long));
        int idx = 0;

        for (int kkk = 0; kkk < event_list->ph_ids_array_size1; kkk++) {
            for (int kkk1 = 0; kkk1 < event_list->ph_ids_array_size2; kkk1++) {
                flat_ph_ids_array[idx] = event_list->ph_ids_array[kkk][kkk1];
                idx++;
            }
        }

        fits_write_col(file->fptr, TLONG, file->phIDCol,
                       file->row, 1, event_list->ph_ids_array_size1 * event_list->ph_ids_array_size2, flat_ph_ids_array, status);
        free(flat_ph_ids_array);
        flat_ph_ids_array = NULL;
    }

	file->row = file->row + event_list->index;
	file->nrows+= event_list->index;
}

/** Allocates memory for a TesEventList structure for the triggering stage:
 *  only event_index, pulse_height */
void allocateTesEventListTriggerSIRENA(TesEventListSIRENA* event_list,int size, int size_phid, int* const status){
	event_list->event_indexes = malloc(size*sizeof*(event_list->event_indexes));
	if (NULL == event_list->event_indexes){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for event_index array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->pulse_heights = malloc(size*sizeof*(event_list->pulse_heights));
	if (NULL == event_list->pulse_heights){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for pulse_height array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->grades1 = malloc(size*sizeof*(event_list->grades1));
	if (NULL == event_list->grades1){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for grade1 array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->ph_ids_array = malloc(size*sizeof*(event_list->ph_ids_array));
	if (NULL == event_list->ph_ids_array){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for ph_ids_array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}
	// Rows in the matrix
    for (int i = 0; i < size_phid; i++) {
        event_list->ph_ids_array[i] = NULL;
    }
    event_list->ph_ids_array_size1 = 0;
    event_list->ph_ids_array_size2 = size_phid;

	event_list->size = size;
}


/** Populates a TesRecord structure with the next record */
int getNextRecordSIRENA(TesTriggerFile* const file,TesRecord* record,int *lastRecord,int *startRecordGroup,long ph_id_column_dim,int* const status){
  int anynul=0;
  char tform2ADC[20];
  LONGLONG rec_trigsize;

  (*lastRecord) = 0;
  (*startRecordGroup) = 0;

  if (NULL==file || NULL==file->fptr) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("No opened trigger file to read from");
    CHECK_STATUS_RET(*status,0);
  }

  if (file->row<=file->nrows) {
    // get length of this record
    // (although unlikely, we might have a very large file, so we best
    // use the LONGLONG interface to the descriptor

    // read TFORM for ADC to know if it is FIXED or VARIABLE length
    fits_read_key(file->fptr,TSTRING, "TFORM2", &tform2ADC, NULL, status);
    if(strstr(tform2ADC, "(") != NULL){
      LONGLONG offset;
      fits_read_descriptll(file->fptr,file->trigCol,file->row,&rec_trigsize,&offset,status);
      CHECK_STATUS_RET(*status,0);

    }else{
      LONGLONG col_width;
      int adc_col_typecode;
      fits_get_coltypell(file->fptr,file->trigCol,&adc_col_typecode,
			 &rec_trigsize,&col_width,status);
      CHECK_STATUS_RET(*status,0);
    }

    fits_read_col(file->fptr, TLONG, file->pixIDCol,
                  file->row,1,1,0,&(record->pixid), &anynul,status);
    CHECK_STATUS_RET(*status,0);
    fits_read_col(file->fptr, TDOUBLE, file->timeCol,
                  file->row,1,1,0,&(record->time), &anynul,status);
    CHECK_STATUS_RET(*status,0);

    fits_read_col(file->fptr, TLONG, file->extendCol,
		  file->row,1,1,0,&(record->extend), &anynul,status);
    if (*status != 0) // No EXTEND column
    {
        *status = 0;
        //// when ADC is integer
        ////fits_read_col(file->fptr, TUSHORT, file->trigCol,
        ////		  file->row,1,record->trigger_size,0,record->adc_array, &anynul,status);
        //// when ADC is DOUBLE
        fits_read_col(file->fptr, TDOUBLE, file->trigCol,
                      file->row,1,record->trigger_size,0,record->adc_double, &anynul,status);
        CHECK_STATUS_RET(*status,0);

        fits_read_col(file->fptr, TLONG, file->ph_idCol,
                      file->row,1,1,0,&(record->phid_list->phid_array[0]), &anynul,status);
        if (*status != 0)	// Simulated files have the PH_ID column but empty
        {
            record->phid_list->phid_array[0] = 0;
            *status = 0;
        }

        (*startRecordGroup) = file->row;

        file->row++;
    }
    else // EXTEND column
    {
        // Starting a new record => record->trigger_size must begin in the ADC column size
        if (record->trigger_size!=(unsigned long) rec_trigsize) {
            resizeTesRecord(record,(unsigned long) rec_trigsize,status);
            CHECK_STATUS_RET(*status,0);
        }

        // It is not necessary to resize the record
        int nrowsTogether = 0;
        if (record->extend == 0)
        {
            // when ADC is integer
            //fits_read_col(file->fptr, TUSHORT, file->trigCol,
            //	 	  file->row,1,record->trigger_size,0,record->adc_array, &anynul,status);
            // when ADC is DOUBLE
            fits_read_col(file->fptr, TDOUBLE, file->trigCol,
                          file->row,1,record->trigger_size,0,record->adc_double, &anynul,status);
            CHECK_STATUS_RET(*status,0);
            (*startRecordGroup) = file->row;
        }
        else // It is necessary to resize the record
        {
            int sizeADC = record->trigger_size;
            int sizeToWrite;

            int sizeLastRow = 0;
            div_t nrowsTogether_div_t = div(record->trigger_size+record->extend,record->trigger_size);

            sizeLastRow = nrowsTogether_div_t.rem;
            if (sizeLastRow == 0)
            {
                nrowsTogether = nrowsTogether_div_t.quot;
            }
            else
            {
                nrowsTogether = nrowsTogether_div_t.quot+1;
            }
            //printf("%s %d %s","nrowsTogether: ",nrowsTogether,"\n");
            resizeTesRecord(record,(unsigned long) record->extend + record->trigger_size,status);
            double *singleRecord = NULL;
            int rowToRead = file->row;
            int index = 0;
            (*startRecordGroup) = file->row;
            for (int i=0;i<nrowsTogether;i++)
            {
                singleRecord = malloc(sizeADC*sizeof(singleRecord));
                // when ADC is integer
                //fits_read_col(file->fptr, TUSHORT, file->trigCol,
                //	 	  file->row,1,record->trigger_size,0,record->adc_array, &anynul,status);
                // when ADC is DOUBLE
                fits_read_col(file->fptr, TDOUBLE, file->trigCol,
                              rowToRead,1,sizeADC,0,singleRecord, &anynul,status);
                CHECK_STATUS_RET(*status,0);

                if (i != nrowsTogether-1)
                {
                    sizeToWrite = sizeADC;
                }
                else
                {
                    //if (sizeLastRow==0) sizeToWrite = sizeADC;
                    //else sizeToWrite = sizeLastRow;
                    sizeToWrite = sizeLastRow;
                }
                for (int j=0;j<sizeToWrite;j++)
                {
                    record->adc_double[j+index] = singleRecord[j];
                }
                free(singleRecord);
                index = index + sizeADC;
                rowToRead++;
            }
        }

        free(record->phid_list->phid_array);
        record->phid_list->phid_array = malloc(ph_id_column_dim*sizeof(*(record->phid_list->phid_array)));
        record->phid_list->size = ph_id_column_dim;
        fits_read_col(file->fptr, TLONG, file->ph_idCol,
                      file->row,1,ph_id_column_dim,0,record->phid_list->phid_array, &anynul,status);
        CHECK_STATUS_RET(*status,0);

        if (record->extend == 0)
        {
            file->row++;
        }
        else
        {
            file->row = file->row + nrowsTogether;
        }
    }

    if (file->row>file->nrows)  (*lastRecord) = 1;

    return(1);
  } else {
    return(0);
  }
}
