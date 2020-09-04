.. SIRENA functions description

.. _FUNCTIONS:

.. role:: pageblue
.. role:: red

.. To refer to one of this functions use :cpp:func:`functionName` and to one of its members :cpp:member:`functionName::memberName`


#######################
SIRENA functions
#######################

.. ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 
   :ref:`A <A>` :ref:`B <B>` :ref:`C <C>` :ref:`D <D>` :ref:`E <E>` :ref:`F <F>` :ref:`G <G>` :ref:`H <H>` :ref:`I <I>` :ref:`J <J>` :ref:`K <K>` :ref:`L <L>` :ref:`M <M>`
   ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 

   ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 
   :ref:`N <N>` :ref:`O <O>` :ref:`P <P>` :ref:`Q <Q>` :ref:`R <R>` :ref:`S <S>` :ref:`T <T>` :ref:`U <U>` :ref:`V <V>` :ref:`W <W>` :ref:`X <X>` :ref:`Y <Y>` :ref:`Z <Z>` 
   ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 


============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============
:ref:`A <A>` :ref:`B <B>` :ref:`C <C>` :ref:`D <D>` :ref:`E <E>` :ref:`F <F>` :ref:`G <G>` :ref:`H <H>` :ref:`I <I>` :ref:`J <J>` :ref:`K <K>` :ref:`L <L>` :ref:`M <M>` :ref:`N <N>`
============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 

============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 
:ref:`O <O>` :ref:`P <P>` :ref:`Q <Q>` :ref:`R <R>` :ref:`S <S>` :ref:`T <T>` :ref:`U <U>` :ref:`V <V>` :ref:`W <W>` :ref:`X <X>` :ref:`Y <Y>` :ref:`Z <Z>` 
============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ ============ 

Search functions by name at :ref:`genindex`.


.. _A: 
    
.. cpp:function:: int addFirstRow(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int runF0orB0val, gsl_vector *E, gsl_vector *PHEIGHT, gsl_matrix *PULSE, gsl_matrix *PULSEB0, gsl_matrix *MF, gsl_matrix *MFB0, gsl_matrix *COVAR, gsl_matrix *WEIGHT, gsl_matrix *PULSEMaxLengthFixedFilter)

    Located in file: *tasksSIRENA.cpp*
    
    This function writes the first row of the library (without intermediate AB-related values, because it would be necessary to have at least two rows=energies in the library). It also writes the *FIXFILTT* and *FIXFILTF* HDUs with the optimal filters in the time and frequency domain with fixed legnths (base-2 values) and the *PRCLOFWM* HDU with the precalculated values for optimal filtering and :option:`EnergyMethod` = **WEIGHTM**.
    
    - Declare variables
    - Write in the first row of the library FITS file some columns with the info provided by the input GSL vectors :cpp:member:`E`, :cpp:member:`PHEIGHT`, :cpp:member:`PULSE`, :cpp:member:`PULSEB0`,             :cpp:member:`MF` and :cpp:member:`MFB0` (and :cpp:member:`COVAR` and :cpp:member:`WEIGHT` if :option:`hduPRCLOFWM` =1) (and :cpp:member:`PULSEMaxLengthFixedFilter` if :option:`largeFilter` > :option:`PulseLength`)
    - Writing HDUs with fixed filters in time (*FIXFILTT*) and frequency (*FIXFILTF*), **Tx** and **Fx** columns respectively (calculating the optimal filters, :cpp:func:`calculus_optimalFilter`).
      In time domain **Tx** columns are real numbers but in frequency domain **Fx** columns are complex numbers (so real parts are written in the first half of the column and imaginary parts in the second one)
    - Calculate and write the pre-calculated values by using the noise weight matrix from noise intervals (M'WM)^{-1}M'W for different lengths, **OFWx** columns in *PRCLOFWM*
   
    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
     
    fitsfile** **inLibObject**

        FITS object containing information of the library FITS file  
        
    double **samprate**

        Sampling rate
        
    int **runF0orB0val**
    
        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    gsl_vector* **E**
    
        First energy to be included in the library 
        
    gsl_vector* **PHEIGHT**
    
        Pulse height associated to the first energy to be included in the library
        
    gsl_matrix* **PULSE**

        Pulse template associated to the first energy to be included in the library
        
    gsl_matrix* **PULSEB0** 
    
        Pulse template without baseline associated to the first energy to be included in the library
    
    gsl_matrix* **MF**
    
        Matched filter associated to the first energy to be included in the library
        
    gsl_matrix* **MFB0**
    
        Matched filter (baseline subtracted) associated to the first energy to be included in the library
        
    gsl_matrix* **COVAR**
    
        Covariance matrix associated to the first energy to be included in the library
    
    gsl_matrix* **WEIGHT**
    
        Weight matrix associated to the first energy to be included in the library

    gsl_matrix* **PULSEMaxLengthFixedFilter**

        Pulse template whose length is :option:`largeFilter` associated to the first energy to be included in the library

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
     
    .. cpp:member:: fitsfile** inLibObject

        FITS object containing information of the library FITS file  
        
    .. cpp:member:: double samprate

        Sampling rate
        
    .. cpp:member:: int runF0orB0val
    
        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    .. cpp:member:: gsl_vector* E
    
        First energy to be included in the library 
        
    .. cpp:member:: gsl_vector* PHEIGHT
    
        Pulse height associated to the first energy to be included in the library
        
    .. cpp:member:: gsl_matrix* PULSE

        Pulse template associated to the first energy to be included in the library
        
    .. cpp:member:: gsl_matrix* PULSEB0 
    
        Pulse template without baseline associated to the first energy to be included in the library
    
    .. cpp:member:: gsl_matrix* MF
    
        Matched filter associated to the first energy to be included in the library
        
    .. cpp:member:: gsl_matrix* MFB0
    
        Matched filter (baseline subtracted) associated to the first energy to be included in the library
        
    .. cpp:member:: gsl_matrix* COVAR
    
        Covariance matrix associated to the first energy to be included in the library
    
    .. cpp:member:: gsl_matrix* WEIGHT
    
        Weight matrix associated to the first energy to be included in the library
	
    .. cpp:member:: gsl_matrix* PULSEMaxLengthFixedFilter

        Pulse template whose length is :option:`largeFilter` associated to the first energy to be included in the library

        
.. cpp:function:: int align(double samprate, gsl_vector **vector1, gsl_vector **vector2)
        
    Located in file: *tasksSIRENA.cpp*
    
    Based on :cite:`GilPita2005`
    
    This function aligns :cpp:member:`vector1` with :cpp:member:`vector2` (by delaying or moving forward :cpp:member:`vector2`) assuming that :cpp:member:`vector1` and :cpp:member:`vector2` are shifted replicas of the same function.

    From the discrete function :math:`x[n] (n=0,...,N-1,N)` and according to the time shifting property of the Fourier transform:

    .. math::

        & x[n]    <------> X[f]\\
        & x[n-m]  <------> X[f] exp(-j2\cdot\pi\cdot m/N)

    If :math:`\mathit{Shift} = m` then :math:`\mathit{PhaseDueToTheShift}= 2\pi m/N` and thus, :math:`m = \mathit{PhaseDueToTheShift}\cdot N/(2\pi)`

    1) Declare variables
    
    2) FFT of :cpp:member:`vector1`
    
    3) FFT of :cpp:member:`vector2`
    
    4) (Phases of the *FFT_vector1* and *FFT_vector2*) :math:`*size/(2\pi)`

    5) Shift between the input vectors
    
    6) *shiftdouble* into *shiftint* (because we are working with samples)

    7) Move forward or delay :cpp:member:`vector1` depending on positive or negative shift

    **Members/Variables**

    double **samprate**

        Sampling rate

    gsl_vector** **vector1**

        GSL vector with input vector

    gsl_vector** **vector2** 

        GSL with input vector which is delayed or moved forward to be aligned with :cpp:member:`vector1`
        
    .. cpp:member:: double samprate

        Sampling rate

    .. cpp:member:: gsl_vector** vector1

        GSL vector with input vector

    .. cpp:member:: gsl_vector** vector2 

        GSL with input vector which is delayed or moved forward to be aligned with :cpp:member:`vector1`
        
            
.. _B:

.. _C:

.. cpp:function:: int calculateEnergy(gsl_vector *vector, int pulseGrade, gsl_vector *filter, gsl_vector_complex *filterFFT,int runEMethod, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, int domain, double samprate, gsl_vector *Pab, gsl_matrix *PRCLWN, gsl_matrix *PRCLOFWM, double *calculatedEnergy, int numlags, double *tstartNewDev, int productSize, int tooshortPulse_NoLags)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the energy of a pulse (:cpp:member:`vector`) depending on the :option:`EnergyMethod`, :option:`OFNoise`, and the :option:`FilterDomain` selected from input parameters.

    a) **OPTFILT** and **NSD** (= **I2R** or **I2RALL**, **I2RNOL** or **I2RFITTED**): Optimal filter = Wiener filter  (see :ref:`optimalFilter_NSD`)

    Once the filter template has been created (:cpp:member:`filter` or :cpp:member:`filterFFT`), pulse height analysis is performed by aligning the template with a pulse and multiplying each point in the template by the corresponding point in the pulse. The sum of these products is the energy.

    In the practice, the alignment of the pulse relative to the trigger is not completely accurate, so a number of *n* lags could be used in order to find the peak value of the energy. The *n* peak values are fitted to a parabola to find the most accurate energy (:option:`LagsOrNot`) and a corrected starting time.
    
    a) **OPTFILT** and **WEIGHTM** (= **I2R** or **I2RALL**, **I2RNOL** or **I2RFITTED**) (see :ref:`optimalFilter_WEIGHTM`)

    c) **WEIGHT** and **WEIGHTN** (see :ref:`covMatrices`)


    **Members/Variables**
    
    gsl_vector* **vector**
    
        Pulse whose energy has to be determined
        
    int **pulseGrade** 
    
        Grade of the input pulse (to decide whether a full or only a rough estimation of energy is required). 
        
    gsl_vector* **filter**
    
        Optimal filter in time domain
    
    gsl_vector_complex* **filterFFT**
    
        Optimal filter in frequency domain
        
    int **runEMethod**
    
        - :option:`EnergyMethod` = **OPTFILT** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2R** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RALL** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RNOL** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RFITTED** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **WEIGHT** :math:`\Rightarrow` :cpp:member:`runEMethod` = 1
        - :option:`EnergyMethod` = **WEIGHTN** :math:`\Rightarrow` :cpp:member:`runEMethod` = 2

    int **indexEalpha**
    
        Index of the energy lower than the energy of the pulse which is being analyzed
        
    int **indexEbeta** 
    
        Index of the energy higher than the energy of the pulse which is being analyzed
        
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    int **domain**
    
        - :option:`FilterDomain` = **T** :math:`\Rightarrow` :cpp:member:`domain` = 0
        - :option:`FilterDomain` = **F** :math:`\Rightarrow` :cpp:member:`domain` = 1

    double **samprate**
    
        Sampling rate in Hz
        
    gsl_vector* **Pab** 
    
        **PAB** column in the library 
        
    gsl_vector* **PRCLWN** 
    
        **PCLx** column in the library

    gsl_vector* **PRCLOFWM** 
    
        **OFWx** column in the library
        
    double* **calculatedEnergy**
    
        Calculated energy in eV.
        
    int **numlags**
    
        Number of lags (if option:`EnergyMethod` = **OPTFILT** or **I2R** or **I2RALL** or **I2RNOL** or **I2RFITTED** and :option:`OFNoise` = **NSD**)
        
    double **tstartNewDev**
    
        Addional deviation of the starting time (if :option:`LagsOrNot` = 1)    
        
    int **productSize**
    
        Size of the scalar product to be calculated
        
    int **tooshortPulse_NoLags**
    
        Pulse too short to apply lags (1) or not (0)
    
    .. cpp:member:: gsl_vector* vector
    
        Pulse whose energy has to be determined
        
    .. cpp:member:: int pulseGrade 
    
        Grade of the input pulse (to decide whether a full or only a rough estimation of energy is required). 
        
    .. cpp:member:: gsl_vector* filter
    
        Optimal filter in time domain
    
    .. cpp:member:: gsl_vector_complex* filterFFT
    
        Optimal filter in frequency domain
        
    .. cpp:member:: int runEMethod 
    
        - :option:`EnergyMethod` = **OPTFILT** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2R** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RALL** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RNOL** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **I2RFITTED** :math:`\Rightarrow` :cpp:member:`runEMethod` = 0
        - :option:`EnergyMethod` = **WEIGHT** :math:`\Rightarrow` :cpp:member:`runEMethod` = 1
        - :option:`EnergyMethod` = **WEIGHTN** :math:`\Rightarrow` :cpp:member:`runEMethod` = 2

    .. cpp:member:: int indexEalpha
    
        Index of the energy lower than the energy of the pulse which is being analyzed
        
    .. cpp:member:: int indexEbeta 
    
        Index of the energy higher than the energy of the pulse which is being analyzed
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    .. cpp:member:: int domain
    
        - :option:`FilterDomain` = **T** :math:`\Rightarrow` :cpp:member:`domain` = 0
        - :option:`FilterDomain` = **F** :math:`\Rightarrow` :cpp:member:`domain` = 1

    .. cpp:member:: double samprate
    
        Sampling rate in Hz
        
    .. cpp:member:: gsl_vector* Pab 
    
        **PAB** column in the library 
        
    .. cpp:member:: gsl_vector* PRCLWN 
    
        **PCLx** column in the library

    .. cpp:member:: gsl_vector* PRCLOFWM 
    
        **OFWx** column in the library
        
    .. cpp:member:: double* calculatedEnergy
    
        Calculated energy in eV.
        
    .. cpp:member:: int numlags
    
        Number of lags (if option:`EnergyMethod` = **OPTFILT** or **I2R** or **I2RALL** or **I2RNOL** or **I2RFITTED** and :option:`OFNoise` = **NSD**)
        
    .. cpp:member:: double tstartNewDev
    
        Addional deviation of the starting time (if :option:`LagsOrNot` = 1)    
        
    .. cpp:member:: int productSize
    
        Size of the scalar product to be calculated
        
    .. cpp:member:: int tooshortPulse_NoLags
    
        Pulse too short to apply lags (1) or not (0)
        
      
.. cpp:function:: int calculateIntParams(ReconstructInitSIRENA *reconstruct_init, int indexa, int indexb, double samprate, int runF0orB0val, gsl_matrix *modelsaux, gsl_matrix *covarianceaux, gsl_matrix *weightaux, gsl_vector *energycolumn, gsl_matrix **Wabaux, gsl_matrix **TVaux, gsl_vector **tEcolumn, gsl_matrix **XMaux, gsl_matrix **YVaux, gsl_matrix **ZVaux, gsl_vector **rEcolumn, gsl_matrix **Pabaux, gsl_matrix **Dabaux, gsl_matrix **PrecalWMaux, gsl_matrix **optimalfiltersabFREQaux, gsl_matrix **optimalfiltersabTIMEaux, gsl_matrix *modelsMaxLengthFixedFilteraux, gsl_matrix **PabMaxLengthFixedFilteraux)
    
    Located in file: *tasksSIRENA.cpp*

    This function calculates some intermediate scalars, vectors and matrices (WAB, TV, tE, XM, YV, ZV, rE, PAB and DAB) for the interpolation and covariance methods. See :ref:`covMatrices` reconstruction method. It is used in :cpp:func:`readAddSortParams` .

    - Declare variables and allocate GSL vectors and matrices
    - Calculate intermediate scalars, vectors and matrices 
    - Free allocated GSL vectors and matrices
    
    **Members/Variables**
        
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    int **indexa** 
    
        Lower index of the library to calculate the intermediate params (:math:`\alpha`)
        
    int **indexb**
    
        Higher index of the library to calculate the intermediate params (:math:`\beta`)
     
    double **samprate**

        Sampling rate
        
    int **runF0orB0val**

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    gsl_matrix* **modelsaux**
        
        GSL input matrix with model template 
    
    gsl_matrix* **covarianceaux**
    
        GSL input matrix with covariance matrix
        
    gsl_matrix* **weightaux**
    
        GSL input matrix with weight matrix 
        
    gsl_vector* **energycolumn**
    
        GSL input vector with list of energies
    
    gsl_matrix** **WAB**
    
        Input/output intermediate parameter
        
    gsl_matrix** **TVaux**
    
        Input/output intermediate parameter
        
    gsl_vector** **tEcolumn**
        
        Input/output intermediate parameter
        
    gsl_matrix** **XMaux**
    
        Input/output intermediate parameter
        
    gsl_matrix** **YVaux**
        
        Input/output intermediate parameter
        
    gsl_matrix** **ZVaux**
        
        Input/output intermediate parameter
        
    gsl_vector** **rEcolumn** 
        
        Input/output intermediate parameter
        
    gsl_matrix** **Pabaux**
    
        Input/output intermediate parameter
        
    gsl_matrix** **Dabaux** 

        Input/output intermediate parameter
        
    gsl_matrix** **precalWMaux**
    
        Input/output intermediate parameter
        
    gsl_matrix** **optimalfiltersabFREQaux**
    
        Input/output intermediate parameter
        
    gsl_matrix** **optimalfiltersabTIMEaux**
    
        Input/output intermediate parameter

    gsl_matrix* **modelsMaxLengthFixedFilteraux**
        
        Input/output intermediate parameter
        
    gsl_matrix** **PabMaxLengthFixedFilteraux**
    
        Input/output intermediate parameter
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: int indexa 
    
        Lower index of the library to calculate the intermediate params (:math:`\alpha`)
        
    .. cpp:member:: int indexb
    
        Higher index of the library to calculate the intermediate params (:math:`\beta`)
     
    .. cpp:member:: double samprate

        Sampling rate
        
    .. cpp:member:: int runF0orB0val

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    .. cpp:member:: gsl_matrix* modelsaux
        
        GSL input matrix with model template 
    
    .. cpp:member:: gsl_matrix* covarianceaux
    
        GSL input matrix with covariance matrix
        
    .. cpp:member:: gsl_matrix* weightaux
    
        GSL input matrix with weight matrix 
        
    .. cpp:member:: gsl_vector* energycolumn
    
        GSL input vector with list of energies
    
    .. cpp:member:: gsl_matrix** WAB
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix** TVaux
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_vector** tEcolumn
        
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **XMaux
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **YVaux
        
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **ZVaux
        
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_vector **rEcolumn 
        
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **Pabaux
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **Dabaux 

        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **precalWMaux
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **optimalfiltersabFREQaux
    
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **optimalfiltersabTIMEaux
    
        Input/output intermediate parameter

    .. cpp:member:: gsl_matrix* modelsMaxLengthFixedFilteraux
        
        Input/output intermediate parameter
        
    .. cpp:member:: gsl_matrix **PabMaxLengthFixedFilteraux
    
        Input/output intermediate parameter


.. cpp:function:: int calculateTemplate(ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight, int inputPulseLength, gsl_vector **pulseaverageMaxLengthFixedFilter)    
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the template (**PULSE** column in the library) of non piled-up pulses.
    Just in case in the detection process some piled-up pulses have not been distinguished as different pulses, a pulseheights histogram is built. This function uses the pulseheights histogram (built by using the **PHEIGHT** column of the library), **Tstart** and **quality** to select the non piled-up pulses.

    1) Declare and initialize variables
    
    2) Before building the histogram, select the pulseheights of the pulses well separated from other pulses whose *quality* = 0
    
    3) Create the pulseheights histogram
    
    4) Calculate the pulseaverage only taking into account the valid pulses:
    
        * Check if the pulse is piled-up or not
        
        * Non piled-up pulses :math:`\Rightarrow` Average them 
    
    5) Calculate covariance and weight matrices
    
    6) Free allocated GSL vectors
    
    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  
        
    PulsesCollection* **pulsesAll** 

        Collection of pulses found in the previous records
        
    PulsesCollection* **pulsesInRecord**

        Collection of pulses found in the current record
        
    double **samprate**

        Sampling rate
        
    gsl_vector** **pulseaverage**

        GSL vector with the pulseaverage (template) of the non piled-up pulses
        
    double* **pulseaverageHeight**

        Height value of the pulseaverage
        
    gsl_matrix** **covariance**

        GSL matrix with covariance matrix
        
    gsl_matrix** **weight**

        GSL matrix with weight matrix (inverse of covariance matrix) 
        
    int **inputPulseLength**

        :option:`PulseLength` input parameter
    
    gsl_vector** **pulseaverageMaxLengthFixedFilter**

        GSL vector with the pulseaverage (template) whose length is :option:`largeFilter` of the non piled-up pulses

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  
        
    .. cpp:member:: PulsesCollection* pulsesAll 

        Collection of pulses found in the previous records
        
    .. cpp:member:: PulsesCollection* pulsesInRecord

        Collection of pulses found in the current record
        
    .. cpp:member: double samprate

        Sampling rate
        
    .. cpp:member:: gsl_vector** pulseaverage

        GSL vector with the pulseaverage (template) of the non piled-up pulses
        
    .. cpp:member:: double* pulseaverageHeight

        Height value of the pulseaverage
        
    .. cpp:member:: gsl_matrix** covariance

        GSL matrix with covariance matrix
        
    .. cpp:member:: gsl_matrix** weight

        GSL matrix with weight matrix (inverse of covariance matrix) 
        
    .. cpp:member: int inputPulseLength

        :option:`PulseLength` input parameter
	
    .. cpp:member:: gsl_vector** pulseaverageMaxLengthFixedFilter

        GSL vector with the pulseaverage (template) whose length is :option:`largeFilter` of the non piled-up pulses
        
        
.. cpp:function::  int calculus_optimalFilter(int TorF, int intermediate, int opmode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex)
    
    Located in file: *tasksSIRENA.cpp*
    
    See description also at :ref:`optimal filter chapter <optimalFilter_NSD>`
    
    This function calculates the optimal filter for a pulse whose matched filter (normalized template) is provided as input 
    parameter, :cpp:member:`matchedfiltergsl`. An optimal filter is just a matched filter that has been adjusted based on the 
    noise spectrum of the system.
    
    It is assumed that all pulses are scaled versions of a template. In the frequency domain (as noise can be frequency dependent), the raw data
    can be expressed as :math:`P(f)=E\cdot S(f)+N(f)`, where :math:`S(f)` is the normalized model pulse shape in the frequency domain, 
    :math:`N(f)` is the power spectrum of the noise and :math:`E` is the scalar amplitude for the photon energy.  
    
    The second assumption is that the noise is stationary, i.e., it does not vary with time. The amplitude of each pulse can then be estimated by 
    minimizing (weighted least-squares sense) the difference between the noisy data and the model pulse shape, being the :math:`\chi^2` condition 
    to be minimized:
    
    .. math::

        \chi^2 = \int \frac{(P(f)-E \cdot S(f))^2}{\langle\lvert N(f)\lvert ^2\rangle} df
         
    In the time domain, the amplitude is the best weighted (optimally filtered) sum of the values in the pulse

    .. math::

       E = k \int p(t)\cdot of(t)

    where :math:`of(t)` is the time domain expression of optimal filter which in frequency domain

    .. math::

        OF(f) = \frac{S^*(f)}{\langle\lvert N(f)\lvert ^2\rangle}

    and :math:`k` is the normalization factor to give :math:`E` in units of energy

    .. math:: 

        k = \int \frac{S(f)\cdot S^{*}(f)}{\langle\lvert N(f)\lvert ^2\rangle} df
             
    Steps:
    
    - FFT calculus of the matched filter (filter template)
    
        - Declare variables
        - Complex FFT values for positive and negative frequencies
        - FFT calculus
        - Generation of the frequencies (positive and negative)
        - Magnitude and argument for positive and negative frequencies
        - Free allocated GSL vectors
        
    - :math:`N(f)`
    - To divide :math:`MatchedFilter(f)/N^2(f)` :math:`\Rightarrow` :math:`MatchedFilter(f)` and :math:`N(f)` must have the same number of points
    
        - *if* (:cpp:member:`mf_size` < *freqgsl->size*) 
    
            - *if* ((*freqgsl->size)%mf_size* == 0) :math:`\Rightarrow` Decimate noise samples
            - *else* :math:`\Rightarrow` It is necessary to work only with the positive frequencies so as not to handle the :math:`f=0` :math:`\Rightarrow` :math:`N(f)` interpolation (:cpp:func:`interpolatePOS`)
                             
        - *else if* (:cpp:member:`mf_size` > *freqgsl->size*) :math:`\Rightarrow` Error: Noise spectrum must have more samples than pulse spectrum
        - *else if* (:cpp:member:`mf_size` == *freqgsl->size*) :math:`\Rightarrow` It is not necessary to do anything
    - :math:`OptimalFilter = MatchedFilter'(f)/N^2(f)`
    - Calculus of the normalization factor
    - Apply the normalization factor
    - Inverse FFT (to get the expression of the optimal filter in time domain)
    
        - Complex :math:`OptimalFilter(f)` :math:`\Rightarrow` Taking into account magnitude :math:`MatchedFilter(f)/N^2(f)` and phase given by :math:`MatchedFilter(f)`
    - Free allocated GSL vectors
        
    **Members/Variables**
    
    int **TorF**

        If :option:`FilterDomain` = **T** :math:`\Rightarrow` :cpp:member:`TorF` = 0; If :option:`FilterDomain` = **F** :math:`\Rightarrow` :cpp:member:`TorF` = 1
        
    int **intermediate** 

        If :option:`intermediate` = 0 :math:`\Rightarrow` Do not write an intermediate file; If :option:`intermediate` = 1 :math:`\Rightarrow` Write an intermediate file
        
    int **opmode**

        If :option:`opmode` = 0 :math:`\Rightarrow` CALIBRATION run (library creation); If :option:`opmode` = 1 :math:`\Rightarrow` RECONSTRUCTION run (energy determination)

    gsl_vector* **matchedfiltergsl** 

        Matched filter associated to the pulse (in general, from the interpolation between two matched filters of the library)
        
    long **mf_size** 

        Matched filter size (samples)
        
    double **samprate**

        Sampling rate
        
    int **runF0orB0val**

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0.
        
    gsl_vector* **freqgsl**

        Frequency axis of the current noise spectral density (input)
        
    gsl_vector* **csdgsl** 

        Current noise spectral density (input)
        
    gsl_vector* * **optimal_filtergsl** 

        Optimal filter in time domain (output)
        
    gsl_vector** **of_f** 

        Frequency axis of the optimal filter spectrum (output)
        
    gsl_vector** **of_FFT** 

        Optimal filter spectrum (absolute values) (output)

    gsl_vector_complex** **of_FFT_complex** 

        Optimal filter spectrum (complex values) (output)

    .. cpp:member:: int TorF

        If :option:`FilterDomain` = **T** :math:`\Rightarrow` :cpp:member:`TorF` = 0; If :option:`FilterDomain` = **F** :math:`\Rightarrow` :cpp:member:`TorF` = 1
        
    .. cpp:member:: int intermediate 

        If :option:`intermediate` = 0 :math:`\Rightarrow` Do not write an intermediate file; If :option:`intermediate` = 1 :math:`\Rightarrow` Write an intermediate file
        
    .. cpp:member:: int opmode

        If :option:`opmode` = 0 :math:`\Rightarrow` CALIBRATION run (library creation); If :option:`opmode` = 1 :math:`\Rightarrow` RECONSTRUCTION run (energy determination)

    .. cpp:member:: gsl_vector* matchedfiltergsl 

        Matched filter associated to the pulse (in general, from the interpolation between two matched filters of the library)
        
    .. cpp:member:: long mf_size 

        Matched filter size (samples)
        
    .. cpp:member:: double samprate

        Sampling rate
        
    .. cpp:member:: int runF0orB0val

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0.
        
    .. cpp:member:: gsl_vector* freqgsl

        Frequency axis of the current noise spectral density (input)
        
    .. cpp:member:: gsl_vector* csdgsl 

        Current noise spectral density (input)
        
    .. cpp:member:: gsl_vector* * optimal_filtergsl 

        Optimal filter in time domain (output)
        
    .. cpp:member:: gsl_vector** of_f 

        Frequency axis of the optimal filter spectrum (output)
        
    .. cpp:member:: gsl_vector** of_FFT 

        Optimal filter spectrum (absolute values) (output)

    .. cpp:member:: gsl_vector_complex** of_FFT_complex 

        Optimal filter spectrum (complex values) (output)
        
        
.. cpp:function:: int convertI2R (char* EnergyMethod, double R0, double Ibias, double Imin, double Imax, double TTR, double LFILTER, double RPARA, double samprate, gsl_vector **invector)
    
    Located in file: *tasksSIRENA.cpp*
    
    This funcion converts the current space into a quasi-resistance space (see :ref:`rSpace` for **I2R**, **I2RALL**, **I2RNOL** and **I2RFITTED** modes). The input :cpp:member:`invector` filled in with current values is filled in here with *I2R*, *I2RALL*, *I2RNOL* or **I2RFITTED** quasi-resistances at the output.

    If :cpp:member:`invector` contains the **ADC** column data from the input FITS file and :math:`I = I_{bias}-(invector*ADUCNV+I_{min})`: 
    
    - Conversion according to :option:`EnergyMethod` = **I2R**:
        
        :math:`DeltaI = I-I_{bias}`    :math:`R/R0 = 1 - (abs(DeltaI)/I_{bias})/(1+abs(DeltaI)/I_{bias})`
        
    - Conversion according to :option:`EnergyMethod` = **I2RALL**:
    
        :math:`R = (V_0-I*R_L-LdI/dt)/I`

    - Conversion according to :option:`EnergyMethod` = **I2RNOL** (*I2RALL* neglecting the circuit inductance):
    
        :math:`R = (V_0-I*R_L)/I`
    
    - Conversion according to :option:`EnergyMethod` = **I2RFITTED**  

        :math:`R = V_0/(I_{fit}+I)`
    
    
    **Members/Variables**

    char* **EnergyMethod**
    
        Quasi-resistance energy calculation method: **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED**, :option:`EnergyMethod`

    double **R0**

        Operating point resistance
        
    double **Ibias**

        Initial bias current
        
    double **Imin**

        Current corresponding to 0 ADU
    
    double **Imax**

        Current corresponding to maximum ADU
        
    double **TTR**

        Transformer Turns Ratio
        
    double **LFILTER**

        Filter circuit inductance
        
    double **RPARA**

        Parasitic resistor value
        
    double **samprate**

        Sampling rate
        
    gsl_vector* **invector**

        GSL vector with input signal values (**ADC** column of the input FITS file)  

    .. cpp:member:: char* EnergyMethod
    
        Quasi-resistance energy calculation method: **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED**, :option:`EnergyMethod`

    .. cpp:member:: double R0

        Operating point resistance
        
    .. cpp:member:: double Ibias

        Initial bias current
        
    .. cpp:member:: double Imin

        Current corresponding to 0 ADU
    
    .. cpp:member:: double Imax

        Current corresponding to maximum ADU
        
    .. cpp:member:: double TTR

        Transformer Turns Ratio
        
    .. cpp:member:: double LFILTER

        Filter circuit inductance
        
    .. cpp:member:: double RPARA

        Parasitic resistor value
        
    .. cpp:member:: double samprate

        Sampling rate
        
    .. cpp:member:: gsl_vector* invector

        GSL vector with input signal values (**ADC** column of the input FITS file)  
        
        
.. cpp:function:: int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, fitsfile **dtcObject, int inputPulselength)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function creates an intermediate FITS file with some useful info (during the development phase) if the :option:`intermediate` input parameter is set to 1.

    The intermediate FITS file will contain 2 HDUs:
    
        * *PULSES* HDU will contain some info about the found pulses: **TSTART**, **I0** (the pulse itself), **TEND**, **TAURISE**, **TAUFALL** and **QUALITY**
        
        * *TESTINFO* HDU will contain columns **FILDER** (the low-pass filtered and differentiated records) and **THRESHOLD**

    If file exists :math:`\Rightarrow` Check :option:`clobber` for overwritting. If it does not, then create it.
        
    **Members/Variables**
        
    ReconstructInitSIRENA** **reconstruct_init**
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
                
    double **samprate**
        
        Sampling rate 
            
    fitsfile **dtcObject**
        
        Object which contains information of the intermediate FITS file (used also by :cpp:func:`writeTestInfo` and :cpp:func:`writePulses`).

    int **inputPulseLength**
        
        :option:`PulseLength` input parameter

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
                
    .. cpp:member:: double samprate
        
        Sampling rate 
            
    .. cpp:member:: fitsfile dtcObject
        
        Object which contains information of the intermediate FITS file (used also by :cpp:func:`writeTestInfo` and :cpp:func:`writePulses`).
	
    .. cpp:member:: int inputPulseLength
        
        :option:`PulseLength` input parameter    

        
.. cpp:function:: int createHisto(gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
        
    Located in file: *tasksSIRENA.cpp*
    
    This function builds the histogram of the input vector.

      - Histogram x-axis values are the different input vector values (pulseheights)
      
      - Histogram y-axis values are the the number of cases per unit of the variable on the horizontal axis

    1) Declare variables
    
    2) It will work with the positive elements of the input vector :math:`\Rightarrow` *invectoraux2*

    3) Check if all the values of :cpp:member:`invector` are the same :math:`\Rightarrow` Histogram of only one bin

    4) Obtain *invector_max* and *invector_min*
    
    5) Obtain *binSize*
    
    6) Create histogram axis

    7) Free allocated GSL vectors

    **Members/Variables**
    
    gsl_vector* **invector**

        GSL input vector

    int **nbins**

        Number of bins to build the histogram

    gsl_vector** **xhistogsl**

        GSL vector with output histogram x-axis

    gsl_vector** **yhistogsl**

        GSL vector with output histogram y-axis
        
    .. cpp:member:: gsl_vector* invector

        GSL input vector

    .. cpp:member:: int nbins

        Number of bins to build the histogram

    .. cpp:member:: gsl_vector** xhistogsl

        GSL vector with output histogram x-axis

    .. cpp:member:: gsl_vector** yhistogsl

        GSL vector with output histogram y-axis
            
                            
.. cpp:function:: int createLibrary(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject, int inputPulseLength)
    
    Located in file: *tasksSIRENA.cpp*

    This function creates the pulse templates library FITS file, if it does not exist yet. Otherwise, it opens it (to add a new row).

        1) If it exists :math:`\Rightarrow` Open it and set *appendToLibrary = true*
        
        2) If it does not exist :math:`\Rightarrow` Create it and set *appendToLibrary = false*

            - Write keyword ``EVENTCNT`` = 1 in the *LIBRARY* extension
            - Write the whole list of input parameters in ``HISTORY`` in the *Primary* extension (by usin 'HDpar_stamp')
            
    **Members/Variables**
            
    ReconstructInitSIRENA** **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    bool **appendToLibrary**

        Used by the function :cpp:func:`writeLibrary`
        
    fitsfile** **inLibObject**

        Object which contains information of the library FITS file (used also by :cpp:func:`writeLibrary`)
	
    int **inputPulseLength**
        
        :option:`PulseLength` input parameter
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    .. cpp:member:: bool appendToLibrary

        Used by the function :cpp:func:`writeLibrary`
        
    .. cpp:member:: fitsfile** inLibObject

        Object which contains information of the library FITS file (used also by :cpp:func:`writeLibrary`)
	
    .. cpp:member:: int inputPulseLength
        
        :option:`PulseLength` input parameter
        
.. cpp:function:: int createTPSreprFile ()

    Located in file: *gennoisespec.cpp*

    This function creates the gennoisespec output FITS file.
    
    Steps:
    
    - Create the noise representation file (if it does not exist already)
    - Create the extensions *NOISE*, *NOISEALL* and *WEIGHTMS*
    - Write keywords
    
    
.. _D:

.. cpp:function:: int differentiate(gsl_vector **invector, int szVct)

    Located in file: *pulseprocess.cpp*

    This function applies the derivative method :math:`x_i-x_{i-1}` to the input vector.

    The derivative method provides more sensitivity to handle with piled-up pulses.
    Moreover, little variations of the baseline will not affect.

    **Members/Variables**
    
    gsl_vector** **invector** 
    
        Input/Ouput GSL vector (non-differentiate input vector/differentiate input vector)
        
    int **szVct**
    
        Size of :cpp:member:`invector`
        
    .. cpp:member:: gsl_vector** invector 
    
        Input/Ouput GSL vector (non-differentiate input vector/differentiate input vector)
        
    .. cpp:member:: int szVct
    
        Size of :cpp:member:`invector`
        

.. _E:

.. cpp:function:: int eigenVV (gsl_matrix *matrixin, gsl_matrix **eigenvectors, gsl_vector **eigenvalues)
    
    Located in file: *tasksSIRENA.cpp*
    
    This funcion provides the principal eigenvectors and eigenvalues of the input matrix (at the moment, the first two eigenvalues and eigenvectors). The eigenvalues and eigenvectors are sorted in descending order and only the principal components are provided.

    - Calculate the eigenvectors and the eigenvalues
    - Sort the eigenvectors and the eigenvalues in descending order
    - Choose the main eigenvectors and eigenvalues (the principal components analysis). At the moment, the first two eigenvectors and eigenvalues

    **Members/Variables**

    gsl_matrix* **matrixin**

        Input GSL matrix

    gsl_matrix** **eigenvectors**

        Subset of eigenvectors of 'matrixin' chosen by PCA (the first two ones)
        
    gsl_vector** **eigenvalues**

        Subset of eigenvalues of 'matrixin' chosen by PCA (the first two ones)
        
    .. cpp:member:: gsl_matrix* matrixin

        Input GSL matrix

    .. cpp:member:: gsl_matrix** eigenvectors

        Subset of eigenvectors of 'matrixin' chosen by PCA (the first two ones)
        
    .. cpp:member:: gsl_vector** eigenvalues

        Subset of eigenvalues of 'matrixin' chosen by PCA (the first two ones)
        

.. cpp:function:: void exit_error(const char* const func, string msg, int status)

    Located in file: *genutils.cpp*

    This function prints out error messages and exits program.

    **Members/Variables**
    
    const char* const **func**
    
        Function name whose error is printed 
        
    string **msg**
    
        Error message to be printed

    int **status**
        
        Status
        
    .. cpp:member:: const char* const func
    
        Function name whose error is printed 
        
    .. cpp:member:: string msg
    
        Error message to be printed

    .. cpp:member::  int status
        
        Status 

        
.. _F:

.. cpp:function:: int FFT(gsl_vector *invector, gsl_vector_complex *outvector, double STD)
    
    Located in file: *genutils.cpp*
    
    This function calculates the FFT of the elements of a vector.

    GSL library (overview of FFTs):

    For physical applications it is important to remember that the index appearing in the DFT does not correspond directly to a physical frequency. If the time-step of the
    DFT is :math:`\Delta` then the frequency domain includes both positive and negative frequencies, ranging from :math:`-1/(2\Delta)` through 0 to :math:`+1/(2\Delta)`. The positive frequencies are stored from the beginning of the array up to the middle, and the negative frequencies are stored backwards from the end of the array.

    Here is a table which shows the layout of the array data, and the correspondence between the time domain data z, and the frequency domain data x.

    =======   ==================   =========================================
     index         z                        x = FFT(z)
    =======   ==================   =========================================
     0        :math:`z(t = 0)`     :math:`x(f = 0)`
     1        :math:`z(t = 1)`     :math:`x(f = 1/(n\Delta))`
     2        :math:`z(t = 2)`     :math:`x(f = 2/(n\Delta))`
     [...]        [........]             [..................]
     n/2      :math:`z(t = n/2)`   :math:`x(f = +1/(2\Delta),-1/(2\Delta))`
     [...]        [........]             [..................]
     n-3      :math:`z(t = n-3)`   :math:`x(f = -3/(n\Delta))`
     n-2      :math:`z(t = n-2)`   :math:`x(f = -2/(n\Delta))`
     n-1      :math:`z(t = n-1)`   :math:`x(f = -1/(n\Delta))`
    =======   ==================   =========================================
    
    The frequency axis will be built as *f = i/STD = i/(size/samprate)* with *i* varying from 0 to *size/2-1* (*n=size* and :math:`\Delta=1/samprate`  sec/sample).

    **Members/Variables**
    
    gsl_vector* **invector**
    
        Input GSL vector
        
    gsl_vector_complex* **outvector**
    
        Output GSL complex vector with the FFT of :cpp:member:`invector`

    double **STD**
        
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*
        
    .. cpp:member:: gsl_vector* invector
    
        Input GSL vector
        
    .. cpp:member:: gsl_vector_complex* outvector
    
        Output GSL complex vector with the FFT of :cpp:member:`invector`

    .. cpp:member::  double STD
        
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*
      
      
.. cpp:function:: int FFTinverse(gsl_vector_complex *invector, gsl_vector *outvector, double STD)
    
    Located in file: *genutils.cpp*
    
    This function calculates the inverse FFT of the elements of a vector.
    
    **Members/Variables**
    
    gsl_vector_complex* **invector**
    
        Input GSL complex vector
        
    gsl_vector* **outvector**
    
        Output GSL vector with the inverse FFT of :cpp:member:`invector`
    
    double **STD**
    
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*
        
    .. cpp:member:: gsl_vector_complex* invector
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector* outvector
    
        Output GSL vector with the inverse FFT of :cpp:member:`invector`
    
    .. cpp:member:: double STD
    
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*

        
.. cpp:function:: int filderLibrary (ReconstructInitSIRENA** reconstruct_init, double samprate)
        
    Located in file: *tasksSIRENA.cpp*

    This function calculates the (low-pass filter and) derivative of the models (*pulse_templates*) in the library (only necessary if first record), 
    and it stores the *pulse_templates_filder* and the *maxDERs* and *samp1DERs* in the :cpp:member:`reconstruct_init` structure.

    The maximum of the (low-pass filtered and) differentiated pulse has to be compared to the *maxDERs* to select the appropriate model. Or, the 1st sample out of the differentiated pulse has to be compared to the *samp1DERs* to select the appropriate model.

    1) Check if it is the first record
    
    2) (Low-pass filter and) differentiate the models (*pulse_templates*) of the library
    
    3) Store the (low-pass filtered) derivatives in *pulse_templates_filder*
    
    4) Calculate the maximum of the (low-pass filtered and) differentiated models (*maxDERs*)
    
    5) Locate the 1st sample of the (low-pass filtered and) differentiated models (*samp1DERs*)

    **Members/Variables**

    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 

    double **samprate**

        Sampling rate
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 

    .. cpp:member:: double samprate

        Sampling rate
            
            
.. cpp:function:: bool fileExists(const std::string& name)
    
    Located in file: *genutils.cpp*
    
    This function checks for file existence returning a boolean value.
    
    **Members/Variables**
    
    const std::string& **name**
         
        File name
        
    .. cpp:member:: const std::string& name
         
        File name
  
  
.. cpp:function:: int filterByWavelets (ReconstructInitSIRENA* reconstruct_init, gsl_vector **invector, int length, int *onlyOnce)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function filters the input/output signal :cpp:member:`invector`, reducing the noise level.
    
    Steps:
    
    - It is only going to work with *n* elements of :cpp:member:'invector'
    - Discrete Wavelet Transform 
    - Sorting coefficients
    - Hard thresholding: *n-nc* coefficients are deleted (those with low energy)
    - Inverse DWT
        
    **Members/Variables**
        
    ReconstructInitSIRENA** **reconstruct_init**
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
                
    gsl_vector** **invector**
        
        Input/output signal 

    int **length**
        
        Length of the wavelet transform
        
    int* **onlyOnce**
        
        In order to control the times to be executed
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
                
    .. cpp:member:: gsl_vector** invector
        
        Input/output signal 

    .. cpp:member:: int length
        
        Length of the wavelet transform
        
    .. cpp:member:: int *onlyOnce
        
        In order to control the times to be executed
        
        
.. cpp:function:: int findInterval(int tail_duration, gsl_vector *invector, gsl_vector *startpulse, int npin, int pulse_length, int nPF, int interval, int *ni, gsl_vector **startinterval)

    Located in file: *gennoisespec.cpp*
    
    This function finds the pulse-free intervals when the input vector has pulses.
    The pulse-free intervals must have a minimum length (*intervalMinBins*).
    The interval with pulse is :math:`Tstart,Tend+nPF*pulse \_ length` (being :math:`Tend=n*pulse \_ length`).

    Steps:
    
    - Declare variables
    - Processing if the input vector has more pulses
        - It looks for pulse-free intervals between pulses
    - Processing if there are no more pulses in the input vector
        - It looks for pulse-free intervals at the end of the event and the search for more pulse-free intervals is finished

    **Members/Variables**
    
    int **tail_duration** 
    
        Length of the tail of a previous pulse
        
    gsl_vector* **invector**
    
        Input vector WITH pulses
        
    gsl_vector* **startpulse** 
    
        Vector with the Tstart of all the pulses of the input vector (samples)
        
    int **npin** 
    
        Number of pulses in the input vector
        
    int **pulse_length** 
    
        Pulse length (samples)
        
    int **nPF** 
    
        Number of pulse lengths after ending the pulse to start the pulse-free interval
        
    int **interval** 
    
        Minimum length of the interval (samples)
        
    int **ni** 
    
        Number of pulse-free intervals in the input vector
        
    gsl_vector** **startinterval** 
    
        Vector with the starting time of each pulse-free interval (samples)
        
    .. cpp:member:: int tail_duration 
    
        Length of the tail of a previous pulse
        
    .. cpp:member:: gsl_vector* invector
    
        Input vector WITH pulses
        
    .. cpp:member:: gsl_vector* startpulse 
    
        Vector with the Tstart of all the pulses of the input vector (samples)
        
    .. cpp:member:: int npin
    
        Number of pulses in the input vector
        
    .. cpp:member:: int pulse_length 
    
        Pulse length (samples)
        
    .. cpp:member:: int nPF
    
        Number of pulse lengths after ending the pulse to start the pulse-free interval
        
    .. cpp:member:: int interval 
    
        Minimum length of the interval (samples)
        
    .. cpp:member:: int ni 
    
        Number of pulse-free intervals in the input vector
        
    .. cpp:member:: gsl_vector** startinterval 
    
        Vector with the starting time of each pulse-free interval (samples)
        
.. cpp:function:: int findIntervalN (gsl_vector *invector, int interval, int *ni, gsl_vector **startinterval)

    Located in file: *gennoisespec.cpp*
    
    This function finds the pulse-free intervals when the input vector has NO pulses.
    The pulse-free intervals must have a minimum length (*intervalMinBins*).

    **Members/Variables**
    
    gsl_vector* **invector** 
    
        Input vector WITHOUT pulses
        
    int **interval** 
    
        Minimum length of the interval (samples)
        
    int* **ni** 
    
        Number of pulse-free intervals in the input vector
        
    gsl_vector** **startinterval** 
    
        Vector with the starting time of each pulse-free interval (samples)
        
    .. cpp:member:: gsl_vector* invector
    
        Input vector WITHOUT pulses
        
    .. cpp:member:: int interval 
    
        Minimum length of the interval (samples)
        
    .. cpp:member:: int* ni 
    
        Number of pulse-free intervals in the input vector
        
    .. cpp:member:: gsl_vector** startinterval
    
        Vector with the starting time of each pulse-free interval (samples)   
        
        
.. cpp:function:: int findMeanSigma(gsl_vector *invector, double *mean, double *sigma)
    
    Located in file: *pulseprocess.cpp*
    
    This function calculates the mean and the standard deviation of the input vector.
    
    **Members/Variables**
    
    gsl_vector* **invector**
         
        Input GSL vector

    double* **mean**
         
        Mean of the elements of :cpp:member:`invector`
        
    double* **sigma**
         
        Standard deviation of the elements of :cpp:member:`invector`
    
    .. cpp:member:: gsl_vector* invector
         
        Input GSL vector

    .. cpp:member:: double* mean
         
        Mean of the elements of :cpp:member:`invector`
        
    .. cpp:member:: double* sigma
         
        Standard deviation of the elements of :cpp:member:`invector`
        
    
.. cpp:function:: int findPulsesCAL(gsl_vector *vectorin, gsl_vector *vectorinDER, gsl_vector **tstart, gsl_vector **quality, gsl_vector **pulseheight, gsl_vector **maxDERgsl, int *nPulses, double *threshold, double scalefactor, double samplingRate, int samplesup, double nsgms, double lb, double lrs, ReconstructInitSIRENA *reconstruct_init, double stopcriteriamkc, double kappamkc)

    Located in file: *pulseprocess.cpp*

    This function is going to find the pulses in a record (in the *CALibration* mode) by using the function :cpp:func:`findTstartCAL`.

    Steps:
    
    - Declare variables
    - Establish the threshold (call :cpp:func:`medianKappaClipping`)
    - Find pulses (call :cpp:func:`findTstartCAL`)
    - If at least a pulse is found
    
      - Get :cpp:member:`pulseheight` of each found pulse (in order to be used to build the pulse templates library)
      
    - Free allocated GSL vectors

    **Members/Variables**
    
    gsl_vector* **vectorin**
    
        Not filtered record
        
    gsl_vector* **vectorinDER**
    
        Derivative of the (low-pass filtered) :cpp:member:`vectorin`

    gsl_vector** **tstart**
        
        Starting time of the found pulses into the record (in samples)
        
    gsl_vector** **quality**
    
        Quality of the found pulses into the record
        
    gsl_vector** **pulseheight**
    
        Pulse height of the found pulses into the record

    gsl_vector** **maxDERgsl**
        
        Maximum of the derivative of the found (low-pass filtered) pulses into the record
        
    int* **nPulses**
    
        Number of found pulses
        
    double* **threshold**
    
        Threshold used to find the pulses (output parameter because it is necessary out of the function)

    double **scalefactor**
    
        Scale factor to calculate the LPF box-car length (:option:`scaleFactor`)
        
    double **samplingRate**
    
        Sampling rate

    int **samplesup**
        
        Number of consecutive samples over the threshold to locate a pulse (:option:`samplesUp`)
        
    double **nsgms**
    
        Number of Sigmas to establish the threshold (:option:`nSgms`)
        
    double **lb**
    
        Vector containing the baseline averaging length used for each pulse

    double **lrs**
        
        Running sum length (:option:`LrsT` in samples)
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    double **stopcriteriamkc**
    
        Used in :cpp:func:`medianKappaClipping` (%)

    double **kappamkc**
        
        Used in :cpp:func:`medianKappaClipping`

    .. cpp:member:: gsl_vector* vectorin
    
        Not filtered record
        
    .. cpp:member:: gsl_vector* vectorinDER
    
        Derivative of the (low-pass filtered) :cpp:member:`vectorin`

    .. cpp:member:: gsl_vector** tstart
        
        Starting time of the found pulses into the record (in samples)
        
    .. cpp:member:: gsl_vector** quality
    
        Quality of the found pulses into the record
        
    .. cpp:member:: gsl_vector** pulseheight
    
        Pulse height of the found pulses into the record

    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the derivative of the found (low-pass filtered) pulses into the record
        
    .. cpp:member:: int* nPulses
    
        Number of found pulses
        
    .. cpp:member:: double* threshold
    
        Threshold used to find the pulses (output parameter because it is necessary out of the function)

    .. cpp:member:: double scalefactor
    
        Scale factor to calculate the LPF box-car length (:option:`scaleFactor`)
        
    .. cpp:member:: double samplingRate
    
        Sampling rate

    .. cpp:member:: int samplesup
        
        Number of consecutive samples over the threshold to locate a pulse (:option:`samplesUp`)
        
    .. cpp:member:: double nsgms
    
        Number of Sigmas to establish the threshold (:option:`nSgms`)
        
    .. cpp:member:: double lb
    
        Vector containing the baseline averaging length used for each pulse

    .. cpp:member:: double lrs
        
        Running sum length (:option:`LrsT` in samples)
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    .. cpp:member:: double stopcriteriamkc
    
        Used in :cpp:func:`medianKappaClipping` (%)

    .. cpp:member:: double kappamkc
        
        Used in :cpp:func:`medianKappaClipping`
        
        
.. cpp:function:: int findPulsesNoise(gsl_vector *vectorin, gsl_vector *vectorinDER, gsl_vector **tstart, gsl_vector **quality, gsl_vector **energy, int *nPulses, double *threshold, double scalefactor, int sizepulsebins, double samplingRate, int samplesup, double nsgms, double lb, double lrs, double stopcriteriamkc, double kappamkc)

    Located in file: *gennoisespec.cpp*
    
    This function is going to find the pulses in a record by using the function :cpp:func:`findTstartNoise`
    
    Steps:
    
    - Declare variables
    - Establish the threshold (call :cpp:func:`medianKappaClipping`)
    - Find pulses (call :cpp:func:`findTstartNoise`)
    - If at least a pulse is found
        - Get the *pulseheight* of each found pulse
    - Free allocated GSL vectors
    
    **Members/Variables**
    
    gsl_vector* **vectorin**
    
        Not filtered record
        
    gsl_vector* **vectorinDER** 
    
        Derivative of the low-pass filtered :cpp:member:`vectorin`
        
    gsl_vector** **tstart**
    
        Starting time of the found pulses into the record (samples)
    
    gsl_vector** **quality**
    
        Quality of the found pulses into the record
        
    int* **nPulses**
    
        Number of found pulses
        
    double* **threshold**
    
        Threshold used to find the pulses (output parameter because it is necessary out of the function)
        
    double **scalefactor**
    
        Scale factor to calculate the LPF box-car length
        
    int **sizepulsebins**
    
        Size of the pulse (samples)
    
    double **samplingRate**
    
        Sampling rate

    int **samplesup**
    
        Number of consecutive samples over the threshold to locate a pulse
        
    double **nsgms**
    
        Number of Sigmas to establish the threshold

    double **lb**
    
        Vector containing the baseline averaging length used for each pulse
        
    double **lrs**
    
        Running sum length
        
    double **stopCriteriamkc**
    
        Used in :cpp:func:`medianKappaClipping_noiseSigma` (%)

    double **kappamkc**
    
        Used in :cpp:func:`medianKappaClipping_noiseSigma`
        
    .. cpp:member:: gsl_vector* vectorin
    
        Not filtered record
        
    .. cpp:member:: gsl_vector* vectorinDER 
    
        Derivative of the low-pass filtered 'vectorin'
        
    .. cpp:member:: gsl_vector** tstart
    
        Starting time of the found pulses into the record (samples)
    
    .. cpp:member:: gsl_vector** quality
    
        Quality of the found pulses into the record
        
    .. cpp:member:: int* nPulses
    
        Number of found pulses
        
    .. cpp:member:: double* threshold
    
        Threshold used to find the pulses (output parameter because it is necessary out of the function)
        
    .. cpp:member:: double scalefactor
    
        Scale factor to calculate the LPF box-car length
        
    .. cpp:member:: int sizepulsebins
    
        Size of the pulse (samples)
    
    .. cpp:member:: double samplingRate
    
        Sampling rate
  
    .. cpp:member:: int samplesup
        
        Number of consecutive samples over the threshold to locate a pulse (*samplesUp*)
        
    .. cpp:member:: double nsgms
    
        Number of Sigmas to establish the threshold
  
    .. cpp:member:: double lb
    
        Vector containing the baseline averaging length used for each pulse
        
    .. cpp:member:: double lrs
    
        Running sum length (equal to the *Lrs* input parameter)
  
    .. cpp:member:: double stopCriteriamkc
    
        Used in :cpp:func:`medianKappaClipping_noiseSigma` (%)
  
    .. cpp:member:: double kappamkc
    
        Used in :cpp:func:`medianKappaClipping_noiseSigma`
        
  
.. cpp:function:: int FindSecondaries(int maxPulsesPerRecord, gsl_vector *adjustedDerivative, double adaptativethreshold, ReconstructInitSIRENA *reconstruct_init, int tstartFirstEvent, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **lagsgsl)

    Located in file: *pulseprocess.cpp*

    This function runs after :cpp:func:`InitialTriggering` to find all the events (except the first one) in the first derivative of the (low-pass filtered) record by using the Adjusted Derivative detection method.

    Steps: 
    
    - Declare variables
    - Establishing the criteria of the slope of the derivative depending on the sampling rate
    - It is necessary to find the tstarts... 
    
      It looks for an event and if a pulse is found, it looks for another event
    
        - It looks for an event since the beginning (or the previous event) to the end of the record. 
          The first condition to detect an event is that the :cpp:member:`adjustedDerivative` was over the :cpp:member:`threshold`
          
            - Select the model of the found pulse from the libary by using the 1st sample of the derivative (*samp1DER*)
            - Dot product between the detected pulse and the pulse template in 3 different lags
            
                - If maximum of the dot product found :math:`\Rightarrow` Stop calculating dot products in more lags
                - If maximum of the dot product not found :math:`\Rightarrow` Calculate dot products in more lags (number of lags is limited to 5)
            
            - If maximum of the dot product not found :math:`\Rightarrow` tstart is the first sample crossing above the threshold (without jitter)
              
                - Average of the first 4 samples of the derivative
                - Find model in order to subtract
              
            - If maximum of the dot product found :math:`\Rightarrow` Parabola analytically defined :math:`\Rightarrow` Locate the maximum :math:`\Rightarrow` New tstart (with jitter)
            
                - Iterative process in order to extract the best template from the library:
                    - *samp1DER* correction
                    - Find the model from the libary by using the corrected *samp1DER*
                    - Dot product in 3 lags
                    - Locate the maximum of the parabola
                - *samp1DER* correction
                - Find model in order to subtract
                - Template correction
                - Average of the first 4 samples of the derivative
                
            - The second condition to detect an event is meeting the criteria of the slope of the derivative

        - Subtract the model from the adjusted derivative

            - Select the model of the found event from the libary by using the first sample of the derivative 
            - Subtract
    
    - ... Or to use the tstart provided as input parameters
    
        - Obtain the *maxDERs* of the events whose tstarts have been provided (by using the maximum of the derivative to find the model)
        
    - Free allocated GSL vectors

    **Members/Variables**
    
    int **maxPulsesPerRecord**
    
        Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    gsl_vector* **adjustedDerivative**
    
        First derivative of the (low-pass filtered) record

    double **adaptativethreshold**
        
        Threshold
        
    double **samprate** 

        Sampling rate
    
    ReconstructInitSIRENA* **reconstruct_init**
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    int **tstartFirstEvent**
        
        Tstart of the first event of the record (in samples) found by :cpp:func:`InitialTriggering`
    
    int* **numberPulses**
        
        Number of found events
    
    gsl_vector** **tstartgsl**
        
        Starting time of the found events (in samples)
    
    gsl_vector** **flagTruncated**
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    gsl_vector** **maxDERgsl**
        
        Maximum of the derivative of the event 
        
    gsl_vector** **samp1DERgsl**
        
        Average of the first 4 samples of the derivative of the event
        
    gsl_vector** **lagsgsl**
        
        Number of necessary lags to establish the tstart (currently limited to 5)

    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    .. cpp:member:: gsl_vector* adjustedDerivative
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double adaptativethreshold
        
        Threshold
        
    .. cpp:member:: double samprate 

        Sampling rate
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: int tstartFirstEvent
        
        Tstart of the first event of the record (in samples) found by :cpp:func:`InitialTriggering`
    
    .. cpp:member:: int* numberPulses
        
        Number of found events
    
    .. cpp:member:: gsl_vector** tstartgsl
        
        Starting time of the found events (in samples)
    
    .. cpp:member:: gsl_vector** flagTruncated
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the derivative of the event 
        
    .. cpp:member:: gsl_vector** samp1DERgsl
        
        Average of the first 4 samples of the derivative of the event
        
    .. cpp:member:: gsl_vector** lagsgsl
        
        Number of necessary lags to establish the tstart (currently limited to 5)
        
        
.. cpp:function:: int FindSecondariesSTC(int maxPulsesPerRecord, gsl_vector *adjustedDerivative, double adaptativethreshold, ReconstructInitSIRENA *reconstruct_init, int tstartFirstEvent, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **lagsgsl)

    Located in file: *pulseprocess.cpp*

    This function runs after :cpp:func:`InitialTriggering` to find all the events (except the first one) in the first derivative of the (low-pass filtered) record by using the Single Threshold Crossing method.

    Steps: 
    
    - Declare variables
    - It is necessary to find the tstarts... 
    
      It looks for an event and if a pulse is found, it looks for another event
    
        - It looks for an event since the beginning (or the previous event) to the end of the record. 
          The condition to detect an event is that the :cpp:member:`adjustedDerivative` was over the :cpp:member:`threshold` at least :option:`samplesUp` consecutive samples
    
    - ... Or to use the tstart provided as input parameters
    
        - Obtain the *maxDERs* of the events whose tstarts have been provided
        
    - Free allocated GSL vectors

    **Members/Variables**
    
    int **maxPulsesPerRecord**
    
        Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    gsl_vector* **adjustedDerivative**
    
        First derivative of the (low-pass filtered) record

    double **adaptativethreshold**
        
        Threshold
        
    double **samprate** 

        Sampling rate
    
    ReconstructInitSIRENA* **reconstruct_init**
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    int **tstartFirstEvent**
        
        Tstart of the first event of the record (in samples) found by :cpp:func:`InitialTriggering`
    
    int* **numberPulses**
        
        Number of found events
    
    gsl_vector** **tstartgsl**
        
        Starting time of the found events (in samples)
    
    gsl_vector** **flagTruncated**
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    gsl_vector** **maxDERgsl**
        
        Maximum of the derivative of the event 
        
    gsl_vector** **samp1DERgsl**
        
        Average of the first 4 samples of the derivative of the event

    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    .. cpp:member:: gsl_vector* adjustedDerivative
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double adaptativethreshold
        
        Threshold
        
    .. cpp:member:: double samprate 

        Sampling rate
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: int tstartFirstEvent
        
        Tstart of the first event of the record (in samples) found by :cpp:func:`InitialTriggering`
    
    .. cpp:member:: int* numberPulses
        
        Number of found events
    
    .. cpp:member:: gsl_vector** tstartgsl
        
        Starting time of the found events (in samples)
    
    .. cpp:member:: gsl_vector** flagTruncated
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the derivative of the event 
        
    .. cpp:member:: gsl_vector** samp1DERgsl
        
        Average of the first 4 samples of the derivative of the event
        
        
.. cpp:function:: int findTstartCAL(int maxPulsesPerRecord, gsl_vector *der, double adaptativethreshold, int nSamplesUp, ReconstructInitSIRENA *reconstruct_init, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl)

    Located in file: *pulseprocess.cpp*

    This function scans all the values of the derivative of the (low-pass filtered) record until it finds :cpp:member:`nSamplesUp` consecutive 
    values (due to the noise more than 1 value is required) over the threshold. To look for more pulses, it finds :cpp:member:`nSamplesUp` consecutive values
    (due to the noise) under the threshold and then, it starts to scan again.

    Steps: 
    
    - Declare variables
    
    - Allocate GSL vectors
    
    - It is possible to find the tstarts...
    
      - Obtain tstart of each pulse in the derivative:
        
        - If :math:`der_i>threshold` and *foundPulse=false*, it looks for :cpp:member:`nSamplesUp` consecutive samples over the threshold

          - If not, it looks again for a pulse crossing over the threshold
          
          - If yes, a pulse is found (truncated if it is at the beginning)
          
        - If :math:`der_i>threshold` and *foundPulse=true*, it looks for a sample under the threshold
        
          - If not, it looks again for a sample under the threshold
          
          - If yes, it looks for :cpp:member:`nSamplesUp` consecutive samples under the threshold and again it starts to look for a pulse
          
      
    - ... Or to use the tstart provided as input parameters
    
      Obtain the *maxDERs* of the pulses whose tstarts have been provided

    **Members/Variables**
    
    int **maxPulsesPerRecord**
    
        Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    gsl_vector* **der**
    
        First derivative of the (low-pass filtered) record

    double **adaptativethreshold**
        
        Threshold
        
    int **nSamplesUp**
    
        Number of consecutive samples over the threshold to 'find' a pulse (:option:`samplesUp`)
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    int* **numberPulses**
    
        Number of found pulses
        
    gsl_vector** **tstartgsl**
    
        Pulses tstart (in samples)
        
    gsl_vector** **flagTruncated**
    
        Flag indicating if the pulse is truncated 

    gsl_vector** **maxDERgsl**
        
        Maximum of the first derivative of the (low-pass filtered) record inside each found pulse

    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    .. cpp:member:: gsl_vector* der
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double adaptativethreshold
        
        Threshold
        
    .. cpp:member:: int nSamplesUp
    
        Number of consecutive samples over the threshold to 'find' a pulse (:option:`samplesUp`)
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    .. cpp:member:: int* numberPulses
    
        Number of found pulses
        
    .. cpp:member:: gsl_vector** tstartgsl
    
        Pulses tstart (in samples)
        
    .. cpp:member:: gsl_vector** flagTruncated
    
        Flag indicating if the pulse is truncated 

    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the first derivative of the (low-pass filtered) record inside each found pulse
        
        
.. cpp:function:: int findTstartNoise(int maxPulsesPerRecord, gsl_vector *der, double adaptativethreshold, int nSamplesUp, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl)
    
    Located in file: *gennoisespec.cpp*.
    
    This function finds the pulses tstarts in the input vector (first derivative of the filtered record).

    This function scans all values the derivative of the (low-pass filtered) record until it finds *nSamplesUp* consecutive values (due to the noise more than 1 value is 
    required) over the threshold. To look for more pulses, it finds *nSamplesUp* consecutive values (due to the noise) under the threshold and then, it starts to scan again.

    Steps: 
    
    - Declare variables
    - Allocate GSL vectors
    - Obtain tstart of each pulse in the derivative:
        - If :math:`der_i > threshold` and *foundPulse=false*, it looks for *nSamplesUp* consecutive samples over the threshold
            - If not, it looks again for a pulse crossing over the threshold
            - If yes, a pulse is found (truncated if it is at the beginning)
        - If :math:`der_i > threshold` and *foundPulse=true*, it looks for a sample under the threshold
            - If not, it looks again for a sample under the threshold
            - If yes, it looks for *nSamplesUp* consecutive samples under the threshold and again it starts to look for a pulse

    **Members/Variables**
    
    int **maxPulsesPerRecord**
    
        Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
        
    gsl_vector* **der**
    
        First derivative of the (low-pass filtered) record
        
    double **adaptativethreshold**
    
        Threshold
        
    int **nSamplesUp**
    
        Number of consecutive samples over the threshold to 'find' a pulse
        
    int* **numberPulses**
    
        Number of found pulses
        
    gsl_vector** **tstartgsl** 
        
        Pulses tstart (samples)
        
    gsl_vector** **flagTruncated**
    
        Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
        
    gsl_vector** **maxDERgsl** 
    
        Maximum of the first derivative of the (low-pass filtered) record inside each found pulse
        
    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector
        
    .. cpp:member:: gsl_vector* der
    
        First derivative of the (low-pass filtered) record
        
    .. cpp:member:: double adaptativethreshold
    
        Threshold
        
    .. cpp:member:: int nSamplesUp
    
        Number of consecutive samples over the threshold to 'find' a pulse
        
    .. cpp:member:: int* numberPulses
    
        Number of found pulses
        
    .. cpp:member:: gsl_vector** tstartgsl 
        
        Pulses tstart (samples)
        
    .. cpp:member:: gsl_vector** flagTruncated
    
        Flag indicating if the pulse is truncated (inside this function only initial truncated pulses are classified)
        
    .. cpp:member:: gsl_vector** maxDERgsl 
    
        Maximum of the first derivative of the (low-pass filtered) record inside each found pulse

 
.. cpp:function:: int find_Esboundary(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, int *indexEalpha, int *indexEbeta, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*.
    
    This function provides the indexes of the two energies which straddle the pulse energy, by  comparing the maximum value of the pulse derivative
    (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`).

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library:
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` :cpp:member:`indexEalpha` = :cpp:member:`indexEbeta` = 0
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` :cpp:member:`indexEalpha` = :cpp:member:`indexEbeta` = Number of templates-1
        
    **Members/Variables**

    double **maxDER**
    
        Max value of the derivative of the (filtered) pulse whose embracing energies are being sought
        
    gsl_vector* **maxDERs**
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info in the library about the energies
    
    int* **indexEalpha** 
    
        Index of the energy lower than the energy of the pulse which is being analyzed
        
    int* **indexEbeta** 
    
        Index of the energy higher than the energy of the pulse which is being analyzed
        
    double* **Ealpha** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    double* **Ebeta** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose embracing energies are being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info in the library about the energies
    
    .. cpp:member:: int* indexEalpha 
    
        Index of the energy lower than the energy of the pulse which is being analyzed
        
    .. cpp:member:: int* indexEbeta 
    
        Index of the energy higher than the energy of the pulse which is being analyzed
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
        
.. cpp:function:: int find_matchedfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function selects the proper matched filter (normalized template) from the calibration library from column **DAB** (or from column **MF** if only one energy included in                                                the library) by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`) for the *DAB* interpolation method (see :ref:`optimal filter chapter <optimalFilter_NSD>`). It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library:
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**
    
    int **runF0orB0val**

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    double **maxDER**
    
        Max value of the derivative of the (filtered) pulse whose matched filter is being sought
        
    gsl_vector* **maxDERs**
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    ReconstructInitSIRENA* **reconstruct_init** 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    gsl_vector** **matchedfilterFound**
    
        GSL vector with the matched filter selected

    gsl_vector** **PabFound** 
    
        **PAB** column from the library 
        
    double* **Ealpha** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    double* **Ebeta**
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
    .. cpp:member:: int runF0orB0val

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose matched filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: gsl_vector** matchedfilterFound
    
        GSL vector with the matched filter selected

    .. cpp:member:: gsl_vector** PabFound 
    
        **PAB** column from the library 
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
        
.. cpp:function:: int find_model_energies(double energy, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)

    Located in file: *pulseprocess.cpp*
    
    This function uses :cpp:member:`energy` in order to choose the proper pulse template (*pulse_templates_B0*) of the calibration library.

    In general, it finds the two energies wich straddle :cpp:member:`energy` in the calibration library and interpolates (:cpp:func:`interpolate_model`):
    
      - If :cpp:member:`energy` is lower than the lowest energy in the library :math:`\Rightarrow` The model with the lowest energy in the library is chosen
      - If :cpp:member:`energy` is higher than the highest energy in the library :math:`\Rightarrow` The model with the highest energy in the library is chosen

    **Members/Variables**
    
    double **energy**
    
        Energy of the pulse whose pulse template is being sought
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the energies of the models (*energies*)  
        and their templates (*pulse_templates*), the number of templates in the library (*ntemplates*), the template duration (*template_duration*) and
        the *pulse_templates_B0*.

    gsl_vector** **modelFound**
        
        Found template of the pulse whose energy is :cpp:member:`energy`
        
    .. cpp:member:: double energy
    
        Energy of the pulse whose pulse template is being sought
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the energies of the models (*energies*)  
        and their templates (*pulse_templates*), the number of templates in the library (*ntemplates*), the template duration (*template_duration*) and
        the *pulse_templates_B0*.

    .. cpp:member:: gsl_vector** modelFound
        
        Found template of the pulse whose energy is :cpp:member:`energy`
        
        
.. cpp:function:: int find_model_maxDERs(double maxDER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)

    Located in file: *pulseprocess.cpp*
    
    This function uses the maximum of the derivative of the (filtered) pulse (:cpp:member:`maxDER`) in order to choose the proper
    pulse template (*pulse_templates_filder*) of the calibration library.

    In general, it finds the two *maxDER* which straddle :cpp:member:`maxDER` in the calibration library and interpolates (:cpp:func:`interpolate_model`):
    
      - If :cpp:member:`maxDER` is lower than the lowest *maxDERs* in the library :math:`\Rightarrow` The model with
        the lowest *maxDERs* in the library is chosen
      - If :cpp:member:`maxDER` is higher than the highest *maxDERs* in the library :math:`\Rightarrow` The model with
        the highest *maxDERs* in the library is chosen

    **Members/Variables**
    
    double **maxDER**
    
        Maximum of the derivative of the (filtered) pulse whose pulse template is being sought
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the number of templates in 
        the library (*ntemplates*), the template duration (*template_duration*), the filtered and differentiated templates (*pulse_templates_filder*)
        and the *maxDERs* of the templates

    gsl_vector** **modelFound**
        
        Found template of the pulse whose maximum of the derivative of the filtered version is :cpp:member:`maxDER`
        
    .. cpp:member:: double maxDER
    
        Maximum of the derivative of the (filtered) pulse whose pulse template is being sought
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the number of templates in 
        the library (*ntemplates*), the template duration (*template_duration*), the filtered and differentiated templates (*pulse_templates_filder*)
        and the *maxDERs* of the templates

    .. cpp:member:: gsl_vector** modelFound
        
        Found template of the pulse whose maximum of the derivative of the filtered version is :cpp:member:`maxDER`
        
        
.. cpp:function:: int find_model_samp1DERs(double samp1DER, ReconstructInitSIRENA *reconstruct_init, gsl_vector **modelFound)

    Located in file: *pulseprocess.cpp*
    
    This function uses the 1st sample of the derivative of the filtered pulse (:cpp:member:`samp1DER`) in order to choose the proper pulse template (*pulse_templates_filder*) of the calibration library.

    It finds the two :cpp:member:`samp1DER` closer in the calibration library and interpolates (:cpp:func:`interpolate_model`)
        
      - If :cpp:member:`samp1DER` is lower than the lowest samp1DER in the library :math:`\Rightarrow` The model with the lowest samp1DER in the library is chosen
      - If :cpp:member:`samp1DER` is higher than the highest samp1DER in the library :math:`\Rightarrow` The model with the highest samp1DER in the library is chosen

    **Members/Variables**
    
    double **samp1DER**
    
        1st sample of the derivative of the filtered pulse whose pulse template is being sought
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the 1st samples of the derivative of the models (*samp1DERs*) and their derived templates (*pulse_templates_filder*), the number of templates in the library (*ntemplates*) and the template duration (*template_duration*).

    gsl_vector** **modelFound**
        
        Found template of the pulse whose 1st sample of the derivative of the filtered pulse is :cpp:member:`samp1DER`
        
    .. cpp:member:: double samp1DER
    
        1st sample of the derivative of the filtered pulse whose pulse template is being sought
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the 1st samples of the derivative of the models (*samp1DERs*) and their derived templates (*pulse_templates_filder*), the number of templates in the library (*ntemplates*) and the template duration (*template_duration*).

    .. cpp:member:: gsl_vector** modelFound
        
        Found template of the pulse whose 1st sample of the derivative of the filtered pulse is :cpp:member:`samp1DER`
        
                
.. cpp:function:: int find_optimalfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    (or 'Tx'or 'Fx'columns if only one energy included in                                                *                        the library) 
    This function selects the proper optimal filter from the calibration library columns **ABTx** or **ABFx** (or from **Tx** or **Fx**columns if only one energy included in                              the library) by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`). It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library:
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**

    double **maxDER**
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    gsl_vector* **maxDERs**
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    ReconstructInitSIRENA* **reconstruct_init** 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info
        in the library (*optimal_filters*)
    
    gsl_vector** **optimalfilterFound**
    
        GSL vector with the optimal filter selected
        
    gsl_vector** **PabFound** 
    
        **PAB** column from the library
        
    double* **Ealpha** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    double* **Ebeta**
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info
        in the library (*optimal_filters*)
    
    .. cpp:member:: gsl_vector** optimalfilterFound
    
        GSL vector with the optimal filter selected
        
    .. cpp:member:: gsl_vector** PabFound 
    
        **PAB** column from the library
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
       
       
.. cpp:function:: int find_prclofwm(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **PRCLOFWMFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    When :option:`EnergyMethod` = **OPTFILT** and :option:`OFNoise` = **WEIGHTM** this function selects the proper precalculated values (**OFWx**) from the calibration *PRCLOFWM* HDU of the library by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library (:cpp:member:`maxDERs`) for the :option:`OFLib` =yes.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library:
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**

    double **maxDER**
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    gsl_vector* **maxDERs**
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    gsl_vector** **PRCLOFWMFound**
    
        GSL vector with some precalculated selected
        
    double* **Ealpha** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    double* **Ebeta**
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: gsl_vector** PRCLOFWMFound
    
        GSL vector with some precalculated selected
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
	
	
.. cpp:function:: int find_prclwn(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **PRCLWNFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    When :option:`EnergyMethod` = **WEIGHTN** this function selects the proper precalculated values (**PCLx**) from the *PRECALWN* HDU of the  calibration library by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library (:cpp:member:`maxDERs`) for the :option:`OFLib` =yes. It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library:
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**

    double **maxDER**
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    gsl_vector* **maxDERs**
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    ReconstructInitSIRENA* **reconstruct_init** 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    gsl_vector** **PRCLWNFound**
    
        GSL vector with the precalculated values selected
        
    gsl_vector** **PabFound** 
    
        **PAB** column from the library
        
    double* **Ealpha** 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    double* **Ebeta**
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
    
    .. cpp:member:: gsl_vector** PRCLWNFound
    
        GSL vector with the precalculated values selected
        
    .. cpp:member:: gsl_vector** PabFound 
    
        **PAB** column from the library
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
        
.. cpp:function:: extern_C_void freeOptimalFilterSIRENA(OptimalFilterSIRENA* OFilterColl)
    
    Located in file: *integraSIRENA.cpp*
    
    Destructor of *OptimalFilterSIRENA* structure.
    
    **Members/Variables**
    
    OptimalFilterSIRENA* **OFilterColl**
    
        Instance of *OptimalFilterSIRENA* structure
        
    .. cpp:member:: OptimalFilterSIRENA* OFilterColl
    
        Instance of *OptimalFilterSIRENA* structure
    
    
.. cpp:function:: extern_C_void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init)    
    
    Located in file: *integraSIRENA.cpp*
    
    Destructor of *ReconstructInitSIRENA* structure.
    
    **Members/Variables**
    
    ReconstructInitSIRENA* **reconstruct_init**
    
        Instance of *ReconstructInitSIRENA* structure
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Instance of *ReconstructInitSIRENA* structure
    
    
.. cpp:function:: extern_C_void freePulsesCollection(PulsesCollection* PulsesColl)
    
    Located in file: *integraSIRENA.cpp*
    
    Destructor of *PulsesCollection* structure.
    
    **Members/Variables**
    
    PulsesCollection* **PulsesColl**
        
        Instance of *PulsesCollection* structure

    .. cpp:member:: PulsesCollection* PulsesColl
        
        Instance of *PulsesCollection* structure
        
    
.. cpp:function:: int fromGslMatrix(void **buffer, gsl_matrix **matrix, int type)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input GSL matrix into an output buffer.
    
    **Members/Variables**
    
    void** **buffer**
    
        Output buffer
    
    gsl_matrix** **matrix**
    
        Input GSL matrix
        
    int **type**

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
        
    .. cpp:member:: void** buffer
    
        Output buffer
    
    .. cpp:member:: gsl_matrix** matrix
    
        Input GSL matrix
        
    .. cpp:member:: int type

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
    
.. cpp:function:: int fromGslVector(void **buffer, gsl_vector **array, int type)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input GSL vector into an output buffer.
    
    **Members/Variables**
    
    void** **buffer**
    
        Output buffer
        
    gsl_vector** **array**
    
        Input GSL vector
        
    int **type**

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
        
    .. cpp:member:: void** buffer
    
        Output buffer
        
    .. cpp:member:: gsl_vector** array
    
        Input GSL vector
        
    .. cpp:member:: int type

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
    
.. _G:

.. cpp:function:: int gennoisespec_main()

    Located in file: *gennoisespec.cpp*
    
    This function calculates the current noise spectral density.
    If there are pulses in a record, the pulses are rejected and it is going to look for pulse-free intervals of a given size (*intervalMinBins*).
    If there are no pulses in a record, the event is divided into pulse-free intervals of a given size (*intervalMinBins*).
    It is going to look for pulse-free intervals, calculate their FFT(not filtered data) and average them.
    
    Another facillity is calculate the weight matrix of the noise (in fact, weight matrixes of the noise of different lengths).
    
    The output FITS file (_noisespec) contains three columns in two extensions, *NOISE* and *NOISEALL*:
        - **FREQ**: Frequency
        - **CSD**: Current noise spectral density: Amount of current per unit (density) of frequency (spectral), as a function of the frequency
        - **SIGMACSD**: Standard error of the mean (filled out with 0's at the moment)
        
    There is also other extension, *WEIGHTMS*, where the weight matrices of the noise are stored.
    
    Steps: 
       
    - Reading all programm parameters by using PIL
    - Open input FITS file
    - Check if input FITS file have been simulated with TESSIM or XIFUSIM
    - To calculate *aducnv* (conversion factor between arbitrary units and A)...
    - ...or read ``ADU_CNV``, ``I_BIAS`` and ``ADU_BIAS``
    - Get structure of input FITS file columns
    - Read info to transform to resistance space
    - Read and check other input keywords
    - Read other necessary keywords from ANY HDU
    - Calculate the sampling rate
        - By using keywords in input FITS file (from ``DELTAT`` or ``TCLOCK``+``DEC_FAC`` or ``NUMROW``+``P_ROW``)
        - If necessary read the sampling rate from input FITS file (from the ``HISTORY`` in the *Primary* HDU)
        - If not possible, provide an error message to include DELTAT (inverse of sampling rate) in the input FITS file
    - Initialize variables and transform from seconds to samples
    - Declare variables
    - Create structure to run Iteration
    - Read columns (**TIME** and **ADC**)
        - Called iteration function: :cpp:func:`inDataIterator`
    - Close input FITS file
    - Generate CSD representation
        - Applying :cpp:func:`medianKappaClipping_noiseSigma` in order to remove the noise intervals with a high sigma (if :option:`rmNoiseInterval` = yes)
        - FFT calculus (EventSamplesFFT)
        - Add to mean FFT samples
        - Current noise spectral density
        - Extra normalization (further than the FFT normalization factor,1/n) in order to get the apropriate noise level provided by Peille (54 pA/rHz)
    - Load in noiseIntervals only those intervals with a proper sigma (if :option:`rmNoiseInterval` = yes) and NumMeanSamples = cnt (in order not to change excesively the code when weightMS) 
    - Generate WEIGHT representation
    - Create output FITS File: GENNOISESPEC representation file 
    - Write extensions *NOISE*, *NOISEALL* and *WEIGHTMS* (call :cpp:func:`writeTPSreprExten`)
    - Free allocated GSL vectors
    - Close output FITS file
    - Free memory
    - Finalize the task
    
    The parameters (*struct Parameters* **par**) read by :cpp:func:`getpar_noiseSpec` are:
    
    char **inFile** 
    
        Name of the input FITS file
        
    char **outFile** 
    
        Name of the output FITS file
    
    int **intervalMinSamples**
    
        Minimum length of a pulse-free interval to use (samples) = *intervalMinBins*
    
    int **nplPF**
    
        Number of pulse lengths after ending the pulse (Tend) to start the pulse-free interval
    
    int **nintervals**
    
        Number of pulse-free intervals to use to calculate the Noise Spectral Density
        
    double **scaleFactor**
    
        Scale factor to apply in order to calculate the LPF box-car length
        
    int **samplesUp**
    
        Consecutive samples over the threshold to locate a pulse
        
    double **nSgms**
    
        Number of Sigmas to establish the threshold
        
    int **pulse_length**
    
        Pulse length (samples)
        
    double **LrsT**
    
        Running sum length (seconds) 
        
    double **LbT**
    
        Baseline averaging length (seconds)
        
    char **weightMS**
    
        Calculate and write the weight matrixes if *weightMS=yes*
        
    char **I2R** 
    
        Transform to resistance space (I2R, I2RALL, I2RNOL, I2RFITTED) or not (I)
        
    char **clobber**
    
        Re-write output files if clobber=yes
        
    int **matrixSize**
    
        Size of noise matrix if only one to be created
        
    char **rmNoiseIntervals**
    
        Remove some noise intervals before calculating the noise spectrum if *rmNoiseIntervals=yes*

.. cpp:function:: int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B, gsl_vector **rmsB)
    
    Located in file: *pulseprocess.cpp*
    
    This function calculates the sum, :cpp:member:`B`, of :cpp:member:`lb` digitized data samples of a pulse-free interval immediately
    before each pulse. If the pulse-free interval before the current pulse is lower than :cpp:member:`lb`, :cpp:member:`B` is calculated with the available
    number of samples. If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse. 
    The number of samples of the pulse-free interval used to calculate :cpp:member:`B` is stored in the :cpp:member:`lb` vector.
    
    Steps: 
    
    First of all, the auxiliary variable *Baux* is initialized to -999 and all the elements of the :cpp:member:`lb` vector are equal to the :option:`LbT` input parameter in samples.
    Then, the code is divided into 2 *if* statements:
    
    - When the current pulse is the first pulse into the record:
    
      - :math:`tstart \geq lb` :math:`\Rightarrow` Sum *lb* samples
      - :math:`0<tstart<lb` :math:`\Rightarrow` Sum the available number of samples (although the available number of samples was lower than *lb*)
      - :math:`tstart=0` :math:`\Rightarrow` If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse
    
    - When the current pulse is not the first pulse into the record:
    
      - :math:`tstart_i-tend_{i-1} \geq lb` :math:`\Rightarrow` Sum lb samples
      - :math:`0<tstart_i-tend{i-1}<lb` :math:`\Rightarrow` Sum the available number of samples (although the available number of samples was lower than *lb*)
      - If there is not a pulse-free interval before the pulse, it is looked for it after the current pulse

    If *Baux* is still -999, a pulse-free interval can not be found to apply the running sum filter. This has to be taken into account,
    out of the function, to try to get a usable :cpp:member:`B`.
    
    **Members/Variables**
    
    gsl_vector* **vectorin**
         
        Input record
        
    gsl_vector* **tstart**
         
        Starting time of the pulses into the record
        
    int **nPulses**
         
        Number of pulses into the record
        
    gsl_vector** **lb**
         
        Vector containing the baseline averaging length used for each pulse
        
    int **sizepulse**
         
        Size of the pulse in samples
        
    gsl_vector** **B**
         
        In general, sum of the *Lb* digitized data samples (:option:`LbT` input parameters in samples) of a pulse-free interval immediately before the current pulse

    gsl_vector** **rmsB**
         
        In general, rms of the baseline related to a pulse-free interval immediately before the current pulse
        
    .. cpp:member:: gsl_vector* vectorin
         
        Input record
        
    .. cpp:member:: gsl_vector* tstart
         
        Starting time of the pulses into the record
        
    .. cpp:member:: int nPulses
         
        Number of pulses into the record
        
    .. cpp:member:: gsl_vector** lb
         
        Vector containing the baseline averaging length used for each pulse
        
    .. cpp:member:: int sizepulse
         
        Size of the pulse in samples
        
    .. cpp:member:: gsl_vector** B
         
        In general, sum of the *Lb* digitized data samples (:option:`LbT` input parameters in samples) of a pulse-free interval immediately before the current pulse

    .. cpp:member:: gsl_vector** rmsB
         
        In general, rms of the baseline related to a pulse-free interval immediately before the current pulse
    
.. cpp:function:: LibraryCollection* getLibraryCollection(const char* const filename, int opmode, int hduPRECALWN, int hduPRCLOFWM, int largeFilter, char* filter_domain, int pulse_length, char *energy_method, char *ofnoise, char *filter_method, char oflib, char **ofinterp, double filtEev, int lagsornot, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function creates and retrieves a *LibraryCollection* from a file.
    
    - Create *LibraryCollection* structure
    - Open FITS file in READONLY mode (move to the first HDU) and get number of templates (rows)
    - Allocate library structure
    - Get **PULSE** and **MF** column numbers (depending the different options)
    - Get template duration
    - Allocate library structure (cont.)
    - Get matched filter duration
    - Read different columns and populate the *LibraryCollection* structure
    - Added new code to handle the new HDUs *FIXFILTF*, *FIXFILTT*, *PRECALWN* and *PRCLOFWM*
    - Free allocated GSL vectors and matrices
    
    **Members/Variables**
    
    const char* const **filename**
        
        File with library information
        
    int **opmode**
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
    
     int **hduPRECALWN**
    
        Add or not the *PRECALWN* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRECALWN`
        
    int **hduPRCLOFWM**
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    int **largeFilter**
    
        Length of the longest fixed filters (only for library creation, :option:`opmode` = 0), :option:`largeFilter`

    char* **filter_domain**
    
        Filtering Domain: Time (**T**) or Frequency (**F**), :option:`FilterDomain`
        
    int **pulse_length**
    
        Pulse length, :option:`PulseLength`
        
    char* **energy_method**
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`

    char* **ofnoise**
    
        For optimal filtering : **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    char* **filter_method**
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
    
    char **oflib**
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
        
    char** **ofinterp**
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in :ref:`tesreconstruction` as *DAB* (but it would be possible to work with *MF*)
        
    double **filtEev**  
    
        Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RALL, I2RNOL and I2RFITTED), :option:`filtEeV`
        
    int* const **status**
    
        Input/output status

    .. cpp:member:: const char* const filename
        
        File with library information
        
    .. cpp:member:: int opmode
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
    
     .. cpp:member:: int hduPRECALWN
    
        Add or not the *PRECALWN* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRECALWN`
        
    .. cpp:member:: int hduPRCLOFWM
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    .. cpp:member:: int largeFilter
    
        Length of the longest fixed filters (only for library creation, :option:`opmode` = 0), :option:`largeFilter`

    .. cpp:member:: char* filter_domain
    
        Filtering Domain: Time (**T**) or Frequency (**F**), :option:`FilterDomain`
        
    .. cpp:member:: int pulse_length
    
        Pulse length, :option:`PulseLength`
        
    .. cpp:member:: char* energy_method
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`

    .. cpp:member:: char* ofnoise
    
        For optimal filtering : **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
    
    .. cpp:member:: char oflib
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
        
    .. cpp:member:: char** ofinterp
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in :ref:`tesreconstruction` as *DAB* (but it would be possible to work with *MF*)
        
    .. cpp:member:: double filtEev  
    
        Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RALL, I2RNOL and I2RFITTED), :option:`filtEeV`
        
    .. cpp:member:: int* const status
    
        Input/output status
        
    
.. cpp:function:: NoiseSpec* getNoiseSpec(const char* const filename, int opmode, int hduPRCLOFWM, char *energy_method, char *ofnoise, char *filter_method, int* const status)
    
    Located  in file: *integraSIRENA.cpp*
    
    This function creates and retrieves a *NoiseSpec* from a file.
    
    - Create *NoiseSpec* structure
    - Open FITS file, move to the *NOISE*, *NOISEALL* and *WEIGHTMS* HDUs and get necessary keywords
    - Allocate *NoiseSpec* structure
    - Get noise spectrum (**CSD**), and noise frequencies (**FREQ**) column numbers
    - Read column **CSD** and save it into the structure
    - Read column **FREQ** and save it into the structure
    - Read columns **Wx** with the noise weight matrix from noise intervals and save them into the structure
    - Return noise spectrum
    
    **Members/Variables**
    
    const char* const **filename**
    
        File name with noise
        
    int **opmode**
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
        
    int **hduPRCLOFWM**
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    char* **energy_method**
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`
        
    char* **ofnoise**
    
         For optimal filtering:  **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    char* **filter_method**
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
        
    int* const **status**
    
        Input/Output status
        
    .. cpp:member:: const char* const filename
    
        File name with noise
        
    .. cpp:member:: int opmode
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
        
    .. cpp:member:: int hduPRCLOFWM
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    .. cpp:member:: char* energy_method
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`
        
    .. cpp:member:: char* ofnoise
    
         For optimal filtering:  **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
        
    .. cpp:member:: int* const status
    
        Input/Output status
    
    
.. cpp:function:: int getpar(struct Parameters* const par)

    Located in file: *tesreconstruction.cpp*
    
    This function gets the input parameter from the command line or their default values from the *tesreconstruction.par* file
    
    **Members/Variables**
    
    struct Parameters* const **par**
    
        Structure containing the input parameters specified in *tesreconstruction.par*
        
    .. cpp:member:: struct Parameters* const par
    
        Structure containing the input parameters specified in *tesreconstruction.par*
        
        
.. cpp:function:: int getpar_noiseSpec(struct Parameters* const par)

    Located in file: *gennoisespec.cpp*
    
    This function gets the input parameter from the command line or their default values from the *gennoisespec.par* file
    
    **Members/Variables**
    
    struct Parameters* const **par**
    
        Structure containing the input parameters specified in *gennoisespec.par*
        
    .. cpp:member:: struct Parameters* const par
    
        Structure containing the input parameters specified in *gennoisespec.par*
   
.. cpp:function:: int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight)
    
    Located in file: *pulseprocess.cpp*
    
    This function estimates the pulse height of a pulse by using a running sum filter. It extracts from the record, :cpp:member:`vectorin`, the pulse whose 
    pulse height is going to be estimated by using :cpp:func:`RS_filter`.

    Steps:
    
    - Declare variables
    - Extracting from the record the pulse whose pulse height is going to be estimated
    - Apply the running sum filter
    
    **Members/Variables**
    
    gsl_vector* **vectorin**
    
        Not filtered record
        
    double **tstart**
    
        Starting time of the pulse whose pulse height is going to be estimated
        
    double **tstartnext**
    
        Starting time of the next pulse whose pulse height is going to be estimated
        
    int **lastPulse**
    
        It is 1 if the pulse is the last one into the record or the only one
        
    double **lrs**
    
        Running sum length (equal to the :option:`LrsT` input parameter in samples)

    double **lb**
    
        Baseline averaging length used for the pulse whose pulse height is going to be estimated
        
    double **B**
    
        In general, sum of the *Lb* digitized data samples (:option:`LbT` input parameters in samples) of a pulse-free interval immediately before the current pulse
        
    int **sizepulse**
    
        Size of the pulse in samples
        
    double* **pulseheight**
    
        Estimated pulse height of the pulse

    .. cpp:member:: gsl_vector* vectorin
    
        Not filtered record
        
    .. cpp:member:: double tstart
    
        Starting time of the pulse whose pulse height is going to be estimated
        
    .. cpp:member:: double tstartnext
    
        Starting time of the next pulse whose pulse height is going to be estimated
        
    .. cpp:member:: int lastPulse
    
        It is 1 if the pulse is the last one into the record or the only one
        
    .. cpp:member:: double lrs
    
        Running sum length (equal to the :option:`LrsT` input parameter in samples)

    .. cpp:member:: double lb
    
        Baseline averaging length used for the pulse whose pulse height is going to be estimated
        
    .. cpp:member:: double B
    
        In general, sum of the *Lb* digitized data samples (:option:`LbT` input parameters in samples) of a pulse-free interval immediately before the current pulse
        
    .. cpp:member:: int sizepulse
    
        Size of the pulse in samples
        
    .. cpp:member:: double* pulseheight
    
        Estimated pulse height of the pulse
            
            
.. cpp:function:: void gsl_vector_complex_absIFCA(gsl_vector *cvnew,gsl_vector_complex *cv)
    
    Located in file: *genutils.cpp*
    
    This function calculates the magnitude of the complex elements of a vector (real part).
    
    **Members/Variables**
    
    gsl_vector_complex* **cv**
    
        Input GSL complex vector
        
    gsl_vector* **cvnew**
    
        Output GSL vector with the absolute values of the elements of :cpp:member:`cv`

    .. cpp:member:: gsl_vector_complex* cv
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector* cvnew
    
        Output GSL vector with the absolute values of the elements of :cpp:member:`cv`
    
    
.. cpp:function:: void gsl_vector_complex_argIFCA(gsl_vector *varg, gsl_vector_complex *vin)
    
    Located in file: *genutils.cpp*

    This function calculates the arguments of the complex elements of a vector.
    
    **Members/Variables**
    
    gsl_vector_complex* **vin**
    
        Input GSL complex vector
        
    gsl_vector* **varg**
    
        Output GSL vector with the arguments of the elements of :cpp:member:`vin`
        
    .. cpp:member:: gsl_vector_complex* vin
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector* varg
    
        Output GSL vector with the arguments of the elements of :cpp:member:`vin`
    
    
.. cpp:function:: void gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z)
    
    Located in file: *genutils.cpp*
    
    This function multiplies the complex elements of a vector by a complex number.
    
    **Members/Variables**
    
    gsl_vector_complex* **cv**
    
        Input/Output (scaled) GSL complex vector
        
    gsl_complex **z**
    
        Input GSL complex number

    .. cpp:member:: gsl_vector_complex* cv
    
        Input/Output (scaled) GSL complex vector
        
    .. cpp:member:: gsl_complex z
    
        Input GSL complex number
        
    
.. cpp:function:: void gsl_vector_sqrtIFCA(gsl_vector *cvnew, gsl_vector *cv)
    
    Located in file: *genutils.cpp*
    
    This function calculates the square root of the elements of a vector.
    
    **Members/Variables**

    gsl_vector* **cv**
    
        Input GSL complex vector
        
    gsl_vector* **cvnew**
    
        Output GSL vector with the square root values of the elements of :cpp:member:`cv`
        
    .. cpp:member:: gsl_vector* cv
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector* cvnew
    
        Output GSL vector with the square root values of the elements of :cpp:member:`cv`
    
    
.. cpp:function:: int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum)
    
    Located in file: *genutils.cpp*
    
    This function returns the sum of some elements of the input vector.
    
    The starting element of the sum is :cpp:member:`offset` from the start of the input vector. It will sum up :cpp:member:`n` elements.

    :cpp:member:`offset` can take values from 0 to *invector->size*

    **Members/Variables**
    
    gsl_vector* **invector**
    
        Input GSL vector
        
    long **offset**
    
        It is the first element to be summed
        
    long **n**
        
        Number of elements in the sum

    double* **sum**
    
        Calculated output value (sum of the corresponding elements)   
        
    .. cpp:member:: gsl_vector* invector
    
        Input GSL vector
        
    .. cpp:member:: long offset
    
        It is the first element to be summed
        
    .. cpp:member:: long n
        
        Number of elements in the sum

    .. cpp:member:: double* sum
    
        Calculated output value (sum of the corresponding elements)    
        
            
.. _H:

.. _I:

.. cpp:function:: int inDataIterator(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)

    Located in file: *gennoisespec.cpp*
    
    This function takes the optimum number of rows to read the input FITS file and works iteratively

    Steps:
    
    - Declare variables
    - Allocate input GSL vectors
    - Read iterator
    - Processing each record
        - Information has been read by blocks (with nrows per block)
        - Just in case the last record has been filled out with 0's :math:`\Rightarrow` Last record discarded
        - Convert to the resistance space if necessary
        - To avoid taking into account the pulse tails at the beginning of a record as part of a pulse-free interval
        - Low-pass filtering
   	- Differentiate 
   	- Finding the pulses: Pulses tstarts are found (call :cpp:func:`findPulsesNoise`)
        - Finding the pulse-free intervals in each record
            - If there are pulses :math:`\Rightarrow` Call :cpp:func:`findInterval`
            - No pulses :math:`\Rightarrow` The whole event is going to be used (DIVIDING into intervals of intervalMinBins size) :math:`\Rightarrow` Call :cpp:func:`findIntervalN`
        - Calculating the mean and sigma of the intervals without pulses together :math:`\Rightarrow` *BSLN0* and *NOISESTD*
    - Preparing the CSD calculus (not filtered data)
    - Free allocated GSL vectors
    
    **Members/Variables**
    
    long **totalrows**
    
        Total number of rows processed
    
    long **offset**
    
        If positive, this number of rows at the beginning of the table (or pixels in the image) will be skipped and will not be passed to the work function
    
    long **firstrow**
    
        First row to read
    
    long **nrows**
    
        It specifies the number of table rows that are to be passed to the work function on each iteration. If *nrows = 0* then CFITSIO will calculate the optimum number for greatest efficiency. If *nrows* is negative, then all the rows or pixels will be passed at one time, and the work function will only be called once. If any variable length arrays are being processed, then the *nrows* value is ignored, and the iterator will always process one row of the table at a time
    
    int **ncols**
    
        Number of columns
    
    iteratorCol* **cols**
    
        Structure of iteration
    
    void* **user_strct**
    
        This is a user supplied pointer that can be used to pass ancillary information from the driver routine to the work function. It may point to a single number, an array, or to a structure containing an arbitrary set of parameters
        
    .. cpp:member:: long totalrows
    
        Total number of rows processed
    
    .. cpp:member:: long offset
    
        If positive, this number of rows at the beginning of the table (or pixels in the image) will be skipped and will not be passed to the work function
    
    .. cpp:member:: long firstrow
    
        First row to read
    
    .. cpp:member:: long nrows
    
        It specifies the number of table rows that are to be passed to the work function on each iteration. If *nrows = 0* then CFITSIO will calculate the optimum number for greatest efficiency. If *nrows* is negative, then all the rows or pixels will be passed at one time, and the work function will only be called once. If any variable length arrays are being processed, then the *nrows* value is ignored, and the iterator will always process one row of the table at a time
    
    .. cpp:member:: int ncols
    
        Number of columns
    
    .. cpp:member:: iteratorCol* cols
    
        Structure of iteration
    
    .. cpp:member:: void* user_strct
    
        This is a user supplied pointer that can be used to pass ancillary information from the driver routine to the work function. It may point to a single number, an array, or to a structure containing an arbitrary set of parameters

.. cpp:function:: extern_C_void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file, fitsfile *fptr, char* const library_file, char* const event_file, int pulse_length, double scaleFactor, int samplesUp, int samplesDown, double nSgms, int detectSP, int opmode, char *detectionMode, double LrsT, double LbT, char* const noise_file, char* filter_domain, char* filter_method, char* energy_method, double filtEev, char *ofnoise, int lagsornot, int nLags, int Fitting35, int ofiter, char oflib, char *ofinterp, char* oflength_strategy, int oflength, int preBuffer, double monoenergy, char hduPRECALWN, char hduPRCLOFWM, int largeFilter, int interm, char* const detectFile, int errorT, int Sum0Filt, char clobber, int maxPulsesPerRecord, double SaturationValue, char* const tstartPulse1, int tstartPulse2, int tstartPulse3, double energyPCA1, double energyPCA2, char * const XMLFile, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function initializes the structure *ReconstructInitSIRENA* with the variables required for SIRENA reconstruction. The values are taken from the input parameters.
    
    - Load *LibraryCollection* structure if library file exists
    - Load *NoiseSpec* structure
    - Fill in the matrix *tstartPulse1_i* if *tstartPulse1* = nameFile Start time (in samples) of the first pulse (0  if detection should be performed by the system; greater than 0 if provided by the user) or file name containing the tstart (in seconds) of every pulse, :option:`tstartPulse1`
    - Fill in *reconstruct_init*
    
    **Members/Variables**
    
    ReconstructInitSIRENA* **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    char* const **record_file**

        Filename of input data file with records, :option:`RecordFile`
        
    fitsfile* **fptr**

        FITS object with pointer to data file
        
    char* const **library_file**

        File name of calibration library, :option:`LibraryFile`
        
    char* const **event_file**
    
        File name of output events (with reconstructed energy), :option:`TesEventFile`
        
    int **pulse_length**
    
        Pulse length, :option:`PulseLength`
        
    double **scaleFactor**
    
        Detection scale factor for initial filtering, :option:`scaleFactor`
        
    int **samplesUp**
    
        Number of samples for threshold trespassing, :option:`samplesUp`
        
    int **samplesDown**
    
        Number of samples below the threshold to look for other pulse, :option:`samplesDown`
        
    double **nSgms**
    
        Number of standard deviations in the kappa-clipping process for threshold estimation, :option:`nSgms`
        
    int **detectSP**
    
        Detect secondary pulses (1) or not (0), :option:`detectSP`
        
    int **opmode**
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
        
    char* **detectionMode**
    
        Adjusted Derivative (AD) or Single Threshold Crossing (STC), :option:`detectionMode`
    
    double **LrsT**
    
        Running sum length for the RS raw energy estimation (seconds), :option:`LrsT`
    
    double **LbT**
    
        Baseline averaging length (seconds), :option:`LbT`
    
    char* const **noise_file**
    
        Noise file, :option:`NoiseFile`
    
    char* **filter_domain**
    
        Filtering Domain: Time (**T**) or Frequency (**F**), :option:`FilterDomain`
    
    char* **filter_method**
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **F0** (deleting the baseline), :option:`FilterMethod`
    
    char* **energy_method**
    
         Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`
         
    double **filtEev**
    
         Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RALL, I2RNOL and I2RFITTED), :option:`filtEeV`

    char* **ofnoise**
    
         For optimal filtering:  **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    int **lagsornot**
    
        Lags (1) or no lags (0), :option:`LagsOrNot`
        
    int **nLags**
    
        Number of lags (positive odd number)
    
    int **Fitting35**
    
        Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
        
    int **ofiter**
    
        Iterate (1) or not iterate (0), :option:`OFIter`
    
    char **oflib**
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
    
    char* **ofinterp**
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in :ref:`tesreconstruction` as *DAB*
    
    char* **oflength_strategy**
    
        Optimal Filter length Strategy: **FREE**, **BYGRADE** or **FIXED**, :option:`OFStrategy`
    
    int **oflength**
    
        Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**), :option:`OFLength`
        
    int **preBuffer**
    
        Some samples added before the starting time of a pulse
    
    double **monoenergy**
    
        Monochromatic energy of input file in eV (only for library creation, :option:`opmode` = 0), :option:`monoenergy`
        
    int **hduPRECALWN**
    
        Add or not the *PRECALWN* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRECALWN`
        
    int **hduPRCLOFWM**
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    int **largeFilter**
    
        Length of the longest fixed filters (only for library creation, :option:`opmode` = 0), :option:`largeFilter`
    
    int **interm**
    
        Write or not intermediate files (1/0), :option:`intermediate`
    
    char* const **detectFile**
    
        Intermediate detections file (if :option:`intermediate` = 1), :option:`detectFile`
        
    int **errorT**
    
        Additional error (in samples) added to the detected time (logically, it changes the reconstructed energies)
        
    int **Sum0Filt**
    
        0-padding: Subtract the sum of the filter (1) or not (0)
        
    char **clobber**
    
        Overwrite or not output files if exist (1/0), :option:`clobber`
    
    int **maxPulsesPerRecord**
    
        Default size of the event list, :option:`EventListSize`
        
    double **SaturationValue**
    
        Saturation level of the ADC curves
    
    int **tstartPulse1**
    
        Start time (in samples) of the first pulse (0 if detection should be performed by the system; greater than 0 if provided by the user) or file name containing the tstart (in seconds) of every pulse, :option:`tstartPulse1`
    
    int **tstartPulse2**
    
        Tstart (samples) of the second pulse, :option:`tstartPulse2`
    
    int **tstartPulse3**
    
        Tstart (samples) of the third pulse (if 0 :math:`\Rightarrow` PAIRS, if not 0 :math:`\Rightarrow` TRIOS), :option:`tstartPulse3`
        
    double **energyPCA1**
    
        First energy (only for :option:`EnergyMethod` = **PCA**)
        
    double **energyPCA2**
    
        Second energy (only for :option:`EnergyMethod` = **PCA**)
        
    char * const **XMLFile**
    
        File name of the XML input file with instrument definition
    
    int* const **status**
    
        Input/Output status

    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: char* const record_file

        Filename of input data file with records, :option:`RecordFile`
        
    .. cpp:member:: fitsfile* fptr

        FITS object with pointer to data file
        
    .. cpp:member:: char* const library_file

        File name of calibration library, :option:`LibraryFile`
        
    .. cpp:member:: char* const event_file
    
        File name of output events (with reconstructed energy), :option:`TesEventFile`
        
    .. cpp:member:: int pulse_length
    
        Pulse length, :option:`PulseLength`
        
    .. cpp:member:: double scaleFactor
    
        Detection scale factor for initial filtering, :option:`scaleFactor`
        
    .. cpp:member:: int samplesUp
    
        Number of samples for threshold trespassing, :option:`samplesUp`
        
    .. cpp:member:: int samplesDown
    
        Number of samples below the threshold to look for other pulse, :option:`samplesDown`
        
    .. cpp:member:: double nSgms
    
        Number of standard deviations in the kappa-clipping process for threshold estimation, :option:`nSgms`
        
    .. cpp:member:: int detectSP
    
        Detect secondary pulses (1) or not (0), :option:`detectSP`
        
    .. cpp:member:: int opmode
    
        Calibration run (0) or energy reconstruction run (1), :option:`opmode`
        
    .. cpp:member:: char* detectionMode
    
        Adjusted Derivative (AD) or Single Threshold Crossing (STC), :option:`detectionMode`
    
    .. cpp:member:: double LrsT
    
        Running sum length for the RS raw energy estimation (seconds), :option:`LrsT`
    
    .. cpp:member:: double LbT
    
        Baseline averaging length (seconds), :option:`LbT`
    
    .. cpp:member:: char* const noise_file
    
        Noise file, :option:`NoiseFile`
    
    .. cpp:member:: char* filter_domain
    
        Filtering Domain: Time (**T**) or Frequency (**F**), :option:`FilterDomain`
    
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **F0** (deleting the baseline), :option:`FilterMethod`
    
    .. cpp:member:: char* energy_method
    
         Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`
         
    .. cpp:member:: double filtEev
    
         Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RALL, I2RNOL and I2RFITTED), :option:`filtEeV`

    .. cpp:member:: char* ofnoise
    
         For optimal filtering:  **NSD** or **WEIGHTM**, :option:`OFNoise`
        
    .. cpp:member:: int lagsornot
    
        Lags (1) or no lags (0), :option:`LagsOrNot`
        
    .. cpp:member:: int nLags
    
        Number of lags (positive odd number)
    
    .. cpp:member:: int Fitting35
    
        Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
        
    .. cpp:member:: int ofiter
    
        Iterate (1) or not iterate (0), :option:`OFIter`
    
    .. cpp:member:: char oflib
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
    
    .. cpp:member:: char* ofinterp
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in :ref:`tesreconstruction` as *DAB*
    
    .. cpp:member:: char* oflength_strategy
    
        Optimal Filter length Strategy: **FREE**, **BYGRADE** or **FIXED**, :option:`OFStrategy`
    
    .. cpp:member:: int oflength
    
        Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**), :option:`OFLength`
        
    .. cpp:member:: int preBuffer
    
        Some samples added before the starting time of a pulse
    
    .. cpp:member:: double monoenergy
    
        Monochromatic energy of input file in eV (only for library creation, :option:`opmode` = 0), :option:`monoenergy`
        
    .. cpp:member:: int hduPRECALWN
    
        Add or not the *PRECALWN* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRECALWN`
        
    .. cpp:member:: int hduPRCLOFWM
    
        Add or not the *PRCLOFWM* HDU in the library file (1/0) (only for library creation, :option:`opmode` = 0), :option:`hduPRCLOFWM`
        
    .. cpp:member:: int largeFilter
    
        Length of the longest fixed filters (only for library creation, :option:`opmode` = 0), :option:`largeFilter`
    
    .. cpp:member:: int interm
    
        Write or not intermediate files (1/0), :option:`intermediate`
    
    .. cpp:member:: char* const detectFile
    
        Intermediate detections file (if :option:`intermediate` = 1), :option:`detectFile`
        
    .. cpp:member:: int errorT
    
        Additional error (in samples) added to the detected time (logically, it changes the reconstructed energies)
        
    .. cpp:member:: int Sum0Filt
    
        0-padding: Subtract the sum of the filter (1) or not (0)
        
    .. cpp:member:: char clobber
    
        Overwrite or not output files if exist (1/0), :option:`clobber`
    
    .. cpp:member:: int maxPulsesPerRecord
    
        Default size of the event list, :option:`EventListSize`
        
    .. cpp:member:: double SaturationValue
    
        Saturation level of the ADC curves
    
    .. cpp:member:: int tstartPulse1
    
        Start time (in samples) of the first pulse (0 if detection should be performed by the system; greater than 0 if provided by the user) or file name containing the tstart (in seconds) of every pulse, :option:`tstartPulse1`
    
    .. cpp:member:: int tstartPulse2
    
        Tstart (samples) of the second pulse, :option:`tstartPulse2`
    
    .. cpp:member:: int tstartPulse3
    
        Tstart (samples) of the third pulse (if 0 :math:`\Rightarrow` PAIRS, if not 0 :math:`\Rightarrow` TRIOS), :option:`tstartPulse3`
        
    .. cpp:member:: double energyPCA1
    
        First energy (only for :option:`EnergyMethod` = **PCA**)
        
    .. cpp:member:: double energyPCA2
    
        Second energy (only for :option:`EnergyMethod` = **PCA**)
        
    .. cpp:member:: char * const XMLFile
    
        File name of the XML input file with instrument definition
    
    .. cpp:member:: int* const status
    
        Input/Output status    
        
.. cpp:function:: int InitialTriggering(gsl_vector *derivative, double nSgms, double scalefactor, double samplingRate, double stopcriteriamkc, double kappamkc, bool *triggerCondition, int *tstart, int *flagTruncated, double *threshold, int tstartProvided)

    Located in file: *pulseprocess.cpp*

    This function finds the first pulse in the input vector, first derivative of the (low-pass filtered) record.

    Steps:
    
    - Declare variables
    - Stablish the :cpp:member:`threshold`
    - It is necessary to find the tstart of the first pulse... 
    
      Obtain tstart of the first pulse in the derivative if :math:`derivative_i>threshold`
      
    - ... Or to use the tstart provided as input parameter

    **Members/Variables**
    
    gsl_vector* **derivative**
    
        First derivative of the (low-pass filtered) record

    double **nSgms**
        
        Number of *Sigmas* to establish the threshold (:option:`nSgms`)
        
    double **scalefactor**
        
        Scale factor to calculate the LPF box-car length (:option:`scaleFactor`)
        
    double **samplingRate**
        
        Sampling rate
        
    double **stopcriteriamkc**
        
        Used in :cpp:func:`medianKappaClipping` (%)
        
    double **kappamkc**
        
        Used in :cpp:func:`medianKappaClipping`
        
    bool* **triggerCondition**
        
        True :math:`\Rightarrow` The algorithm has found the first event
        
        False :math:`\Rightarrow` The algorithm has not found any event
        
    int* **tstart**
        
        First event tstart (in samples)
        
    int* **flagTruncated**
        
        Flag indicating if the event is truncated 
        
    double* **threshold**
        
        Calculated threshold  (output parameter because it is necessary out of the function)
     
    int **tstartProvided**
        
        Tstart of the first pulse provided as input parameter

    .. cpp:member:: gsl_vector* derivative
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double nSgms
        
        Number of *Sigmas* to establish the threshold (:option:`nSgms`)
        
    .. cpp:member:: double scalefactor
        
        Scale factor to calculate the LPF box-car length (:option:`scaleFactor`)
        
    .. cpp:member:: double samplingRate
        
        Sampling rate
        
    .. cpp:member:: double stopcriteriamkc
        
        Used in :cpp:func:`medianKappaClipping` (%)
        
    .. cpp:member:: double kappamkc
        
        Used in :cpp:func:`medianKappaClipping`
        
    .. cpp:member:: bool* triggerCondition
        
        True :math:`\Rightarrow` The algorithm has found the first event
        
        False :math:`\Rightarrow` The algorithm has not found any event
        
    .. cpp:member:: int* tstart
        
        First event tstart (in samples)
        
    .. cpp:member:: int* flagTruncated
        
        Flag indicating if the event is truncated 
        
    .. cpp:member:: double* threshold
        
        Calculated threshold  (output parameter because it is necessary out of the function)
     
    .. cpp:member:: int tstartProvided
        
        Tstart of the first pulse provided as input parameter    
        
        
.. cpp:function:: int interactivePars(inparam *taskPars, int np, string task)
    
    Located in file *inoututils.cpp*
    
    This function reads input parameters interactively (provided by the user or taken as default values).
    Used in tool :ref:`gennoisespec`.
    
    **Members/Variables**
    
    inparam* **taskPars**
    
        Instance of *inparam* structure storing input parameters
        
    int **np**
    
        Number of parameters
    
    string **task**
    
        Tool name
        
    .. cpp:member:: inparam* taskPars
    
        Instance of *inparam* structure storing input parameters
        
    .. cpp:member:: int np
    
        Number of parameters
    
    .. cpp:member:: string task
    
        Tool name
        
.. cpp:function:: int interpolatePOS(gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function interpolates an input vector (:cpp:member:`x_in`, :cpp:member:`y_in`), creating an output vector (:cpp:member:`x_out`, :cpp:member:`y_out`) with the size and frequency step given. *POS* comes from the fact that the input spectrum only has positive frequencies (in order to not handle the f=0 bin).

    - Declare and initialize variables
    - GSL method applied for interpolatation
    - Generate the interpolated output vector
    - Free memory

    **Members/Variables**
    
    gsl_vector* **x_in**
    
        GSL input vector with the abscissas of the vector which is going to be interpolated 
        
    gsl_vector* **y_in**

        GSL input vector with the ordinates of the vector which is going to be interpolated 
        
    long **size** 
    
        Size of the interpolated output vector
        
    double **step**
    
        Frequency step of the interpolated output vector
    
    gsl_vector** **x_out**
    
        GSL output vector with the abscissas of the interpolated vector
        
    gsl_vector** **y_out**

        GSL output vector with the ordinates of the interpolated vector
        
    .. cpp:member:: gsl_vector* x_in
    
        GSL input vector with the abscissas of the vector which is going to be interpolated 
        
    .. cpp:member:: gsl_vector* y_in

        GSL input vector with the ordinates of the vector which is going to be interpolated 
        
    .. cpp:member:: long size 
    
        Size of the interpolated output vector
        
    .. cpp:member:: double step
    
        Frequency step of the interpolated output vector
    
    .. cpp:member:: gsl_vector** x_out
    
        GSL output vector with the abscissas of the interpolated vector
        
    .. cpp:member:: gsl_vector** y_out

        GSL output vector with the ordinates of the interpolated vector
            
            
.. cpp:function:: int interpolate_model(gsl_vector **modelFound, double p_model, gsl_vector *modelIn1, double p_modelIn1, gsl_vector *modelIn2, double p_modelIn2)

    Located in file: *pulseprocess.cpp*

    This function interpolates the pulse model, :math:`p(t,E)`, between two models of the pulse models library,
    :math:`p(t,E_1)` and :math:`p(t,E_2)`, being :math:`E_1<E<E_2`.

    According to the interpolation method:
    
    .. math::
   
        p(t,E)={\frac{E_2-E}{E_2-E_1}}p(t,E_1)+{\frac{E-E_1}{E_2-E_1}}p(t,E_2)
        
        
    **Members/Variables**
    
    gsl_vector** **modelFound**
    
        Found model of the pulse whose *energy* or *maxDER* is :cpp:member:`p_model`
        
    double **p_model**
    
        Parameter (*energy* or *maxDER*) of the pulse whose model is being sought

    gsl_vector* **modelIn1**
        
        Model of the pulse whose parameter (*energy* or *maxDER*) is immediately lower than :cpp:member:`p_model` in the library FITS file
        
    double **p_modelIn1**
        
        Parameter (*energy* or *maxDER*) immediately lower than :cpp:member:`p_model` in the library FITS file
        
    gsl_vector* **modelIn2**
        
        Model of the pulse whose parameter (*energy* or *maxDER*) is immediately greater than :cpp:member:`p_model` in the library FITS file
        
    double **p_modelIn2**
        
        Parameter (*energy* or *maxDER*) immediately greater than :cpp:member:`p_model` in the library FITS file
        
    .. cpp:member:: gsl_vector** modelFound
    
        Found model of the pulse whose *energy* or *maxDER* is :cpp:member:`p_model`
        
    .. cpp:member:: double p_model
    
        Parameter (*energy* or *maxDER*) of the pulse whose model is being sought

    .. cpp:member:: gsl_vector* modelIn1
        
        Model of the pulse whose parameter (*energy* or *maxDER*) is immediately lower than :cpp:member:`p_model` in the library FITS file
        
    .. cpp:member:: double p_modelIn1
        
        Parameter (*energy* or *maxDER*) immediately lower than :cpp:member:`p_model` in the library FITS file
        
    .. cpp:member:: gsl_vector* modelIn2
        
        Model of the pulse whose parameter (*energy* or *maxDER*) is immediately greater than :cpp:member:`p_model` in the library FITS file
        
    .. cpp:member:: double p_modelIn2
        
        Parameter (*energy* or *maxDER*) immediately greater than :cpp:member:`p_model` in the library FITS file
        
        
.. cpp:function:: bool isNumber(string s)

    Located in file: *genutils.cpp*

    This function returns TRUE if the input string is a number or FALSE if not.
        
    **Members/Variables**
    
    string **s**
        
        Input string 

    .. cpp:member:: string s
        
        Input string 
        
        
.. _J:

.. _K:

.. _L:
    
.. cpp:function:: int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double)
        
    Located in file: *tasksSIRENA.cpp*
    
    This fucntion loads the structure :cpp:member:`record` into the :cpp:member:`adc_double` GSL vector.

    It checks if the record has been filled out with 0's :math:`\Rightarrow` It only loads the first values (which are different from 0).

    **Members/Variables**

    TesRecord* **record** 

        Member of *TesRecord* structure that contains the input record 

    double **time_record** 

        Starting time of the record (output)

    gsl_vector** **adc_double** 

        Storage of the record to be processed (input/output)
        
    .. cpp:member:: TesRecord* record 

        Member of *TesRecord* structure that contains the input record 

    .. cpp:member:: double time_record 

        Starting time of the record (output)

    .. cpp:member:: gsl_vector** adc_double 

        Storage of the record to be processed (input/output)

            
.. cpp:function:: int lpf_boxcar (gsl_vector **invector, int szVct, int sampleRate)
        
    Located in file: *pulseprocess.cpp*
    
    This function implements a low pass filtering as a box-car function in time.

    The box-car function is a temporal average window:

    .. math::

        x_{i-1}=\sum_{0}^{n-1}\frac{I_i}{n}

    .. math::

        x_i=\sum_{1}^{n}\frac{I_i}{n}

    If the cut frequency of the filter is :math:`\mathit{f_c}`, the box-car length (*n*) is
    
    .. math::
    
        \frac{1}{f_c}samprate 
    
    Steps:
    
    - Declare variables
    - Define the LPF (frequency domain) and the box-car function (time domain)
    - It is going to work with a longer vector to not have fake results for the last *boxLength* windows
    - Apply the box-car window by shifting it along the (lengthened) input vector
    - Free allocated GSL vectors
    
    The function returns:
    
      - 1: Function cannot run
      - 3: Cut-off frequency too high :math:`\Rightarrow` Equivalent to not filter
      - 4: Cut-off frequency too low
    
    **Members/Variables**

    gsl_vector** **invector**

        Input/Output GSL vector (non-filtered input vector/filtered input vector) 

    int **szVct**

        Size of :cpp:member:`invector`

    int **sampleRate**

        Sampling rate (samples/s)    
        
    .. cpp:member:: gsl_vector** invector

        Input/Output GSL vector (non-filtered input vector/filtered input vector) 

    .. cpp:member:: int szVct

        Size of :cpp:member:`invector`

    .. cpp:member:: int sampleRate

        Sampling rate (samples/s)    
        
            
.. _M:

.. cpp:function:: int matrix2vector(gsl_matrix *matrixin, gsl_vector **vectorout)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function converts an input square matrix :math:`[n \times n]` into an output :math:`n^2` vector. It puts the first row of the matrix (:math:`n` elements) in the first :math:`n` elements of the vector (from :math:`0` to :math:`n-1`), the second row of the matrix in the elements from :math:`n` to :math:`2n-1` of the vector and so on.

    **Members/Variables**
    
    gsl_matrix* **matrixin**

        GSL input square matrix :math:`[n \times n]`
    
    gsl_vector** **vectorout**

        GSL output vector whose length is :math:`n^2`
        
    .. cpp:member:: gsl_matrix* matrixin

        GSL input square matrix :math:`[n \times n]`
    
    .. cpp:member:: gsl_vector** vectorout

        GSL output vector whose length is :math:`n^2`
        
        
.. cpp:function:: int medianKappaClipping (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, int boxLPF, double *threshold)
    
    Located in file: *pulseprocess.cpp*
    
    This function calculates a threshold in the first derivative of the record by using a Kappa-clipping method 
    (replacing points beyond :math:`mean\pm kappa \cdot sigma` with the median).

    Mean and sigma are calculated and values of :cpp:member:`invector` out of :math:`(mean+kappa \cdot sigma,mean-kappa \cdot sigma)` are replaced
    with the median (it is trying to look for the baseline). And this process is iteratively repeated until there are
    no points beyond :math:`mean \pm kappa \cdot sigma`. Finally, the threshold is calculated as :math:`mean+nSigmas \cdot sigma` ('+' is used because
    `if there are pulses in the input invector they are always positive`).

    Steps: 
    
    - Declare variables
    - Calculate the median
    - Iterate until there are no points out of the maximum excursion ( :math:`kappa \cdot sigma`)
    - Establish the threshold as mean+nSigmas*sigma
    
    **Members/Variables**
    
    gsl_vector* **invector**
         
        First derivative of the (filtered) record
        
    double **kappa**
         
        Value to establish the range around of the mean
        
    double **stopCriteria**
         
        It is given in %
        
    double **nSigmas**
         
        Times sigma to calculate threshold as :math:`mean+nSigmas \cdot sigma`
        
    int **boxLPF**
         
        Length of the low-pass filtering box-car
        
    double* **threshold**
         
        Calculated threshold
        
    .. cpp:member:: gsl_vector* invector
         
        First derivative of the (filtered) record
        
    .. cpp:member:: double kappa
         
        Value to establish the range around of the mean
        
    .. cpp:member:: double stopCriteria
         
        It is given in %
        
    .. cpp:member:: double nSigmas
         
        Times sigma to calculate threshold as :math:`mean+nSigmas \cdot sigma`
        
    .. cpp:member:: int boxLPF
         
        Length of the low-pass filtering box-car
        
    .. cpp:member:: double* threshold
         
        Calculated threshold
        

.. cpp:function:: int medianKappaClipping_noiseSigma (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, double *mean, double *sigma)
    
    Located in file: *gennoisespec.cpp*
    
    
    This function provides the mean and the sigma of an input vector (with noise sigmas) by using a Kappa-clipping 
    method (replacing points beyond :math:`mean\pm kappa \cdot sigma` with the median).

    First, mean and sigma are calculated and :cpp:member:`invector` values out of :math:`(mean+kappa \cdot sigma,mean-kappa \cdot sigma)` are replaced
    with the median (it is trying to look for the baseline). And this process is iteratively repeated until there are
    no points beyond :math:`mean \pm kappa \cdot sigma`. Finally, the mean and sigma of the resulting vector are provided.
    
    Steps: 
    
    - Declare variables
    - Calculate the median
    - Iterate until there are no points out of the maximum excursion ( :math:`kappa \cdot sigma`)
    - Calculate mean and sigma
    
    **Members/Variables**
    
    gsl_vector* **invector**
         
        First derivative of the (filtered) record
        
    double **kappa**
         
        Value to establish the range around of the mean
        
    double **stopCriteria**
         
        It is given in %
        
    double **nSigmas**
         
        Times sigma to calculate threshold as :math:`mean+nSigmas \cdot sigma`
        
    double* **mean**
         
        Mean value of the :cpp:member:`invector` (no points beyond :math:`mean \pm kappa \cdot sigma`)
        
    double* **sigma**
         
        Sigma value of the :cpp:member:`invector` (no points beyond :math:`mean \pm kappa \cdot sigma`)
        
    .. cpp:member:: gsl_vector* invector
         
        First derivative of the (filtered) record
        
    .. cpp:member:: double kappa
         
        Value to establish the range around of the mean
        
    .. cpp:member:: double stopCriteria
         
        It is given in %
        
    .. cpp:member:: double nSigmas
         
        Times sigma to calculate threshold as :math:`mean+nSigmas \cdot sigma`
        
    .. cpp:member:: int boxLPF
         
        Length of the low-pass filtering box-car
        
    .. cpp:member:: double* mean
         
        Mean value of the :cpp:member:`invector` (no points beyond :math:`mean \pm kappa \cdot sigma`)
        
    .. cpp:member:: double* sigma
         
        Sigma value of the :cpp:member:`invector` (no points beyond :math:`mean \pm kappa \cdot sigma`)

        
.. cpp:function:: void MyAssert(int expr, char* msg)

    Located in file: *tesreconstruction.c*

    This function displays an error message if the condition in :cpp:member:`expr` is true.
    
    **Members/Variables**
    
    int **expr**
    
        Condition to be true in order to display the error message
    
    char* msg 
    
        Message to be displayed
    

.. _N:

.. cpp:function:: extern_C_OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status)    
    
    Located in file: *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *OptimalFilterSIRENA* data structure.
    
    **Members/Variables**
    
    int* const **status**

        Input/output status
        
    .. cpp:member:: int* const status

        Input/output status
        
        
.. cpp:function:: extern_C_PulsesCollection* newPulsesCollection(int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *PulsesCollection* data structure. 
    
    **Members/Variables**
    
    int* const **status**
        
        Input/output status

    .. cpp:member:: int* const status
        
        Input/output status    
        
.. cpp:function:: extern_C_ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status)
    
    Located in file *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *ReconstructInitSIRENA* data structure.
    
    **Members/Variables**
    
    int* const **status**
        
        Input/output status
        
    .. cpp:member:: int* const status
        
        Input/output status

        
.. cpp:function:: int noDetect(gsl_vector *der, ReconstructInitSIRENA *reconstruct_init, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl, gsl_vector **samp1DERgsl) 
    
    Located in file: *pulseprocess.cpp*
    
    This function runs if the starting time of the pulses are agiven as input parameters (:option:`tstartPulse1` != 0). 
    It looks for the maximum of the derivative of the pulse and the average of the first 4 samples of the derivative of the pulse.
       
    **Members/Variables**
    
    gsl_vector* **der**
    
        First derivative of the (low-pass filtered) record

    ReconstructInitSIRENA* **reconstruct_init**
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    int* **numberPulses**
        
        Number of events
    
    gsl_vector** **tstartgsl**
        
        Starting time of the events (in samples)
    
    gsl_vector** **flagTruncated**
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    gsl_vector** **maxDERgsl**
        
        Maximum of the derivative of the event 
        
    gsl_vector** **samp1DERgsl**
        
        Average of the first 4 samples of the derivative of the event
        
    .. cpp:member:: gsl_vector* der
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    .. cpp:member:: int* numberPulses
        
        Number of events
    
    .. cpp:member:: gsl_vector** tstartgsl
        
        Starting time of the events (in samples)
    
    .. cpp:member:: gsl_vector** flagTruncated
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the derivative of the event 
        
    .. cpp:member:: gsl_vector** samp1DERgsl
        
        Average of the first 4 samples of the derivative of the event
    
    
.. _O:

.. _P:


.. cpp:function:: int obtainRiseFallTimes (gsl_vector *recordNOTFILTERED, double samprate, gsl_vector *tstartgsl, gsl_vector *tendgsl, gsl_vector *Bgsl, gsl_vector *Lbgsl, int numPulses, gsl_vector **tauRisegsl, gsl_vector **tauFallgsl)

    Located in file: *tasksSIRENA.cpp*
    
    This function provides an estimation of the rise and fall time of the detected pulses in a record.
    
    - Find the maximum of each pulse: *amax*
    - Baseline of each pulse: *abase*
    - Find the first sample in the rising part above the 50% (*amax/2*): *t2*
        - Previous and post sample to *t2*: *t1* and *t3*
        - Line by using 3 points: *(t1,a1)*, *(t2,a2)* and *(t3,a3)*
        - *t0 (t0,abase)*
        - *tmax (tmax,amax)*
        - Rise time = *tmax-t0*
    - Find the previous sample in the decreasing part to the first sample below the 50% (*amax/2*): *t2*
        - Previous and post sample to *t2*: *t3* and *t1*
        - Line by using 3 points: *(t1,a1)*, *(t2,a2)* and *(t3,a3)*
        - *t0 (t0,abase)*
        - *tmax (tmax,amax)*
        - Fall time = *t0-tmax*
       
    **Members/Variables**
    
    gsl_vector* **recordNOTFILTERED**
    
        Record neither low-pass filtered nor differentiated
        
    double **samprate**
    
        Sampling rate
        
    gsl_vector* **tstartgsl**
    
        Starting time of the detected pulses in the record (samples)
    
    gsl_vector* **tendgsl**
    
        Ending time of the detected pulses in the record (samples)
    
    gsl_vector* **Bgsl**
    
        In general, sum of the *Lb* digitized data samples of a pulse-free interval immediately before each pulse
    
    gsl_vector* **Lbgsl**
    
        Number of samples added in *Bgsl* for each pulse
    
    int **numPulses**
    
        Number of detected pulses in the record
    
    gsl_vector** **tauRisegsl**
    
        Rise time of the detected pulses in the record (seconds)
    
    gsl_vector** **tauFallgsl**
    
        Fall time of the detected pulses in the record (seconds)
    
    .. cpp:member:: gsl_vector* recordNOTFILTERED
    
        Record neither low-pass filtered nor differentiated
    
    .. cpp:member:: double samprate
    
        Sampling rate
        
    .. cpp:member:: gsl_vector* tstartgsl
    
        Starting time of the detected pulses in the record (samples)
    
    .. cpp:member:: gsl_vector* tendgsl
    
        Ending time of the detected pulses in the record (samples)
    
    .. cpp:member:: gsl_vector* Bgsl
    
        In general, sum of the *Lb* digitized data samples of a pulse-free interval immediately before each pulse
    
    .. cpp:member:: gsl_vector* Lbgsl
    
        Number of samples added in *Bgsl* for each pulse
    
    .. cpp:member:: int numPulses
    
        Number of detected pulses in the record
    
    .. cpp:member:: gsl_vector** tauRisegsl
    
        Rise time of the detected pulses in the record (seconds)
    
    .. cpp:member:: gsl_vector** tauFallgsl
    
        Fall time of the detected pulses in the record (seconds)
    

.. cpp:function:: int parabola3Pts (gsl_vector *x, gsl_vector *y, double *a, double *b, double *c)
    
    Located in file: *genutils.cpp*
    
    This function calculates the equation of a parabola given 3 points.
       
    **Members/Variables**
    
    gsl_vector* **x**
    
        Input GSL with *x* vector
        
    gsl_vector* **y**
    
        Input GSL with *y* vector
        
    double* **a**
    
        Fit coefficient of the quadratic term
        
    double* **b**
    
        Fit coefficient of the linear term    
        
    double* **c**
    
        Fit coefficient (independent term)
        
    .. cpp:member:: gsl_vector* x
    
        Input GSL with *x* vector
        
    .. cpp:member:: gsl_vector* y
    
        Input GSL with *y* vector
        
    .. cpp:member:: double* a
    
        Fit coefficient of the quadratic term
        
    .. cpp:member:: double* b
    
        Fit coefficient of the linear term    
        
    .. cpp:member:: double* c
    
        Fit coefficient (independent term)
        
        
.. cpp:function:: int polyFit(gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b, double *c)
    
    Located in file: *genutils.cpp*
    
    This function makes a polynomial fitting :math:`ax^2+bx+c` using the regression quadratic analysis. To measure how well model agrees with the data, the chi-square merit function is used, which in this case is
    
    .. math::

        \chi^2 (a,b,c)= \sum_{i=1}^{N}\left(\frac{y_i-a-bx_i-c{x_i}^2}{\sigma_i}\right)^2
        
    This equation is minimized to determine *a*, *b* and *c*. Then
    
    .. math::
    
       \begin{array}{ccc}
       S_{(x,x)}=\sum{{x_i}^2}-\frac{(\sum{x_i})^2}{N} & S_{(x^2,y)}=\sum{{x_i}^2y_i}-\frac{\sum{{x_i}^2}\cdot\sum{y_i}}{N}\\
       S_{(x,y)}=\sum{x_iy_i}-\frac{\sum{x_i}\cdot\sum{y_i}}{N} & S_{(x^2,x^2)}=\sum{{x_i}^4}-\frac{\left(\sum{{x_i}^2}\right)^2}{N}\\
       S_{(x,x^2)}=\sum{{x_i}^3}-\frac{\sum{x_i}\cdot\sum{{x_i}^2}}{N} & \\
       \end{array}

    .. math::
    
       a = \frac{S_{(x^2,y)}S_{(x,x)}-S_{(x,y)}S_{(x,x^2)}}{S_{(x,x)}S_{(x^2,x^2)} -{\vert S_{(x,x^2)} \vert}^2}

    .. math::
    
       b = \frac{S_{(x,y)}S_{(x^2,x^2)}-S_{(x^2,y)}S_{(x,x^2)}}{S_{(x,x)}S_{(x^2,x^2)} -{\vert S_{(x,x^2)} \vert}^2}
       
    .. math::
    
       c = \frac{\sum{y_i}}{N}-b\frac{\sum{x_i}}{N}-a\frac{\sum{{x_i}^2}}{N}
       
    **Members/Variables**
    
    gsl_vector* **x_fit**
    
        Input GSL with *x* vector
        
    gsl_vector* **y_fit**
    
        Input GSL with *y* vector
        
    double* **a**
    
        Fit coefficient of the quadratic term
        
    double* **b**
    
        Fit coefficient of the linear term    
        
    double* **c**
    
        Fit coefficient (independent term)
        
    .. cpp:member:: gsl_vector* x_fit
    
        Input GSL with *x* vector
        
    .. cpp:member:: gsl_vector* y_fit
    
        Input GSL with *y* vector
        
    .. cpp:member:: double* a
    
        Fit coefficient of the quadratic term
        
    .. cpp:member:: double* b
    
        Fit coefficient of the linear term    
        
    .. cpp:member:: double* c
    
        Fit coefficient (independent term)
    
    
.. cpp:function:: int polyFitLinear(gsl_vector *x_fit, gsl_vector *y_fit, double *a, double *b)
    
    Located in file: *genutils.cpp*
    
    This function makes a linear fitting :math:`ax+b` using the regression linear analysis. To measure how well model agrees with the data, the chi-square merit function is used, which in this case is
    
    .. math::

        \chi^2 (a,b)= \sum_{i=1}^{N}\left(\frac{y_i-a-bx_i}{\sigma_i}\right)^2
        
    This equation is minimized to determine *a* and *b*. Then
    
    .. math::
    
       a = \frac{N \sum{x_i y_i}- \sum{x_i}\sum{y_i}}{N \sum{x_i^2} - {\left(\sum{x_i}\right)}^2}

    .. math::
    
       b = \frac{\sum{y_i}}{N}-a \frac{\sum{x_i}}{N}
       
       
    **Members/Variables**
    
    gsl_vector* **x_fit**
    
        Input GSL with *x* vector
        
    gsl_vector* **y_fit**
    
        Input GSL with *y* vector
        
    double* **a**
    
        Fit coefficient of the linear term
        
    double* **b**
    
        Fit coefficient (independent term)    
        
    .. cpp:member:: gsl_vector* x_fit
    
        Input GSL with *x* vector
        
    .. cpp:member:: gsl_vector* y_fit
    
        Input GSL with *y* vector
        
    .. cpp:member:: double* a
    
        Fit coefficient of the linear term
        
    .. cpp:member:: double* b
    
        Fit coefficient (independent term)    

        
.. cpp:function:: void print_error( const char* const func, string message, int status)

    Located in file: *genutils.cpp*

    This function prints out error messages.

    **Members/Variables**
    
    const char* const **func**
    
        Function name whose error is printed 
        
    string **msg**
    
        Error message to be printed 

    int **status**
        
        Status
        
    .. cpp:member:: const char* const func
    
        Function name whose error is printed 
        
    .. cpp:member:: string msg
    
        Error message to be printed 

    .. cpp:member:: int status
        
        Status

        
.. cpp:function:: int procRecord (ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, gsl_vector *recordWithoutConvert2R, PulsesCollection *foundPulses)
        
    Located in file: *tasksSIRENA.cpp*
    
    This function processes the input record (detecting the pulses):

    1) Declare and initialize variables
    
    2) Allocate GSL vectors
    
    3) (Low-pass filtering and) differentiation
    
    4) Find the events (pulses) in the record
    
       - If production mode (:option:`opmode` = 1):
       
            - No detect if :option:`tstartPulse1` != 0: 'noDetect' 
            - Detect (:option:`tstartPulse1` != 0): 
            
                - 'InitialTriggering'
                - 'FindSecondaries' (:option:`detectionMode` = AD) or 'FindSecondariesSTC' (:option:`detectionMode` = STC)
                
       - If calibration mode (:option:`opmode` = 0): 'findPulsesCAL'
    
    5) Calculate the end time of the found pulses and check if the pulse is saturated
    
    6) Calculate the baseline (mean and standard deviation) before a pulse (in general *before*) :math:`\Rightarrow` To be written in **BSLN** and **RMSBSLN** columns in the output FITS file
    
    7) Obtain the approximate rise and fall times of each pulse 
    
    8) Load the found pulses data in the input/output *foundPulses* structure
    
    9) Write test info (if *reconstruct_init->intermediate* = 1)
    
    10) Write pulses info in intermediate output FITS file (if *reconstruct_init->intermediate* = 1)
    
    11) Free allocated GSL vectors

    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    double **tstartRecord** 

        Starting time of the record (in order to calculate absolute times)

    double **samprate** 

        Sampling rate (in order to low-pass filter)

    fitsfile* **dtcObject**

        Object which contains information of the intermediate FITS file (to be written if :option:`intermediate` = 1)

    gsl_vector* **record**

        GSL vector with signal values of input record

    gsl_vector* **recordWithoutConvert2R**

        GSL vector with original signal values of input record (without being converted to R space)
        
    PulsesCollection* **foundPulses** 

        Input/output structure where the info about found pulses is stored
        
    long **num_previousDetectedPulses** 

        Number of previous detected pulses (to know the index to get the proper element from *tstartPulse1_i* in case :option:`tstartPulse1` was a file name)
        
    int **pixid** 

        Pixel ID (from the input file) to be propagated
    
    int **phid**

        Photon ID (from the input file) to be propagated
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: double tstartRecord 

        Starting time of the record (in order to calculate absolute times)

    .. cpp:member:: double samprate 

        Sampling rate (in order to low-pass filter)

    .. cpp:member:: fitsfile* dtcObject

        Object which contains information of the intermediate FITS file (to be written if :option:`intermediate` = 1)

    .. cpp:member:: gsl_vector* record

        GSL vector with signal values of input record

    .. cpp:member:: gsl_vector* recordWithoutConvert2R

        GSL vector with original signal values of input record (without being converted to R space)
        
    .. cpp:member:: PulsesCollection* foundPulses 

        Input/output structure where the info about found pulses is stored
        
    .. cpp:member:: long num_previousDetectedPulses 

        Number of previous detected pulses (to know the index to get the proper element from *tstartPulse1_i* in case :option:`tstartPulse1` was a file name)
        
    .. cpp:member:: int pixid 

        Pixel ID (from the input file) to be propagated
    
    .. cpp:member:: int phid

        Photon ID (from the input file) to be propagated
        
                
.. cpp:function:: int pulseGrading(ReconstructInitSIRENA *reconstruct_init, int grade1, int grade2, int OFlength_strategy, int *pulseGrade, long *OFlength)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function provides the pulse grade (Pileup=-2, Rejected=-1, HighRes=1, MidRes=2, LimRes=3, LowRes=4) and the optimal filter length by taking into account the info read from the XML file and the :option:`OFStrategy` (**FREE**, **BYGRADE** or **FIXED**).

    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).
     
    int **grade1**
     
        Pulse duration (length of optimal filter applied)
     
    int **grade2**
    
        Difference between the start time of the pulse and the start time of the previous pulse
        
    int **OFlength_strategy**
     
        Same as :option:`OFStrategy` (input)
        
    int* **pulseGrade**
    
        Pulse grade (output)
    
    long* **OFlength**
    
        Optimal filter length (= :option:`OFLength` only if :option:`OFStrategy` = **FIXED** and :option:`OFLength` <= grade1) (output)
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).
     
    .. cpp:member:: int grade1
     
        Pulse duration (length of optimal filter applied)
     
    .. cpp:member:: int grade2
    
        Difference between the start time of the pulse and the start time of the previous pulse
        
    .. cpp:member:: int OFlength_strategy
     
        Same as :option:`OFStrategy` (input)
        
    .. cpp:member:: int* pulseGrade
    
        Pulse grade (output)
    
    .. cpp:member:: long* OFlength
    
        Optimal filter length (= :option:`OFLength` only if :option:`OFStrategy` = **FIXED** and :option:`OFLength` <= grade1) (output)
        
        
.. _Q:

.. _R:

.. cpp:function:: int readAddOrderParams(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int eventcntLib, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, gsl_vector *pulsetemplateMaxLengthFixedFilter)
    
    Located in file: *tasksSIRENA.cpp*

    This function reads the library data, add new data (a new row) and sort the data according to an energy-ascending order.

    - Declare variables
    - Load values already in the library
    - Add new values
    - Realign
    - Add intermeadiate values
    - Recalculate intermediate values of some new pairs
    - Write values in the library
    - Free allocated GSL vectors
    
    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
     
    fitsfile** **inLibObject**

        FITS object containing information of the library FITS file  
        
    double **samprate**

        Sampling rate

    int **eventcntLib** 

        Number of templates in the library
        
    double **estenergy** 

        Pulse height of the template whose energy is going to be added to the library
        
    gsl_vector* **pulsetemplate**

        GSL vector with the pulse template whose energy is going to be added to the library
        
    gsl_matrix* **covariance**

        GSL matrix with covariance matrix of the energy which is going to be added to the library
    
    gsl_matrix* **weight**

        GSL matrix with weight matrix of the energy which is going to be added to the library
        
    gsl_vector* **pulsetemplateMaxLengthFixedFilter**

        GSL vector with the :option:`largeFilter`-length template whose energy is going to be added to the library
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
     
    .. cpp:member:: fitsfile** inLibObject

        FITS object containing information of the library FITS file  
        
    .. cpp:member:: double samprate

        Sampling rate

    .. cpp:member:: int eventcntLib 

        Number of templates in the library
        
    .. cpp:member:: double estenergy 

        Pulse height of the template whose energy is going to be added to the library
        
    .. cpp:member:: gsl_vector* pulsetemplate

        GSL vector with the pulse template whose energy is going to be added to the library
        
    .. cpp:member:: gsl_matrix* covariance

        GSL matrix with covariance matrix of the energy which is going to be added to the library
    
    .. cpp:member:: gsl_matrix* weight

        GSL matrix with weight matrix of the energy which is going to be added to the library
	
    .. cpp:member:: gsl_vector* pulsetemplateMaxLengthFixedFilter

        GSL vector with the :option:`largeFilter`-length template whose energy is going to be added to the library
        
        
.. cpp:function:: int readFitsComplex(IOData obj, gsl_matrix **result)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a complex column of a FITS file. After that, the function puts them into a GSL matrix for an easier processing.
    
    **Members/Variables**
    
    IOData **obj**
    
        Input object for complex FITS column
    
    gsl_matrix** **result**
    
        Output GSL matrix
        
    .. cpp:member:: IOData obj
    
        Input object for complex FITS column
    
    .. cpp:member:: gsl_matrix** result
    
        Output GSL matrix
        
        
.. cpp:function:: int readFitsSimple(IOData obj, gsl_vector **result)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a simple column of a FITS file. After that, the function puts them into a GSL vector for an easier processing.
    
    **Members/Variables**
    
    IOData **obj**
    
        Input object for simple FITS column
    
    gsl_vector** **result**
    
        Output GSL vector

    .. cpp:member:: IOData obj
    
        Input object for simple FITS column
    
    .. cpp:member:: gsl_vector** result
    
        Output GSL vector
    
.. cpp:function:: extern_C_void reconstructRecordSIRENA(TesRecord* record, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init,  int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function is the main wrapper function to detect, grade and calculate the energy of the pulses in the input records.
    
    - Inititalize *PulsesCollection* structure
    - Check consistency of some input parameters
    - If first record, read the necessary keywords and columns from the input file in order to convert from current to quasi-resistance space
    - In case of running with threading
    - Detect pulses in input record (:cpp:func:`runDetect`). 
    - If reconstruction (:option:`opmode` = 1) and not PCA:
        - Filter and calculate energy of pulses (:cpp:func:`runEnergy`)
    - Fill in the :cpp:member:`pulsesAll` structure
    - Populate output event list with pulses energies, arrival time and grading
    
    **Members/Variables**
    
    TesRecord* **record** 

        Instance of *TesRecord* structure that contains the input record 
        
    int **trig_reclength**
    
        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)
    
    TesEventList* **event_list**
        
        Instance of *TesEventList* structure that contains the information of the reconstructed pulses
        
    ReconstructInitSIRENA* **reconstruct_init**
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    int **lastRecord**
    
        If record being analyzed is the last one, :cpp:member:`lastRecord` = 1. Otherwise it is equal to 0
        
    int **nRecord**
    
        Input record number
    
    PulsesCollection** **pulsesAll**
    
        Member of *PulsesCollection* structure to successively store all the pulses used to create the library. Re-populated after each processed record.
    
    OptimalFilterSIRENA** **optimalFilter**
    
        Optimal filters used in reconstruction
    
    int* const **status**
    
        Input/output status
        
    .. cpp:member:: TesRecord* record 

        Instance of *TesRecord* structure that contains the input record 
        
    .. cpp:member:: int trig_reclength
    
        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)
    
    .. cpp:member:: TesEventList* event_list
        
        Instance of *TesEventList* structure that contains the information of the reconstructed pulses
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    .. cpp:member:: int lastRecord
    
        If record being analyzed is the last one, :cpp:member:`lastRecord` = 1. Otherwise it is equal to 0
        
    .. cpp:member:: int nRecord
    
        Input record number
    
    .. cpp:member:: PulsesCollection** pulsesAll
    
        Member of *PulsesCollection* structure to successively store all the pulses used to create the library. Re-populated after each processed record.
    
    .. cpp:member:: OptimalFilterSIRENA** optimalFilter
    
        Optimal filters used in reconstruction
    
    .. cpp:member:: int* const status
    
        Input/output status
        

.. cpp:function:: int RS_filter(gsl_vector *vector, double lrs, double lb, double B, double *pulseheight)
    
    Located in file: *pulseprocess.cpp*

    This function uses the running sum filter to find the pulse height. It always works in time domain.

    A running sum filter, *RS*, is the sum of :cpp:member:`lrs` digitized data samples. It is continuously updated upon the arrival of
    new data point. Simultaneously a baseline filter, :cpp:member:`B`, is the sum of :cpp:member:`lb` digitized data samples without pulses. The 
    algorithm looks for the time when *RS/lrs* reaches its maximum. At that time *RS* is stored, :math:`RS_{max}`, and the baseline
    is scaled with :cpp:member:`lrs`, *Bp* ( :math:`Bp=B \cdot lrs/lb`). Then, the pulse height related to the pulse pseudoenergy is given by:
    
    .. math::
    
        Pulse height=\frac{RS_{max}-B_p}{lrs}
       
    **Members/Variables**
    
    gsl_vector* **vector**

        Not filtered pulse (extracted from the record in :cpp:func:`getPulseHeight`)
        
    double **lrs**

        Running sum length (samples)
        
    double **lb**

        Baseline averaging length (samples)
        
    double **B**

        In general, sum of the :cpp:member:`lb` digitized data samples of a pulse-free interval immediately before the current pulse
        
    double* **pulseheight**

        Pulseheight of the input pulse    
        
    .. cpp:member:: gsl_vector* vector

        Not filtered pulse (extracted from the record in :cpp:func:`getPulseHeight`)
        
    .. cpp:member:: double lrs

        Running sum length (samples)
        
    .. cpp:member:: double lb

        Baseline averaging length (samples)
        
    .. cpp:member:: double B

        In general, sum of the :cpp:member:`lb` digitized data samples of a pulse-free interval immediately before the current pulse
        
    .. cpp:member:: double* pulseheight

        Pulseheight of the input pulse    
        
        
.. cpp:function:: void runDetect(TesRecord* record, int trig_reclength, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
     
    Located in file: *tasksSIRENA.cpp*

    This function is responsible for the **detection** in SIRENA, record by record. It is used both for library creation (:option:`opmode` = 0) and energy reconstruction (:option:`opmode` = 1) runnings.

    Conditions:   

        - If first record and :option:`opmode` = 1  :math:`\Rightarrow`  Run :cpp:func:`filderLibrary`

        - If last record and :option:`opmode` = 0 :math:`\Rightarrow` Run :cpp:func:`calculateTemplate` and :cpp:func:`writeLibrary`

        - If :option:`intermediate` = 1 :math:`\Rightarrow` :cpp:func:`writeTestInfo` and :cpp:func:`writePulses`

        - If :option:`opmode` = 0 :math:`\Rightarrow` Find pulses by using :cpp:func:`findPulsesCAL`
        
        - If :option:`opmode` = 1 :math:`\Rightarrow` Find pulses by :cpp:func:`InitialTriggering` and :cpp:func:`FindSecondaries` or :cpp:func:`FindSecondariesSTC`
        
    Steps:

        1) Create library file if it is necessary: calibration (:option:`opmode` = 0) and last record (run :cpp:func:`createLibrary`)           

        2) Create intermediate output FITS file if required (:cpp:func:`createDetectFile`)     

        3) (Filter and) differentiate the *models* of the library (only for the first record in :option:`opmode` = 1). Run  (:cpp:func:`filderLibrary`)

        4) Store the input record in *invector* (:cpp:func:`loadRecord`)

        5) Convert *I* into *R* if :option:`EnergyMethod` = **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED** (:cpp:func:`convertI2R`)
        
        6) Process each record (:cpp:func:`proceRecord`): 

                - (Low-pass filter and) differentiate                          
                - Find pulses                                                           
                - Load the found pulses data in the input/output *foundPulses* structure
                - Write test info in intermediate output FITS file if :option:`intermediate` = 1 (:cpp:func:`writeTestInfo`)                            
                - Write pulses info in intermediate output FITS file if :option:`intermediate` = 1 (:cpp:func:`writePulses`)
                
        **From this point forward, I2R, I2RALL, I2RNOL and I2RFITTED are completely equivalent to OPTFILT**     
           
        7) If last record in :option:`opmode` = 0 run:                                         

                * :cpp:func:`calculateTemplate` (and :cpp:func:`weightMatrix`)
                * :cpp:func:`writeLibrary`
                
        8) If last record and PCA:
        
                - In order to not have restrictions when providing (\*reconstruct_init)->energyPCAx
                - Covariance data
                - Eigenvalues and eigenvectors
                - RSxN (S=2)
                - AE straight line: Pto0(x,y) and Pto10(x,y)
                - Calculus of the rotation angle
                - Rotation
                - Histograms of the two clusters (two energies)
                - Conversion factor from arbitrary unit to eV
                - Energy calculation

        9) Close intermediate output FITS file if it is necessary   

    **Members/Variables**

    TesRecord* **record** 

        Member of *TesRecord* structure that contains the input record 
        
    int **trig_reclength** 

        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)

    int **lastRecord**

        Integer to verify whether *record* is the last one (=1) to be read (and thus if library file will be created)
        
    PulsesCollection* **pulsesAll**

        Member of *PulsesCollection* structure to successively store all the pulses used to create the library. Re-populated after each processed record
        
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    PulsesCollection** **pulsesInRecord**

        Member of *PulsesCollection* structure to store all the pulses found in the input record
        
    .. cpp:member:: TesRecord* record 

        Member of *TesRecord* structure that contains the input record 
        
    .. cpp:member:: int trig_reclength 

        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)

    .. cpp:member:: int lastRecord

        Integer to verify whether *record* is the last one (=1) to be read (and thus if library file will be created)
        
    .. cpp:member:: PulsesCollection* pulsesAll

        Member of *PulsesCollection* structure to successively store all the pulses used to create the library. Re-populated after each processed record
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: PulsesCollection** pulsesInRecord

        Member of *PulsesCollection* structure to store all the pulses found in the input record

                
.. cpp:function:: void runEnergy(TesRecord* record, int trig_reclength, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the pulse energy applying different methods (from :option:`EnergyMethod` and :option:`OFNoise`).
    It only runs in RECONSTRUCTION mode (:option:`opmode` = *1*) (except to :option:`EnergyMethod` = **PCA**).

    - Declare variables
    - Store the :cpp:member:`record` in *invector* (:cpp:func:`loadRecord`)
    - Subtract the baseline if :option:`EnergyMethod` = **OPTFILT** and *runF0orB0val* = 1 (:option:`FilterMethod` = **B0**)
    - Subtract the baseline if :option:`EnergyMethod` = **WEIGHT**
    - Check Quality
    - For each pulse:

        - Establish the pulse grade (HighRes=1, MidRes=2, LimRes=3, LowRes=4, Rejected=-1, Pileup=-2) and the optimal filter length
        - Pulse: Load the proper piece of the record in *pulse*
        - Get the low resolution energy estimator by filtering with a 4-samples-length filter:
            - Load the low resolution pulse in *pulse_lowres*
            - Get the filter
            - Calculate the low resolution estimator
        - If :option:`OFIter` = 1, in the first iteration ( *numiteration* = 0) the values of *maxDER* and *maxDERs* are used in 
          :cpp:func:`find_matchedfilterDAB`, :cpp:func:`find_optimalfilterDAB` or :cpp:func:`find_Esboundary` getting the values of the *energies* which straddle the *maxDER* (*Ealpha* and *Ebeta*). It will have more iterations if the calculated *energy* is out of *[Ealpha, Ebeta]*. If *energy* is in *[Ealpha, Ebeta]* the iterative process stops.
            
                - If :option:`EnergyMethod` = **OPTFILT** (or **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED**) and :option:`OFLib` = 0 and :option:`OFNoise` = **NSD**:
                
                    - Find the matched filter and load it in *filter* (:cpp:func:`find_matchedfilterDAB`) 
                    - Calculate the optimal filter
                    
                - If :option:`EnergyMethod` = **OPTFILT** (or **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED**) and :option:`OFLib` = 1 and :option:`OFNoise` = **NSD**:
                
                    - If it is necessary, choose the base-2 system value closest (lower than or equal) to the pulse length
                    - Find the optimal filter and load it in *filter* (:cpp:func:`find_optimalfilterDAB`)
                    
                - If :option:`EnergyMethod` = **WEIGHT** or **WEIGHTN**:
                
                    - Get the indexes of the two energies which straddle the pulse (:cpp:func:`find_Esboundary`)
                    - If :option:`EnergyMethod` = **WEIGHTN** and :option:`OFLib` = 1:
                    
                       - Choose the base-2 system value closest (lower than or equal) to the pulse length
                       - :cpp:func:`find_prclwn` to find the appropriate values of the *PRECALWN* HDU (**PCLx** columns)
                
                - If :option:`EnergyMethod` = **OPTFILT** (or **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED**) and :option:`OFLib` = 1 and :option:`OFNoise` = **WEIGHTM**:
                
                    - Choose the base-2 system value closest (lower than or equal) to the pulse length
                    - :cpp:func:`find_prclofwm` to find the appropriate values of the *PRCLOFWM* HDU (**OFWx** columns)
                
                - Subtract the sum of the filter if :option:`EnergyMethod` = **OPTFILT**, :option:`OFNoise` = **NSD**, :option:`FilterDomain` = **T**, 0-padding and :option:`Sum0Filt` =1 
                - Calculate the energy of each pulse
                - If using lags, it is necessary to modify the tstart of the pulse 
        - In order to subtract the pulse model, it has to be located in the tstart with jitter and know its values in the digitized samples
        - Subtract the pulse model from the record
        - Write info of the pulse in the output intemediate file if :option:`intermediate` = 1
    - Free allocated GSL vectors

    **Members/Variables**

    TesRecord** **record**

        Structure that contains the input ADC record
        
    int **trig_reclength** 

        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)
        
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    PulsesCollection* **pulsesInRecord**

        Collection of pulses found in the current record
        
    OptimalFilterSIRENA** **optimalFilter**
    
        Optimal filters used in reconstruction
        
     PulsesCollection* **pulsesAll**
    
        Member of *PulsesCollection* structure to store all the pulses found in the input FITS file. To know the index to get the proper element from *tstartPulse1_i* in case :option:`tstartPulse1` was a file name
        
    .. cpp:member:: TesRecord** record

        Structure that contains the input ADC record
        
    .. cpp:member:: int trig_reclength 

        Record size (just in case threading and input files with different **ADC** lengths but the same record size indeed)
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    .. cpp:member:: PulsesCollection* pulsesInRecord

        Collection of pulses found in the current record
        
    .. cpp:member:: OptimalFilterSIRENA** optimalFilter
    
        Optimal filters used in reconstruction
        
    .. cpp:member::  PulsesCollection* pulsesAll
    
        Member of *PulsesCollection* structure to store all the pulses found in the input FITS file. To know the index to get the proper element from *tstartPulse1_i* in case :option:`tstartPulse1` was a file name
        
        
.. cpp:function:: void th_runEnergy(TesRecord* record, int trig_reclength, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
    
    Located in file: *tasksSIRENA.cpp*


    This function is responsible for the **reconstruction** in SIRENA (instead of :cpp:func:`runEnergy`) when the **THREADING** running option has been chosen (hardcoded at this moment). 

.. _S:

.. cpp:function:: int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
    
    Located in file: *tasksSIRENA.cpp*

    This function returns as :cpp:member:`vectorout` the :cpp:member:`vectorin` delayed :cpp:member:`m` samples.

    **Members/Variables**

    int **m**

        Delay in samples

    gsl_vector* **vectorin**

        GSL vector with input vector
        
    gsl_vector* **vectorout** 

        GSL with input vector (:cpp:member:`vectorin`) delayed :cpp:member:`m` samples

    .. cpp:member:: int m

        Delay in samples

    .. cpp:member:: gsl_vector* vectorin

        GSL vector with input vector
        
    .. cpp:member:: gsl_vector* vectorout 

        GSL with input vector (:cpp:member:`vectorin`) delayed :cpp:member:`m` samples
        
        
.. cpp:function:: int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m)
        
    Located in file: *tasksSIRENA.cpp*

    This function returns as :cpp:member:`vectorout` the :cpp:member:`vectorin` moved forward by :cpp:member:`m` samples.
    
    **Members/Variables**
    
    int **m**

        Advance in samples

    gsl_vector* **vectorin**

        GSL vector with input vector

    gsl_vector* **vectorout** 

        GSL with input vector (:cpp:member:`vectorin`) moved forward :cpp:member:`m` samples
        
    .. cpp:member:: int m

        Advance in samples

    .. cpp:member:: gsl_vector* vectorin

        GSL vector with input vector

    .. cpp:member:: gsl_vector* vectorout 

        GSL with input vector (:cpp:member:`vectorin`) moved forward :cpp:member:`m` samples
        
            
.. _T:

.. cpp:function:: int tesreconstruction_main()

    Located in file: *tesreconstruction.c*
    
    This function is mainly a wrapper to pass a data file to the SIRENA tasks in order to reconstruct the energy of the incoming X-ray photons after their detection.
    
    It can run the SIRENA tasks or the Philippe Peille's tasks depending on the 'Rcmethod' selected.
    
    Steps:
    
    - Register HEATOOL
    - Reading all programm parameters by using PIL
    - Read XML info
    - Read the grading info from the input XML file
    - Obtain the samplig rate and the 'trig_reclength':
        - If Rcmethod starts with '@' :math:`\Rightarrow` List of record input FITS files. For every FITS file:
            - Open FITS file
            - Check if input FITS file have been simulated with TESSIM or XIFUSIM
            - If it is a xifusim simulated file
                - Obtain 'trig_reclength' from the ``HISTORY`` block
        - If Rcemethod doesn't start with '@' :math:`\Rightarrow` Single record input FITS file
            - Open FITS file
            - Check if input FITS file have been simulated with TESSIM or XIFUSIM
            - If it is a xifusim simulated file
                - Obtain 'trig_reclength' from the ``HISTORY`` block
    - Sixt standard keywords structure
    - Open output FITS file
    - Initialize PP data structures needed for pulse filtering
    - Initialize SIRENA data structures needed for pulse filtering
    - Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
    - Build up TesEventList to recover the results of the reconstruction
    - Reconstruct the input record FITS file:
        - If Rcmethod starts with '@' :math:`\Rightarrow` List of record input FITS files. For every FITS file:
            - Open record file
            - Initialize: initializeReconstruction or initializeReconstructionSIRENA
            - Build up TesRecord to read the file
            - Iterate of records and do the reconstruction
                - Reconstruct: reconstructRecord or reconstructRecordSIRENA
                - Save events to the event_list
                - Copy trigger keywords to event file
                - Close file
        - If Rcemethod doesn't start with '@' :math:`\Rightarrow` Single record input FITS file
            - Open record file
            - Initialize: initializeReconstruction or initializeReconstructionSIRENA
            - Build up TesRecord to read the file
            - Iterate of records and do the reconstruction
                - Reconstruct: reconstructRecord or reconstructRecordSIRENA
                - Save events to the event_list
                - Copy trigger keywords to event file
                - Close file
    - Save GTI extension to event file
    - Free memory
    
    The user must supply the following input parameters (*tesreconstruction.par* file).
    
    Common parameters:

    char **Rcmethod** 
        
        Reconstruction method (**PP** or **SIRENA**). If SIRENA :math:`\Rightarrow` If Rcmethod starts with '@' it provides a file text containing several record input FITS files
        
    char **RecordFile**
    
        Record FITS file
        
    char **TesEventFile**
    
        Output event list file
        
    int **PulseLength**
    
        Pulse length
        
    int **EventListSize**
    
        Default size of the event list
        
    char **clobber**
    
        Overwrite or not output files if exist (1/0)
        
    char **history**
    
        Write program parameters into output file
        
    PP parameters:
    
    double **SaturationValue**
    
        Saturation level of the ADC curves
        
    char **OptimalFilterFile**
    
        Optimal filters file
        
    char **PulseTemplateFile**
    
        Pulse template file
        
    double **Threshold**
    
        Threshold level
    
    double **Calfac**
    
        Calibration factor (should be read from the xml file)
    
    int **NormalExclusion**
    
        Minimal distance before using OFs after a misreconstruction
    
    int **DerivateExclusion**
    
        Minimal distance before reconstructing any event after a misreconstruction

    SIRENA parameters:
    
    char **LibraryFile**
    
        File with calibration library
        
    double **scaleFactor**
    
        Detection scale factor for initial filtering
        
    int **samplesUp**
    
        Number of consecutive samples up for threshold trespassing (only used in calibration run, and in production run with **STC** detection mode)
    
    int **samplesDown**
    
        Number of consecutive samples below the threshold to look for other pulse (only used in production run with **STC** detection mode)
    
    double **nSgms**
    
        Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm
        
    int **detectSP**
    
        Detect secondary pulses (1) or not (0)
    
    double **LrsT**
    
        Running sum length for the RS raw energy estimation (seconds) (only for library creation)

    double **LbT**
    
        Baseline averaging length (seconds)
    
    double **monoenergy**
    
        Monochromatic energy of the pulses in the input FITS file in eV (only for library creation)

    char **hduPRECALWN**
    
        Add or not the PRECALWN HDU in the library file (1/0) (only for library creation)

    char **hduPRCLOFWM**
    
        Add or not the PRCLOFWM HDU in the library file (1/0) (only for library creation)
        
    int **largeFilter**
    
        Length of the longest fixed filter (only for library creation)

    int **opmode**
    
        Calibration run (0) or energy reconstruction run (1)
        
    char **detectionMode**
    
        Adjusted Derivative (**AD**) or Single Threshold Crossing (**STC**)

    char **NoiseFile**
    
        Noise FITS file with noise spectrum
        
    char **FilterDomain**
    
        Filtering Domain: Time (**T**) or Frequency (**F**)
        
    char **FilterMethod** 
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline)
        
    char **EnergyMethod**
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**
        
    double **filtEeV**
    
        Energy of the filters of the library to be used to calculate energy (only for **OPTFILT**, **I2R**, **I2RALL**, **I2RNOL** and **I2RFITTED**)
        
    char **OFNoise**
    
        Noise to use with Optimal Filtering: **NSD** or **WEIGHTM**
        
    int **LagsOrNot**
    
        Lags or no lags (1/0)
        
    int **nLags**
    
        Number of lags (positive odd number)
        
    int **Fitting35**
    
        Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
    
    int **OFIter**
    
        Iterate or not iterate (1/0)
    
    int **OFLib**
    
        Work or not with a library with optimal filters (1/0)
        
    char **OFStrategy** 
    
        Optimal Filter length Strategy: **FREE**, **BYGRADE** or **FIXED**
        
    int **OFLength**
    
        Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**)
        
    int **preBuffer**
    
        Some samples added before the starting time of a pulse
        
    int **intermediate**
    
        Write or not intermediate files (1/0)
        
    char **detectFile**
    
        Intermediate detections file (if intermediate*=1)
        
    int **errorT**
    
        Additional error (in samples) added to the detected time (Logically, it changes the reconstructed energies )
        
    int **Sum0Filt**
    
        0-padding: Subtract the sum of the filter (1) or not (0)
        
    char **tstartPulse1** 
    
        Integer number: Sample where the first pulse starts or nameFile: File where the tstart (seconds) of every pulse is
    
    int **tstartPulse2**
    
        Tstart (samples) of the second pulse
        
    int **tstartPulse3**
    
        Tstart (samples) of the third pulse (if 0 :math:`\Rightarrow` PAIRS, if not 0 :math:`\Rightarrow` TRIOS)
        
    double **energyPCA1**
    
        First energy (only for PCA)
        
    double **energyPCA2**
    
        Second energy (only for PCA)
    
    char **XMLFile**
    
        XML input FITS file with instrument definition
        
    .. cpp:member:: char Rcmethod
        
        Reconstruction method (**PP** or **SIRENA**)
        SIRENA: If Rcmethod starts with '@' it provides a file text containing several record input FITS files
        
    .. cpp:member:: char RecordFile
    
        Record FITS file
        
    .. cpp:member:: char TesEventFile
    
        Output event list file
        
    .. cpp:member:: int PulseLength
    
        Pulse length
        
    .. cpp:member:: int EventListSize
    
        Default size of the event list
        
    .. cpp:member:: char clobber
    
        Overwrite or not output files if exist (1/0)
        
    .. cpp:member:: char history
    
        Write program parameters into output file
        
    PP parameters:
    
    .. cpp:member:: double SaturationValue
    
        Saturation level of the ADC curves
        
    .. cpp:member:: char OptimalFilterFile 
    
        Optimal filters file
        
    .. cpp:member:: char PulseTemplateFile
    
        Pulse template file
        
    .. cpp:member:: double Threshold
    
        Threshold level
    
    .. cpp:member:: double Calfac
    
        Calibration factor (should be read from the xml file)
    
    .. cpp:member:: int NormalExclusion
    
        Minimal distance before using OFs after a misreconstruction
    
    .. cpp:member:: int DerivateExclusion
    
        Minimal distance before reconstructing any event after a misreconstruction

    SIRENA parameters:
    
    .. cpp:member:: char LibraryFile
    
        File with calibration library
        
    .. cpp:member:: double scaleFactor
    
        Detection scale factor for initial filtering
        
    .. cpp:member:: int samplesUp
    
        Number of consecutive samples up for threshold trespassing (only used in calibration run, and in production run with **STC** detection mode)
    
    .. cpp:member:: int samplesDown
    
        Number of consecutive samples below the threshold to look for other pulse (only used in production run with **STC** detection mode)
    
    .. cpp:member:: double nSgms
    
        Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm
        
    .. cpp:member:: int detectSP
    
        Detect secondary pulses (1) or not (0)
    
    .. cpp:member:: double LrsT
    
        Running sum length for the RS raw energy estimation (seconds) (only for library creation)

    .. cpp:member:: double LbT
    
        Baseline averaging length (seconds)
    
    .. cpp:member:: double monoenergy
    
        Monochromatic energy of the pulses in the input FITS file in eV (only for library creation)

    .. cpp:member:: char hduPRECALWN
    
        Add or not the PRECALWN HDU in the library file (1/0) (only for library creation)

    .. cpp:member:: char hduPRCLOFWM
    
        Add or not the PRCLOFWM HDU in the library file (1/0) (only for library creation)
        
    .. cpp:member:: int largeFilter
    
        Length of the longest fixed filter (only for library creation)

    .. cpp:member:: int opmode
    
        Calibration run (0) or energy reconstruction run (1)
        
    .. cpp:member:: char detectionMode
    
        Adjusted Derivative (**AD**) or Single Threshold Crossing (**STC**)

    .. cpp:member:: char NoiseFile
    
        NoiseFile: Noise FITS file with noise spectrum
        
    .. cpp:member:: char FilterDomain
    
        Filtering Domain: Time (**T**) or Frequency (**F**)
        
    .. cpp:member:: char FilterMethod 
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline)
        
    .. cpp:member:: char EnergyMethod
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**
        
    .. cpp:member:: double filtEeV
    
        Energy of the filters of the library to be used to calculate energy (only for **OPTFILT**, **I2R**, **I2RALL**, **I2RNOL** and **I2RFITTED**)
        
    .. cpp:member:: char OFNoise
    
        Noise to use with Optimal Filtering: **NSD** or **WEIGHTM**
        
    .. cpp:member:: int LagsOrNot
    
        Lags or no lags (1/0)
        
    .. cpp:member:: int nLags
    
        Number of lags (positive odd number)
        
    .. cpp:member:: int Fitting35
    
        Number of lags to analytically calculate a parabola (3) or to fit a parabola (5)
    
    .. cpp:member:: int OFIter
    
        Iterate or not iterate (1/0)
    
    .. cpp:member:: int OFLib
    
        Work or not with a library with optimal filters (1/0)
        
    .. cpp:member:: char OFStrategy 
    
        Optimal Filter length Strategy: **FREE**, **BYGRADE** or **FIXED**
        
    .. cpp:member:: int OFLength
    
        Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**)
        
    .. cpp:member:: int preBuffer
    
        Some samples added before the starting time of a pulse
        
    .. cpp:member:: int intermediate
    
        Write or not intermediate files (1/0)
        
    .. cpp:member:: char detectFile 
    
        Intermediate detections file (if intermediate*=1)
        
    .. cpp:member:: int errorT
    
        Additional error (in samples) added to the detected time (Logically, it changes the reconstructed energies )
        
    .. cpp:member:: int Sum0Filt
    
        0-padding: Subtract the sum of the filter (1) or not (0)
        
    .. cpp:member:: char tstartPulse1 
    
        Integer number: Sample where the first pulse starts or nameFile: File where the tstart (seconds) of every pulse is
    
    .. cpp:member:: int tstartPulse2
    
        Tstart (samples) of the second pulse
        
    .. cpp:member:: int tstartPulse3
    
        Tstart (samples) of the third pulse (if 0 :math:`\Rightarrow` PAIRS, if not 0 :math:`\Rightarrow` TRIOS)
        
    .. cpp:member:: double energyPCA1
    
        First energy (only for PCA)
        
    .. cpp:member:: double energyPCA2
    
        Second energy (only for PCA)
    
    .. cpp:member:: char XMLFile
    
        XML input FITS file with instrument definition
        
    
.. cpp:function:: int th_runDetect (TesRecord* record, int trig_reclength, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
     
    Located in file: *tasksSIRENA.cpp*

    This function is responsible for the **detection** in SIRENA (instead of :cpp:func:`runDetect`) when the **THREADING** running option has been chosen (hardcoded at this moment). It is used both for library creation (:option:`opmode` = 0) and energy reconstruction (:option:`opmode` = 1) runnings.
    
.. cpp:function:: int toGslMatrix(void **buffer, gsl_matrix **matrix, long numCol, int numRow, int type, int eventini)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input buffer into an output GSL matrix. Columns and rows are input parameters.
    
    **Members/Variables**
    
    void** **buffer**
    
        Input buffer with data
        
    gsl_matrix** **matrix**
    
        Output GSL matrix 
        
    long **numCol**
    
        Number of columns
        
    int **numRow**
    
        Number of rows
        
    int **type**
    
        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
    int **eventini**
    
        Initial event to start writing
        
    .. cpp:member:: void** buffer
    
        Input buffer with data
        
    .. cpp:member:: gsl_matrix** matrix
    
        Output GSL matrix 
        
    .. cpp:member:: long numCol
    
        Number of columns
        
    .. cpp:member:: int numRow
    
        Number of rows
        
    .. cpp:member:: int type
    
        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
    .. cpp:member:: int eventini
    
        Initial event to start writing
        
    
.. cpp:function:: int toGslVector(void **buffer, gsl_vector **array, long nevent, int eventini, int type)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input buffer into an output GSL vector.
    
    **Members/Variables**
    
    void** **buffer**
    
        Input buffer with data
        
    gsl_vector** **array**
    
        Output GSL vector
        
    long **nevent**
    
        Number of elements to store
        
    int **eventini**
    
        Initial element number
        
    int **type**

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
        
    .. cpp:member:: void** buffer
    
        Input buffer with data
        
    .. cpp:member:: gsl_vector** array
    
        Output GSL vector
        
    .. cpp:member:: long nevent
    
        Number of elements to store
        
    .. cpp:member:: int eventini
    
        Initial element number
        
    .. cpp:member:: int type

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
        
    
.. _U:

.. _V:

.. cpp:function:: int vector2matrix(gsl_vector *vectorin, gsl_matrix **matrixout)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function converts an input :math:`n^2` vector into an output square matrix :math:`[n \times n]`. It puts the first :math:`n` elements of the vector in the first row of the matrix, the second group of :math:`n` elements (from :math:`n` to :math:`2n-1`) of the vector in the second row and so on.
    
    **Members/Variables**
    
    gsl_vector** **vectorin**

        GSL input vector whose length is :math:`n^2`
        
    gsl_matrix* **matrixout**

        GSL output square matrix :math:`[n \times n]`
        
    .. cpp:member:: gsl_vector** vectorin

        GSL input vector whose length is :math:`n^2`
        
    .. cpp:member:: gsl_matrix* matrixout

        GSL output square matrix :math:`[n \times n]`
        
    
.. _W:

.. cpp:function:: int weightMatrix(ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, long nonpileupPulses, gsl_vector *nonpileup, gsl_vector *pulseaverage, gsl_matrix **covariance, gsl_matrix **weight)    
        
    Located in file: *tasksSIRENA.cpp*
    
    :cite:`Fixsen2004`
    
    This function calculates the weight matrix by using the non piled-up pulses found in all the records, stored in *pulsesAll* (previous records) and *pulsesInRecord* (current record). The weight matrix of each energy (and other intermediate values) will be stored in the library by the function :cpp:func:`fillInLibraryData`.

    Definitions:

        :math:`S_i^p`: Value of the ith-sample of the pulse number *p*
        
        :math:`M_i^p`: Value of the ith-sample of the model number *p* (model= *pulseaverage*):
        
        .. math::

            M_i = <S_i> = (1/N)\sum_{p=1}^{N}S_i^p
            
        N: number of non piled-up pulses
        
        .. math::

            & D_i = S_i - M_i \\
            & V_{ij} = <D_iD_j> = E[(S_i-M_i)(S_j-M_j)] = (1/N)\sum_{p=1}^{N}(S_i^p-M_i^p)(S_j^p-M_j^p) \\
            & V = \left[\begin{matrix} <D_1D_1> & <D_1D_2> & ... & <D_1D_n> \\
            <D_2D_1> & <D_2D_2> & ... & <D_2D_n> \\
            ....  &  ....  & ... &  ....  \\
            <D_nD_1> & <D_nD_2> & ... & <D_nD_n>\end{matrix}\right]

        where *n* is the :option:`PulseLength` and thus :math:`V = [n \times n]`.
        
        The weight matrix :math:`W = [V]^{-1}`.

    Steps:

        - Calculate the elements of the diagonal of the covariance matrix
        - Calculate the elements out of the diagonal of the covariance matrix
        - If saturated pulses :math:`\Rightarrow` Covariance matrix is a singular matrix :math:`\Rightarrow` Non invertible 
        
          In order to allow the covariance matrix to be inverted :math:`\Rightarrow` Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
          
          - Elements of the diagonal: Generating a random double :math:`f_1` between a range *(fMin,fMax)* (-NoiseStd,NoiseStd) to replace 0's with :math:`f_1^2`
          - Elements out of the diagonal: Generating two random doubles :math:`f_1` and :math:`f_2` between a range *(fMin,fMax)* (-NoiseStd,NoiseStd) to replace 0's with :math:`f_1 \cdot f_2`
     
        - Calculate the weight matrix

    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    bool **saturatedPulses**

        If *true*, all the pulses ( :option:`opmode` = 0 :math:`\Rightarrow` all the pulses have the same energy) are saturated
        
    PulsesCollection* **pulsesAll** 

        Collection of pulses found in the previous records
        
    PulsesCollection* **pulsesInRecord**

        Collection of pulses found in the current record
        
    long **nonpileupPulses**

        Number of non piled-up pulses
        
    gsl_vector* **nonpileup**

        GSL vector containing info about all the pulses informing if they are piled-up or not
        
    gsl_vector** **pulseaverage**

        GSL vector with the pulseaverage (= template = model) of the non piled-up pulses

    gsl_matrix** **covariance**

        GSL matrix with covariance matrix
        
    gsl_matrix** **weight**

        GSL matrix with weight matrix
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    .. cpp:member:: bool saturatedPulses

        If *true*, all the pulses ( :option:`opmode` = 0 :math:`\Rightarrow` all the pulses have the same energy) are saturated
        
    .. cpp:member:: PulsesCollection* pulsesAll 

        Collection of pulses found in the previous records
        
    .. cpp:member:: PulsesCollection* pulsesInRecord

        Collection of pulses found in the current record
        
    .. cpp:member:: long nonpileupPulses

        Number of non piled-up pulses
        
    .. cpp:member:: gsl_vector* nonpileup

        GSL vector containing info about all the pulses informing if they are piled-up or not
        
    .. cpp:member:: gsl_vector** pulseaverage

        GSL vector with the pulseaverage (= template = model) of the non piled-up pulses

    .. cpp:member:: gsl_matrix** covariance

        GSL matrix with covariance matrix
        
    .. cpp:member:: gsl_matrix** weight

        GSL matrix with weight matrix
    
    
.. cpp:function:: int weightMatrixNoise (gsl_matrix *intervalMatrix, gsl_matrix **weight)

    Located in file: *gennoisespec.cpp*
    
    This function calculates the weight matrix of the noise

        :math:`D_i`: Pulse free interval
        :math:`V`: Covariance matrix

            :math:`V_{ij} = E[DiDj]-E[Di]E[Dj]` 

        :math:`Di^p`: Value of the pth-sample of the pulse-free interval i
        :math:`N`: Number of samples
    
        .. math::

            & V_{ij} =  <D_iD_j> = E[D_iD_j] = (1/N)sum_{p=1}^{N}(Di^p)(Dj^p) \\
            & V = \left[\begin{matrix} <D_1D_1> & <D_1D_2> & ... & <D_1D_n> \\
            <D_2D_1> & <D_2D_2> & ... & <D_2D_n> \\
            ....  &  ....  & ... &  ....  \\
            <D_nD_1> & <D_nD_2> & ... & <D_nD_n>\end{matrix}\right]
            
        where *n* is the :option:`PulseLength` and thus :math:`V = [n \times n]`.
        
        The weight matrix :math:`W = 1/V`.
        
    Steps:

    - Calculate the elements of the diagonal of the covariance matrix
    - Calculate the elements out of the diagonal of the covariance matrix
    - Calculate the weight matrix

    **Members/Variables**
    
    gsl_matrix* **intervalMatrix**
    
        GSL matrix containing pulse-free intervals whose baseline is 0 (baseline previously subtracted) [nintervals x intervalMinSamples]
        
    gsl_matrix** **weight** 
    
        GSL matrix with weight matrix
        
    .. cpp:member:: gsl_matrix* intervalMatrix
    
        GSL matrix containing pulse-free intervals whose baseline is 0 (baseline previously subtracted) [nintervals x intervalMinSamples]
        
    .. cpp:member:: gsl_matrix** weight 
    
        GSL matrix with weight matrix
    
    
.. cpp:function:: int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double energy, gsl_vector *optimalfilter, fitsfile **dtcObject)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function runs in RECONSTRUCTION mode and writes the optimal filter info (in the *FILTER* HDU) for each pulse
    if :option:`intermediate` = 1 and either :option:`OFLib` = no or :option:`OFLib` = yes, :option:`filtEeV` = 0 and the the number of energies in the library FITS file is greater than 1.

    - Declare variables
    - Open intermediate FITS file
    - If (:option:`OFLib` = no) or (option:`OFLib` = yes, :option:`filtEeV` = 0 and the the number of energies in the library FITS file is greater than 1):
        - Create the *FILTER* HDU if it is the first pulse
        - Write data:
            - **OPTIMALF** or **OPTIMALFF** column (in time or frequency domain)
            - **OFLENGTH** column
    - Write **ENERGY** column in *PULSES* HDU
    - Close intermediate output FITS file if it is necessary
    - Free memory

    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    int **pulse_index**
    
        Index of the pulse whose info is going to be written (to know if it is the first pulse)
        
    double **energy**
    
        Estimated energy (eV)
        
    gsl_vector* **optimalfilter** 
    
        Optimal filter (in time or frequency domain)
        
    fitsfile** **dtcObject** 
    
        Fitsfile object for intermeadiate file name
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    .. cpp:member:: int pulse_index
    
        Index of the pulse whose info is going to be written (to know if it is the first pulse)
        
    .. cpp:member:: double energy
    
        Estimated energy (eV)
        
    .. cpp:member:: gsl_vector* optimalfilter 
    
        Optimal filter (in time or frequency domain)
        
    .. cpp:member:: fitsfile** dtcObject 
    
        Fitsfile object for intermeadiate file name
    
    
.. cpp:function:: int writeFitsComplex(IOData obj, gsl_matrix *matrix)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a GSL matrix. After that, the function puts them into a complex column of the output FITS file.
    
    **Members/Variables**
    
    IOData **obj**
    
        Object for FITS column to be written
        
    gsl_matrix* **matrix**
    
        Input GSL matrix with data
        
    .. cpp:member:: IOData obj
    
        Object for FITS column to be written
        
    .. cpp:member:: gsl_matrix* matrix
    
        Input GSL matrix with data
    
    
.. cpp:function:: int writeFitsSimple(IOData obj, gsl_vector *vector)
    
    Located  in file: *inoututils.cpp*
    
    This function reads values of a GSL vector.	After that, the function puts them into a column of the output FITS file.
    
    **Members/Variables** 
    
    IOData **obj**
        
        Object for FITS column to be written
        
    gsl_vector* **vector**
        
        Input GSL vector with data
        
    .. cpp:member:: IOData obj
        
        Object for FITS column to be written
        
    .. cpp:member:: gsl_vector* vector
        
        Input GSL vector with data
        
        
.. cpp:function:: int writeLibrary(ReconstructInitSIRENA *reconstruct_init, double samprate, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject, gsl_vector *pulsetemplateMaxLengthFixedFilter)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function writes the library (reordering if it is necesary and calculating some intermediate parameters)

        - Adding a new row to the library if *appendToLibrary = true* (:cpp:func:`readAddSortParams`)
        - Write the first row of the library if *appendToLibrary = false* (:cpp:func:`addFirstRow`)  
        
        - In both cases, the keywords ``CREADATE`` and ``SIRENAV`` with the date and SIRENA version are written 

    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values) 
        
    double **samprate**

        Sampling rate
        
    double **estenergy**

        Pulse height of the template whose energy is going to be added to the library
        
    gsl_vector* **pulsetemplate**

        GSL vector with the pulse template whose energy is going to be added to the library
    
    gsl_matrix** **covariance**

        GSL matrix with covariance matrix
        
    gsl_matrix** **weight**

        GSL matrix with weight matrix
        
    bool **appendToLibrary**

        *true* if adding a new row to the library and *false* if it is the first row to be added
        
    fitsfile** **inLibObject**

        FITS object containing information of the library FITS file 

    gsl_vector* **pulsetemplateMaxLengthFixedFilter**

        GSL vector with the :option:`largeFilter`-length pulse template whose energy is going to be added to the library
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values) 
        
    .. cpp:member:: double samprate

        Sampling rate
        
    .. cpp:member:: double estenergy

        Pulse height of the template whose energy is going to be added to the library
        
    .. cpp:member:: gsl_vector* pulsetemplate

        GSL vector with the pulse template whose energy is going to be added to the library
    
    .. cpp:member:: gsl_matrix** covariance

        GSL matrix with covariance matrix
        
    .. cpp:member:: gsl_matrix** weight

        GSL matrix with weight matrix
        
    .. cpp:member:: bool appendToLibrary

        *true* if adding a new row to the library and *false* if it is the first row to be added
        
    .. cpp:member:: fitsfile** inLibObject

        FITS object containing information of the library FITS file 
	
    .. cpp:member:: gsl_vector* pulsetemplateMaxLengthFixedFilter

        GSL vector with the :option:`largeFilter`-length pulse template whose energy is going to be added to the library
    
    
.. cpp:function:: void writeLog(FILE *fileRef, string type, int verbosity, string message)

    Located  in file: *inoututils.cpp*
    
    This function includes the processing of the each level of message in the log file and the output screen:
     
      - Verbosity = 0 :math:`\Rightarrow` The log file and the output screen include Errors
      - Verbosity = 1 :math:`\Rightarrow` The log file and the output screen include Errors and Warnings
      - Verbosity = 2 :math:`\Rightarrow` The log file and the output screen include Errors, Warnings and Alerts
      - Verbosity = 3 :math:`\Rightarrow` The log file and the output screen include Errors, Warnings, Alerts and Log messages

    **Members/Variables**
    
    FILE* **fileRef**
    
        File reference to log file 
        
    string **type**
    
        String to indicate error type "Error", "Warning", "Alert","Log" or "OK" 

     int **verbosity**
        
        Integer value for verbosity
        
     string **message**
        
        String message to print 
        
    .. cpp:member:: FILE* fileRef
    
        File reference to log file 
        
    .. cpp:member:: string type
    
        String to indicate error type "Error", "Warning", "Alert","Log" or "OK" 

    .. cpp:member::  int verbosity
        
        Integer value for verbosity
        
    .. cpp:member::  string message
        
        String message to print 
    

.. cpp:function:: int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject)
        
    Located in file: *tasksSIRENA.cpp*

    This function writes the data of the pulses found in the record in the intermediate FITS file (in the *PULSES* HDU). The pulses info given is: **TSTART**, **I0** (the pulse itself), **TEND**, **TAURISE**, **TAUFALL** and **QUALITY**.
    
    **Members/Variables**
    
    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    double **samprate** 

        Sampling rate (to convert samples to seconds)

    double **initialtime**

        Starting time of the record (in order to calculate absolute times)

    gsl_vector* **invectorNOTFIL**

        GSL vector with the original record (neither low-pass filtered nor differentiated)

    int **numPulsesRecord**

        Number of pulses found in the record

    gsl_vector* **tstart**

        GSL vector with the start times of the found pulses

    gsl_vector* **tend**

        GSL vector with the end times of the found pulses
    
    gsl_vector* **quality** 

        GSL vector with the quality of the found pulses
        
        0 :math:`\Rightarrow` Standard (good) pulses
        
        1 :math:`\Rightarrow` Truncated pulses at the beginning
        
        2 :math:`\Rightarrow` Truncated pulses at the end
        
        10 :math:`\Rightarrow` Saturated pulses
        
        11 :math:`\Rightarrow` Truncated and saturated pulses

    gsl_vector* **taurise**

        GSL vector with the rise time constants of the found pulses (to be done)

    gsl_vector* **taufall**

        GSL vector with the fall time constants of the found pulses (to be done)

    fitsfile* **dtcObject**

        Object which contains information of the intermediate FITS file
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
        
    .. cpp:member:: double samprate 

        Sampling rate (to convert samples to seconds)

    .. cpp:member:: double initialtime 

        Starting time of the record (in order to calculate absolute times)

    .. cpp:member:: gsl_vector* invectorNOTFIL 

        GSL vector with the original record (neither low-pass filtered nor differentiated)

    .. cpp:member:: int numPulsesRecord 

        Number of pulses found in the record

    .. cpp:member:: gsl_vector* tstart

        GSL vector with the start times of the found pulses

    .. cpp:member:: gsl_vector* tend

        GSL vector with the end times of the found pulses
    
    .. cpp:member:: gsl_vector* quality 

        GSL vector with the quality of the found pulses
        
        0 :math:`\Rightarrow` Standard (good) pulses
        
        1 :math:`\Rightarrow` Truncated pulses at the beginning
        
        2 :math:`\Rightarrow` Truncated pulses at the end
        
        10 :math:`\Rightarrow` Saturated pulses
        
        11 :math:`\Rightarrow` Truncated and saturated pulses

    .. cpp:member:: gsl_vector* taurise

        GSL vector with the rise time constants of the found pulses (to be done)

    .. cpp:member:: gsl_vector* taufall

        GSL vector with the fall time constants of the found pulses (to be done)

    .. cpp:member:: fitsfile* dtcObject

        Object which contains information of the intermediate FITS file

    
.. cpp:function:: int writeTestInfo(ReconstructInitSIRENA* reconstruct_init, gsl_vector *recordDERIVATIVE, double threshold, fitsfile *dtcObject)    
    
    Located in file: *tasksSIRENA.cpp*

    This function writes the *TESTINFO* HDU in the intermediate FITS file. The written columns are **FILDER** (low-pass filtered and differentiated record) and **THRESHOLD**.

    **Members/Variables**

    ReconstructInitSIRENA** **reconstruct_init**

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    gsl_vector* **recordDERIVATIVE** 

        GSL vector with input record (low-pass filtered and) differentiated
        
    double **threshold**

        Threshold value used to find pulses
        
    fitsfile **dtcObject**

        Object which contains information of the intermediate FITS file
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
        
    .. cpp:member:: gsl_vector* recordDERIVATIVE 

        GSL vector with input record (low-pass filtered and) differentiated
        
    .. cpp:member:: double threshold

        Threshold value used to find pulses
        
    .. cpp:member:: fitsfile dtcObject

        Object which contains information of the intermediate FITS file
   
   
.. cpp:function:: int writeTPSreprExten ()

    Located in file: *gennoisespec.cpp*
    
    This function writes the noisespec output FITS file.

    Steps:
    
    - Allocate GSL vectors
    - Write the data in the output FITS file (print only half of FFT to prevent aliasing)
    - *NOISE* HDU only contains positive frequencies (:math:`\Rightarrow` Multiply by 2 the amplitude)
    - *NOISEALL* HDU contains negative and positive frequencies :math:`\Rightarrow` It is the HDU read to build the optimal filters
    - *WEIGHTMS* HDU
   
   
.. _X:

.. _Y:

.. _Z:



    
    
    
.. |mem| replace:: **Members/Variables**    
    
    
    
    
    
    
