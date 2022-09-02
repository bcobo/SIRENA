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

.. cpp:function:: int addFirstRow(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int runF0orB0val, gsl_vector *E, gsl_vector *PHEIGHT, gsl_matrix *PULSE, gsl_matrix *PULSEB0, gsl_matrix *MF, gsl_matrix *MFB0, gsl_matrix *COVAR, gsl_matrix *WEIGHT)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function writes the first row of the library (without intermediate AB-related values, because it would be necessary to have at least two rows=energies in the library). It also writes the **FIXFILTT** and **FIXFILTF** HDUs with the optimal filters in the time and frequency domain (real and imag parts) with fixed legnth.

    **Members/Variables**

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses 
        :option:`mode` and the noise_spectrum in order to run :cpp:func:`calculus_optimalFilter`.
     
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
    
    .. cpp:member:: double samprate

        Sampling rate

    .. cpp:member:: gsl_vector** vector1

        GSL vector with input vector

    .. cpp:member:: gsl_vector** vector2 

        GSL with input vector which is delayed or moved forward to be aligned with :cpp:member:`vector1`
            
.. _B:

.. _C:

.. cpp:function:: int calculateEnergy(gsl_vector *vector, int pulseGrade, gsl_vector *filter, gsl_vector_complex *filterFFT, int runEMethod, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, int domain, double samprate, gsl_vector *Pab, gsl_matrix *PRCLWN, double *calculatedEnergy)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the energy of a pulse (:cpp:member:`vector`) depending on the :option:`EnergyMethod` and the :option:`FilterDomain` 
    selected from input parameters.

    a) **OPTFILT** (= **I2R** or **I2RALL**, **I2RNOL** or **I2RFITTED**): Optimal filter = Wiener filter  (see :ref:`optimalFilter`)

    Once the filter template has been created (:cpp:member:`filter` or :cpp:member:`filterFFT`), pulse height analysis is performed by aligning the template
    with a pulse and multiplying each point in the template by the corresponding point in the pulse. The sum of these products is the energy.

    In the practice, the alignment of the pulse relative to the trigger is not completely accurate, so a number of *n* lags could be used in order to find 
    the peak value of the energy. The *n* peak values are fitted to a parabola to find the most accurate energy (:option:`LagsOrNot`).

    b) **WEIGHT** and **WEIGHTN** (see :ref:`covMatrices`)


    **Members/Variables**
    
    .. cpp:member:: gsl_vector* vector
    
        Pulse whose energy has to be determined
        
    .. cpp:member:: int pulseGrade 
    
        Grade of the input pulse (to decide whether a full or only a rough estimation of energy is required). If :cpp:member:`pulseGrade` =-1 (rejected) the provided calculated energy is -1
        
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

    .. cpp:member:: double nrmfctr 
    
        Normalization factor to make the estimation of *E* in units of energy
        
    .. cpp:member:: double samprate
    
        Sampling rate in Hz
        
    .. cpp:member:: gsl_vector* Pab 
    
        **PAB** column in the library 
        
    .. cpp:member:: gsl_vector* PRCLWN 
    
        **PRCLx** column in the library
        
    .. cpp:member:: double* calculatedEnergy
    
        Calculated energy in eV. 
        

.. cpp:function:: int calculateIntParams(ReconstructInitSIRENA *reconstruct_init, int indexa, int indexb, double samprate, int runF0orB0val, gsl_matrix *modelsaux, gsl_matrix *covarianceaux, gsl_matrix *weightaux, gsl_vector *energycolumn, gsl_matrix **Wabaux, gsl_matrix **TVaux, gsl_vector **tEcolumn, gsl_matrix **XMaux, gsl_matrix **YVaux, gsl_matrix **ZVaux, gsl_vector **rEcolumn, gsl_matrix **Pabaux, gsl_matrix **Dabaux, gsl_matrix **PrecalWMaux, gsl_matrix **optimalfiltersabFREQaux, gsl_matrix **optimalfiltersabTIMEaux)
    
    Located in file: *tasksSIRENA.cpp*

    This function calculates some intermediate scalars, vectors and matrices (WAB, TV, tE, XM, YV, ZV, rE, PAB and DAB) for the interpolation and covariance methods. See :ref:`covMatrices` reconstruction method.

    - Declare variables and allocate GSL vectors and matrices
    - Calculate intermediate scalars, vectors and matrices 
    - Free allocated GSL vectors and matrices
    
    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function
        uses :option:`PulseLength`
    
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
    
        GSL input matrix with weight matrix
        
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
    

.. cpp:function:: int calculateTemplate(ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight)    
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the template (**PULSE** column in the library) of non piled-up pulses.
    Just in case in the detection process some piled-up pulses have not been distinguished as different pulses, a pulseheights histogram is built. It uses the pulseheights histogram (built by using the **PHEIGHT** column of the library), **Tstart** and **quality** to select the non piled-up pulses.

    1) Declare and initialize variables
    
    2) Before building the histogram, select the pulseheights of the pulses well separated from other pulses and with *quality* = 0
    
    3) Create the pulseheights histogram
    
    4) Calculate the pulseaverage only taking into account the valid pulses:
    
        * Check if the pulse is piled-up or not
        
        * Non piled-up pulses :math:`\Rightarrow` Align and average them (currently :cpp:func:`align` is not used because the starting times of the pulses are provided as input parameter)
    
    5) Calculate covariance and weight matrices
    
    6) Free allocated GSL vectors
    
    **Members/Variables**

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  In particular, this function uses :option:`PulseLength` and :option:`EnergyMethod`
        
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
        
        
.. cpp:function::  int calculus_optimalFilter(int TorF, int intermediate, int mode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex)
    
    Located in file: *tasksSIRENA.cpp*
    
    See description also at :ref:`optimal filter chapter <optimalFilter>`
    
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
    
    - FFT calculus of the matched filter (normalized template)
    
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
    - Inverse FFT (to get the expression of the optimal filter in time domain)
    
        - Complex :math:`OptimalFilter(f)` :math:`\Rightarrow` Taking into account magnitude :math:`MatchedFilter(f)/N^2(f)` and phase given by :math:`MatchedFilter(f)`
    - Free allocated GSL vectors
        
    **Members/Variables**
    
    .. cpp:member:: int TorF

        If :option:`FilterDomain` = **T** :math:`\Rightarrow` :cpp:member:`TorF` = 0; If :option:`FilterDomain` = **F** :math:`\Rightarrow` :cpp:member:`TorF` = 1
        
    .. cpp:member:: int intermediate 

        If :option:`intermediate` = 0 :math:`\Rightarrow` Do not write an intermediate file; If :option:`intermediate` = 1 :math:`\Rightarrow` Write an intermediate file
        
    .. cpp:member:: int mode

        If :option:`mode` = 0 :math:`\Rightarrow` CALIBRATION run (library creation); If :option:`mode` = 1 :math:`\Rightarrow` RECONSTRUCTION run (energy determination)

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
        
.. cpp:function:: int convertI2R(ReconstructInitSIRENA** reconstruct_init, TesRecord **record, gsl_vector **invector)
    
    Located in file: *tasksSIRENA.cpp*
    
    This funcion converts the current space into a quasi-resistance space (see :ref:`rSpace` for **I2R**, **I2RALL**, **I2RNOL** and **I2RFITTED** modes). The input :cpp:member:`invector` filled in with current values is filled in here with *I2R*, *I2RALL*, *I2RNOL* or **I2RFITTED** quasi-resistances at the output.

    - Read ``R0``, ``I0_START``, ``ADUCNV`` and ``IMIN`` from the :option:`RecordFile` by using the pointer *(*reconstruct_init)->record_file_fptr*
    - Conversion according to :option:`EnergyMethod` = **I2R**
    - Conversion according to :option:`EnergyMethod` = **I2RALL**
        
      - Read also ``RL``, ``LFILTER`` and ``TTR`` from the :option:`RecordFile`
    - Conversion according to :option:`EnergyMethod` = **I2RNOL**
    
      - *I2RALL* neglecting the circuit inductance
    - Conversion according to :option:`EnergyMethod` = **I2RFITTED**  
    - Because in :cpp:func:`runEnergy` the record (*TesRecord* type) is used :math:`\Rightarrow` The **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED** transformed record has to be used

    **Members/Variables**

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: TesRecord** record

        Structure containing the input ADC record
        
    .. cpp:member:: gsl_vector** invector

        Input current (ADC) vector & output resistance (*I2R* or *I2RALL* or *I2RNOL* or *I2RFITTED*) vector
    

.. cpp:function:: int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, fitsfile **dtcObject)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function creates an intermediate FITS file with some useful info (during the development phase) if the :option:`intermediate` input parameter is set to 1.

    The intermediate FITS file will contain 2 HDUs:
    
        * **PULSES** HDU will contain some info about the found pulses: **TSTART**, **I0** (the pulse itself), **TEND**, **TAURISE**, **TAUFALL** and **QUALITY**
        
        * **TESTINFO** HDU will contain columns **FILDER** (the low-pass filtered and differentiated records) and **THRESHOLD**

    If file exists :math:`\Rightarrow` Check :option:`clobber` for overwritting. If it does not, then create it.
        
    **Members/Variables**
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses some of its
        values to call the intermediate file (:option:`detectFile`), to write input parameters info in the Primary extension of the intermediate file and also :option:`clobber`
                
    .. cpp:member:: double samprate
        
        Sampling rate (calculated from input FITS keyword ``DELTAT``)
            
    .. cpp:member:: fitsfile dtcObject
        
        Object which contains information of the intermediate FITS file (used also by :cpp:func:`writeTestInfo` and :cpp:func:`writePulses`).
                       

.. cpp:function:: int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
        
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
    
    .. cpp:member:: gsl_vector* invector

        GSL input vector

    .. cpp:member:: int nbins

        Number of bins to build the histogram

    .. cpp:member:: gsl_vector** xhistogsl

        GSL vector with output histogram x-axis

    .. cpp:member:: gsl_vector** yhistogsl

        GSL vector with output histogram y-axis
            
                            
.. cpp:function:: int createLibrary(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject)
    
    Located in file: *tasksSIRENA.cpp*

    This function creates the pulse templates library FITS file, if it does not exist yet. Otherwise, it opens it (to add a new row).

        1) If it exists :math:`\Rightarrow` Open it and set *appendToLibrary = true*
        
        2) If it does not exist :math:`\Rightarrow` Create it and set *appendToLibrary = false*

            - Write keyword ``EVENTCNT`` = 1 in the LIBRARY extension
            - Write info keywords about the running parameters in the Primary extension
            
    **Members/Variables**
            
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses parameters 
        to call the library (:option:`LibraryFile`) and to write input parameters info in the Primary extension of the library file
    
    .. cpp:member:: bool appendToLibrary

        Used by the function :cpp:func:`writeLibrary`
        
    .. cpp:member:: fitsfile** inLibObject

        Object which contains information of the library FITS file (used also by :cpp:func:`writeLibrary`)
        
.. _D:

.. cpp:function:: int differentiate(gsl_vector **invector, int szVct)

    Located in file: *pulseprocess.cpp*

    This function applies the derivative method :math:`x_i-x_{i-1}` to the input vector.

    The derivative method provides more sensitivity to handle with piled-up pulses.
    Moreover, little variations of the baseline will not affect.

    **Members/Variables**
    
    .. cpp:member:: gsl_vector** invector 
    
        Input/Ouput GSL vector (non-differentiate input vector/differentiate input vector)
        
    .. cpp:member:: int szVct
    
        Size of :cpp:member:`invector`

.. _E:

.. cpp:function:: int eigenVV (gsl_matrix *matrixin, gsl_matrix **eigenvectors, gsl_vector **eigenvalues)
    
    Located in file: *tasksSIRENA.cpp*
    
    This funcion provides the principal eigenvectors and eigenvalues of the input matrix (at the moment, the first two eigenvalues and eigenvectors).

    - Calculate the eigenvectors and the eigenvalues
    - Sort the eigenvectors and the eigenvalues in descending order
    - Choose the main eigenvectors and eigenvalues (the principal components analysis). At the moment, the first two eigenvectors and eigenvalues

    **Members/Variables**

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
    
    .. cpp:member:: gsl_vector *invector
    
        Input GSL vector
        
    .. cpp:member:: gsl_vector_complex *outvector
    
        Output GSL complex vector with the FFT of :cpp:member:`invector`

    .. cpp:member::  double STD
        
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*
      
      
.. cpp:function:: int FFTinverse(gsl_vector_complex *invector, gsl_vector *outvector, double STD)
    
    Located in file: *genutils.cpp*
    
    This function calculates the inverse FFT of the elements of a vector.
    
    **Members/Variables**
    
    .. cpp:member:: gsl_vector_complex *invector
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector *outvector
    
        Output GSL vector with the inverse FFT of :cpp:member:`invector`
    
    .. cpp:member:: double STD
    
        SelectedTimeDuration = (Size of :cpp:member:`invector`)/*samprate*
        
        
.. cpp:function:: int filderLibrary (ReconstructInitSIRENA** reconstruct_init, double samprate)
        
    Located in file: *tasksSIRENA.cpp*

    This function calculates the (low-pass filter and) derivative of the models (*pulse_templates*) in the library (only necessary if first record), 
    and it stores the *pulse_templates_filder* and the *maxDERs* (see below) in the :cpp:member:`reconstruct_init` structure.

    The maximum of the (low-pass filtered and) differentiated pulse has to be compared to the *maxDERs* to select the appropriate model.

    1) Check if it is the first record
    
    2) (Low-pass filter and) differentiate the models (*pulse_templates*) of the library
    
    3) Store the (low-pass filtered) derivatives in *pulse_templates_filder*
    
    4) Calculate the maximum of the (low-pass filtered and) differentiated models (*maxDERs*)

    **Members/Variables**

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses some parameters 
        to filter (:option:`scaleFactor`) and other parameters to handle the pulse templates and their derivatives (*ntemplates*, *pulse_templates*, 
        *pulse_templates_filder* and *maxDERs*)

    .. cpp:member:: double samprate

        Sampling rate
            
            
.. cpp:function:: bool fileExists(const std::string& name)
    
    Located in file: *genutils.cpp*
    
    This function checks for file existence returning a boolean value.
    
    **Members/Variables**
    
    .. cpp:member:: const std::string& name
         
        File name 
        
        
.. cpp:function:: int filterByWavelets (ReconstructInitSIRENA* reconstruct_init, gsl_vector **invector, int length, int *onlyOnce)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function filters the input/output signal :cpp:member:`invector`, reducing the noise level.
    
    Steps:
    
    - It is only going to work with *n* elements of :cpp:member:'invector'
    - Discrete Wavelet Treansform 
    - Sorting coefficients
    - Hard thresholding: *n-nc* coefficeints are deleted (those with low energy)
    - Inverse DWT
        
    **Members/Variables**
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). 
                
    .. cpp:member:: gsl_vector **invector
        
        Input/output signal 

    .. cpp:member:: int length
        
        Length of the wavelet transform
        
    .. cpp:member:: int *onlyOnce
        
        In order to control the times to be executed
        
        
.. cpp:function:: int findMeanSigma(gsl_vector *invector, double *mean, double *sigma)
    
    Located in file: *pulseprocess.cpp*
    
    This function calculates the mean and the standard deviation of the input vector.
    
    **Members/Variables**
    
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
    - First step to look for single pulses (call :cpp:func:`medianKappaClipping`)
    
      - Call :cpp:func:`findTstartCAL`
    - If at least a pulse is found
    
      - Get :cpp:member:`pulseheight` of each found pulse (in order to be used to build the pulse templates library)
    - Free allocated GSL vectors

    **Members/Variables**
    
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
        
    .. cpp:member:: double*threshold
    
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
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses 
        *maxPulsesPerRecord* (:option:`EventListSize`) and :option:`PulseLength` (and :cpp:func:`findTstartCAL` uses :option:`tstartPulse1`, :option:`tstartPulse2` and :option:`tstartPulse3`)
        
    .. cpp:member:: double stopcriteriamkc
    
        Used in :cpp:func:`medianKappaClipping` (%)

    .. cpp:member:: double kappamkc
        
        Used in :cpp:func:`medianKappaClipping`
        
  
.. cpp:function:: int FindSecondaries(int maxPulsesPerRecord, gsl_vector *adjustedDerivative, double adaptativethreshold, ReconstructInitSIRENA *reconstruct_init, int tstartFirstEvent, int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl)

    Located in file: *pulseprocess.cpp*

    This function runs after :cpp:func:`InitialTriggering` to find all the events (except the first one) in the first derivative of the (low-pass filtered) record.

    Steps: 
    
    - Declare variables
    - It is possible to find the tstarts... 
    
      It looks for an event and if a pulse is found (foundPulse=true), it looks for another event
    
        - It looks for an event since the beginning (or the previous event) to the end of the record. 
          The first condition to detect an event is that the :cpp:member:`adjustedDerivative` was over the :cpp:member:`threshold`

        - Subtract the model from the adjusted derivative

            - Select the model of the found event from the libary by using the first sample of the derivative 
            - Subtract
    
    - ... Or to use the tstart provided as input parameters
    
        - Obtain the *maxDERs* of the events whose tstarts have been provided (by using the maximum of the derivative to find the model)
    - Free allocated GSL vectors

    **Members/Variables**
    
    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of events per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    .. cpp:member:: gsl_vector* adjustedDerivative
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double adaptativethreshold
        
        Threshold
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
        
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses some parameters
        (:option:`PulseLength`, :option:`tstartPulse1`, :option:`tstartPulse2`, :option:`tstartPulse3`, :option:`EnergyMethod`) and the templates
    
    .. cpp:member:: int tstartFirstEvent
        
        Tstart of the first event of the record (in samples) found by :cpp:func:`InitialTriggering`
    
    .. cpp:member:: int* numberPulses
        
        Number of found events
    
    .. cpp:member:: gsl_vector** tstartgsl
        
        Starting time of the found events (in samples)
    
    .. cpp:member:: gsl_vector** flagTruncated
        
        Flag indicating if the event is truncated (inside this function only initial truncated pulses are classified)
        
    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the first derivative of the (low-pass filtered) record inside each found event 
        
    .. cpp:member:: gsl_vector** samp1DERgsl
        
        First sample of the first derivative of the (low-pass filtered) record inside each found event
        
        
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
    
    .. cpp:member:: int maxPulsesPerRecord
    
        Expected maximum number of pulses per record in order to not allocate the GSL variables with the size of the input vector (:option:`EventListSize`)
        
    .. cpp:member:: gsl_vector* der
    
        First derivative of the (low-pass filtered) record

    .. cpp:member:: double adaptativethreshold
        
        Threshold
        
    .. cpp:member:: int nSamplesUp
    
        Number of consecutive samples over the threshold to 'find' a pulse (:option:`samplesUp`)
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function
        uses :option:`PulseLength`, :option:`tstartPulse1`, :option:`tstartPulse2` and :option:`tstartPulse3`
        
    .. cpp:member:: int* numberPulses
    
        Number of found pulses
        
    .. cpp:member:: gsl_vector** tstartgsl
    
        Pulses tstart (in samples)
        
    .. cpp:member:: gsl_vector** flagTruncated
    
        Flag indicating if the pulse is truncated 

    .. cpp:member:: gsl_vector** maxDERgsl
        
        Maximum of the first derivative of the (low-pass filtered) record inside each found pulse
             
        
.. cpp:function:: int find_Esboundary(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, int *indexEalpha, int *indexEbeta, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*.
    
    This function provides the indexes of the two energies which straddle the pulse energy, by  comparing the maximum value of the pulse derivative
    (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`).

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` :cpp:member:`indexEalpha` = :cpp:member:`indexEbeta` = 0
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` :cpp:member:`indexEalpha` = :cpp:member:`indexEbeta` = Number of templates-1 
 
    **Members/Variables**

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
    
    
.. cpp:function:: int find_matchedfilter(int runF0orB0val, double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function chooses the proper matched filter (normalized template) of the calibration library (**MF** or **MFB0**) by comparing the maximum of their derivatives (:cpp:member:`maxDER` versus :cpp:member:`maxDERs`).

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library and interpolates (:cpp:func:`interpolate_model`)
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The matched filter with the lowest :cpp:member:`maxDERs` (first row) in the library is chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The matched filter with the highest :cpp:member:`maxDERs` (last row) in the library is chosen

    **Members/Variables**
    
    .. cpp:member:: int runF0orB0val

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose matched filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses 
        the info in the library (*matched_filters*)
    
    .. cpp:member:: gsl_vector** matchedfilterFound
    
        GSL vector with the matched filter selected
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
        
.. cpp:function:: int find_matchedfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function selects the proper matched filter (normalized template) from the calibration library (from column **DAB**) by comparing 
    the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`) for the *DAB* interpolation method (see :ref:`optimal filter chapter <optimalFilter>`). It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**
    
    .. cpp:member:: int runF0orB0val

        If :option:`FilterMethod` = **F0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 1. If :option:`FilterMethod` = **B0** :math:`\Rightarrow` :cpp:member:`runF0orB0val` = 0

    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose matched filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info in the library (*matched_filters*)
    
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
    
    .. cpp:member:: double samp1DER
    
        1st sample of the derivative of the filtered pulse whose pulse template is being sought
        
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the 1st samples of the derivative of the models (*samp1DERs*) and their derived templates (*pulse_templates_filder*), the number of templates in the library (*ntemplates*) and the template duration (*template_duration*).

    .. cpp:member:: gsl_vector** modelFound
        
        Found template of the pulse whose 1st sample of the derivative of the filtered pulse is :cpp:member:`samp1DER`
        
        
.. cpp:function:: int find_optimalfilter(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function selects the proper optimal filter from the calibration library (column **OF**) by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`).

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library and interpolates (:cpp:func:`interpolate_model`) 
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The optimal filter with the lowest :cpp:member:`maxDERs` (first row) in the library is chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The optimal filter with the highest :cpp:member:`maxDERs` (last row) in the library is chosen

    **Members/Variables**

    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info 
        in the library (*optimal_filters*)
    
    .. cpp:member:: gsl_vector** optimalfilterFound
    
        GSL vector with the optimal filter selected
        
    .. cpp:member:: double* Ealpha 
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the lower limit
        
    .. cpp:member:: double* Ebeta
    
        Energy (in eV) which straddle the :cpp:member:`maxDER` in the higher limit
        
        
.. cpp:function:: int find_optimalfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function selects the proper optimal filter from the calibration library (column **OFAB**) by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library  (:cpp:member:`maxDERs`). It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**

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
        
        
.. cpp:function:: int find_prclwn(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **PRCLWNFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
    
    Located in file: *tasksSIRENA.cpp*
    
    When :option:`EnergyMethod` = **WEIGHTN** this function selects the proper precalculated values (**PRCLx**) from the calibration library by comparing the maximum value of the pulse derivative (:cpp:member:`maxDER`) to the list of maximums in the library (:cpp:member:`maxDERs`) for the :option:`OFLib`=yes. It also selects the proper row from the column **PAB**.

    It finds the two embracing :cpp:member:`maxDERs` in the calibration library
    
        - If :cpp:member:`maxDER` is lower than the lowest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data with the lowest :cpp:member:`maxDERs` (first row) in the library are chosen
    
        - If :cpp:member:`maxDER` is higher than the highest :cpp:member:`maxDERs` in the library :math:`\Rightarrow` The data of the penultimate row in the library are chosen

    **Members/Variables**

    .. cpp:member:: double maxDER
    
        Max value of the derivative of the (filtered) pulse whose optimal filter is being sought
        
    .. cpp:member:: gsl_vector* maxDERs
    
        GSL vector with the maximum values of the derivatives of the templates in the library to be compared with the pulse being analysed
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init 
    
        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the info in the library (*PRECALWN*)
    
    .. cpp:member:: gsl_vector** PRCLWNFound
    
        GSL vector with the optimal filter selected
        
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
    
    .. cpp:member:: OptimalFilterSIRENA* OFilterColl
    
        Instance of *OptimalFilterSIRENA* structure
    
    
.. cpp:function:: extern_C_void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init)    
    
    Located in file: *integraSIRENA.cpp*
    
    Destructor of *ReconstructInitSIRENA* structure.
    
    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init
    
        Instance of *ReconstructInitSIRENA* structure
    
    
.. cpp:function:: extern_C_void freePulsesCollection(PulsesCollection* PulsesColl)
    
    Located in file: *integraSIRENA.cpp*
    
    Destructor of *PulsesCollection* structure.
    
    **Members/Variables**
    
    .. cpp:member:: PulsesCollection* PulsesColl
        
        Instance of *PulsesCollection* structure
    
    
.. cpp:function:: int fromGslMatrix(void **buffer, gsl_matrix **matrix, int type)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input GSL matrix into an output buffer.
    
    **Members/Variables**
    
    .. cpp:member:: void** buffer
    
        Output buffer
    
    .. cpp::member:: gsl_matrix** matrix
    
        Input GSL matrix
        
    .. cpp::member:: int type

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
    
.. cpp:function:: int fromGslVector(void **buffer, gsl_vector **array, int type)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input GSL vector into an output buffer.
    
    **Members/Variables**
    
    .. cpp:member:: void** buffer
    
        Output buffer
        
    .. cpp:member:: gsl_vector** array
    
        Input GSL vector
        
    . cpp:member:: int type

        FITS type (TINT, TSHORT, TDOUBLE, etc.)
    
.. _G:

.. cpp:function:: int getB(gsl_vector *vectorin, gsl_vector *tstart, int nPulses, gsl_vector **lb, int sizepulse, gsl_vector **B)
    
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
    
    
.. cpp:function:: LibraryCollection* getLibraryCollection(const char* const filename,int mode, char *energy_method, char *filter_method, char oflib, char **ofinterp, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function creates and retrieves a *LibraryCollection* from a file.
    
    - Create *LibraryCollection* structure
    - Open FITS file in READONLY mode (move to the first HDU) and get number of templates (rows)
    - Allocate library structure
    - Get PULSE and MF column numbers (depending the different options)
    - Get template duration
    - Allocate library structure (cont.)
    - Get matched filter duration
    - Read different columns and populate the *LibraryCollection* structure
    - Added new code to handle the new HDUs **FIXFILTF**, **FIXFILTT** and **PRECALWN**
    - Free allocated GSL vectors and matrices
    
    **Members/Variables**
    
    .. cpp:member::const char* const filename
        
        File with library information
        
    .. cpp:member:: int mode
    
        Calibration run (0) or energy reconstruction run (1), :option:`mode`
    
    .. cpp:member:: char* energy_method
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED**, :option:`EnergyMethod`
        
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
    
    .. cpp :member:: char oflib
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
        
    .. cpp:member:: char** ofinterp
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in ':ref:`tesreconstruction` as *DAB* (but it would be possible to work with *MF*)
        
    .. cpp:member:: int* const status
    
        Input/output status
    
    
.. cpp:function:: NoiseSpec* getNoiseSpec(const char* const filename,int mode, char *energy_method, char *filter_method, int* const status)
    
    Located  in file: *integraSIRENA.cpp*
    
    This function creates and retrieves a *NoiseSpec* from a file.
    
    - Create *NoiseSpec* structure
    - Open FITS file, move to the NOISE and NOISEALL HDUs and get interesting keywords
    - Allocate *NoiseSpec* structure
    - Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
    - Read column CSD and save it into the structure
    - Read column FREQ and save it into the structure
    - Return noise spectrum
    
    **Members/Variables**
    
    .. cpp:member:: const char* const filename
    
        File name with noise
        
    .. cpp:member:: int mode
    
        Calibration run (0) or energy reconstruction run (1), :option:`mode`
        
    .. cpp:member:: char* energy_method
    
        Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED** or **PCA**, :option:`EnergyMethod`
        
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **B0** (deleting the baseline), :option:`FilterMethod`
        
    .. cpp:member:: int* const status
    
        Input/Output status
    
    
.. cpp:function:: int getPulseHeight(gsl_vector *vectorin, double tstart, double tstartnext, int lastPulse, double lrs, double lb, double B, int sizepulse, double *pulseheight)
    
    Located  in file: *pulseprocess.cpp*
    
    This function estimates the pulse height of a pulse by using a running sum filter. It extracts from the record, :cpp:member:`vectorin`, the pulse whose 
    pulse height is going to be estimated by using :cpp:func:`RS_filter`.

    Steps:
    
    - Declare variables
    - Extracting from the record the pulse whose pulse height is going to be estimated
    - Apply the running sum filter
    
    **Members/Variables**
    
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
        
    .. cpp:member:: double *pulseheight
    
        Estimated pulse height of the pulse
        
            
.. cpp:function:: void gsl_vector_complex_absIFCA(gsl_vector *cvnew,gsl_vector_complex *cv)
    
    Located in file: *genutils.cpp*
    
    This function calculates the magnitude of the complex elements of a vector (real part).
    
    **Members/Variables**
    
    .. cpp:member:: gsl_vector_complex *cv
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector *cvnew
    
        Output GSL vector with the absolute values of the elements of :cpp:member:`cv`
    
    
.. cpp:function:: void gsl_vector_complex_argIFCA(gsl_vector *varg, gsl_vector_complex *vin)
    
    Located in file: *genutils.cpp*

    This function calculates the arguments of the complex elements of a vector.
    
    **Members/Variables**
    
    .. cpp:member:: gsl_vector_complex *vin
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector* varg
    
        Output GSL vector with the arguments of the elements of :cpp:member:`vin`
    
    
.. cpp:function:: void gsl_vector_complex_scaleIFCA(gsl_vector_complex *cv,gsl_complex z)
    
    Located in file: *genutils.cpp*
    
    This function multiplies the complex elements of a vector by a complex number.
    
    **Members/Variables**
    
    .. cpp:member:: gsl_vector_complex *cv
    
        Input/Output (scaled) GSL complex vector
        
    .. cpp:member:: gsl_complex z
    
        Input GSL complex number
    
    
.. cpp:function:: void gsl_vector_sqrtIFCA(gsl_vector *cvnew, gsl_vector *cv)
    
    Located in file: *genutils.cpp*
    
    This function calculates the square root of the elements of a vector.
    
    **Members/Variables**

    .. cpp:member:: gsl_vector *cv
    
        Input GSL complex vector
        
    .. cpp:member:: gsl_vector *cvnew
    
        Output GSL vector with the square root values of the elements of :cpp:member:`cv`
    
    
.. cpp:function:: int gsl_vector_Sumsubvector(gsl_vector *invector, long offset, long n, double *sum)
    
    Located in file: *genutils.cpp*
    
    This function returns the sum of some elements of the input vector.
    
    The starting element of the sum is :cpp:member:`offset` from the start of the input vector. It will sum up :cpp:member:`n` elements.

    :cpp:member:`offset` can take values from 0 to *invector->size*

    **Members/Variables**
    
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

.. cpp:function:: extern_C_void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file, fitsfile *fptr, char* const library_file, char* const event_file, int pulse_length, double scaleFactor, double samplesUp, double nSgms, int mode, double LrsT, double LbT, char* const noise_file, char* filter_domain, char* filter_method, char* energy_method, int lagsornot, int ofiter, char oflib, char* ofinterp, char* oflength_strategy, int oflength, double monoenergy, int interm, char* const detectFile, char* const filterFile, char clobber, int maxPulsesPerRecord, double SaturationValue, int tstartPulse1, int tstartPulse2, int tstartPulse3, double energyPCA1, double energyPCA2, char * const XMLFile, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function initializes the structure *ReconstructInitSIRENA* with the variables required for SIRENA reconstruction. The values are taken from the input parameters.
    
    - Load *LibraryCollection* structure if library file exists
    - Load *NoiseSpec* structure
    - Fill in *reconstruct_init*
    
    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA* reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: char* const record_file

        Filename of input data file with records, :option:`RecordFile`
        
    .. cpp:member:: fitsfile *fptr

        FITS object with pointer to data file
        
    .. cpp:member:: char* const library_file

        File name of calibration library, :option:`LibraryFile`
        
    .. cpp:member:: char* const event_file
    
        File name of output events (with reconstructed energy), :option:`TesEventFile`
        
    .. cpp:member:: int pulse_length
    
        Pulse length, :option:`PulseLength`
        
    .. cpp:member:: double scaleFactor
    
        Detection scale factor for initial filtering, :option:`scaleFactor`
        
    .. cpp:member:: double samplesUp
    
        Number of samples for threshold trespassing, :option:`samplesUp`
        
    .. cpp:member:: double nSgms
    
        Number of standard deviations in the kappa-clipping process for threshold estimation, :option:`nSgms`
    
    .. cpp:member:: int mode
    
        Calibration run (0) or energy reconstruction run (1), :option:`mode`
    
    .. cpp:member:: double LrsT
    
        Running sum length for the RS raw energy estimation (seconds), :option:`LrsT`
    
    .. cpp:member:: double LbT
    
        Baseline averaging length for the RS raw energy estimation (seconds), :option:`LbT`
    
    .. cpp:member:: char* const noise_file
    
        Noise file, :option:`NoiseFile`
    
    .. cpp:member:: char* filter_domain
    
        Filtering Domain: Time (**T**) or Frequency (**F**), :option:`FilterDomain`
    
    .. cpp:member:: char* filter_method
    
        Filtering Method: **F0** (deleting the zero frequency bin) or **F0** (deleting the baseline), :option:`FilterMethod`
    
    .. cpp:member:: char* energy_method
    
         Energy calculation Method: **OPTFILT**, **WEIGHT**, **WEIGHTN**, **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED**, :option:`EnergyMethod`
        
    .. cpp:member:: int lagsornot
    
        Lags (1) or no lags (0), :option:`LagsOrNot`
    
    .. cpp:member:: int ofiter
    
        Iterate (1) or not iterate (0), :option:`OFIter`
    
    .. cpp:member:: char oflib
    
        Work or not with a library with optimal filters (1/0), :option:`OFLib`
    
    .. cpp:member:: char* ofinterp
    
        Optimal Filter by using the Matched Filter or the DAB as matched filter (*MF*/*DAB*)
        It has been fixed in :ref:`tesreconstruction` as *DAB*
    
    .. cpp:member:: char* oflength_strategy
    
        Optimal Filter length Strategy: **FREE**, **BASE2**, **BYGRADE** or **FIXED**, :option:`OFStrategy`
    
    .. cpp:member:: int oflength
    
        Optimal Filter length (taken into account if :option:`OFStrategy` = **FIXED**), :option:`OFLength`
    
    .. cpp:member:: double monoenergy
    
        Monochromatic energy of input file in eV (only for library creation), :option:`monoenergy`
    
    .. cpp:member:: int interm
    
        Write or not intermediate files (1/0), :option:`intermediate`
    
    .. cpp:member:: char* const detectFile
    
        Intermediate detections file (if :option:`intermediate` = 1), :option:`detectFile`
    
    .. cpp:member:: char* const filterFile
    
        Intermediate filters file (if :option:`intermediate` = 1), :option:`filterFile`
        
    .. cpp:member:: char clobber
    
        Overwrite or not output files if exist (1/0), :option:`clobber`
    
    .. cpp:member:: int maxPulsesPerRecord
    
        Default size of the event list, :option:`EventListSize`
        
    .. cpp:member:: double SaturationValue
    
        Saturation level of the ADC curves, :option:`SaturationValue`
    
    .. cpp:member:: int tstartPulse1
    
        Tstart (samples) of the first pulse (different from 0 if the tstartPulsei input parameters are going to be used), :option:`tstartPulse1`
    
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
    - It is possible to find the tstart of the first pulse... 
    
      Obtain tstart of the first pulse in the derivative if :math:`derivative_i>threshold`
      
    - ... Or to use the tstart provided as input parameter

    **Members/Variables**
    
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
    
    .. cpp:member:: inparam* taskPars
    
        Instance of *inparam* structure storing input parameters
        
    .. cpp:member:: int np
    
        Number of parameters
    
    .. cpp:member:: string task
    
        Tool name
        
        
.. cpp:function:: int interpolatePOS(gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function interpolates an input vector (:cpp:member:`x_in`, :cpp:member:`y_in`), creating an output vector (:cpp:member:`x_out`, :cpp:member:`y_out`) with the size 
    and frequency step given. *POS* comes from the fact that the input spectrum only has positive frequencies (in order to not handle the f=0 bin).

    - Declare and initialize variables
    - GSL method applied for interpolatation
    - Generate the interpolated output vector
    - Free memory

    **Members/Variables**
    
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
        
.. _J:

.. _K:

.. _L:
    
.. cpp:function:: int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double)
        
    Located in file: *tasksSIRENA.cpp*
    
    This fucntion loads the structure :cpp:member:`record` into the :cpp:member:`adc_double` GSL vector.

    It checks if the record has been filled out with 0's :math:`\Rightarrow` It only loads the first values (which are different from 0).

    **Members/Variables**

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

    .. cpp:member:: gsl_vector** invector

        Input/Output GSL vector (non-filtered input vector/filtered input vector) 

    .. cpp:member:: int szVct

        Size of :cpp:member:`invector`

    .. cpp:member:: int sampleRate

        Sampling rate (samples/s)    
            
.. _M:

.. cpp:function::  int matrix2vector(gsl_matrix *matrixin, gsl_vector **vectorout)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function converts an input :math:`[n \times n]` square matrix into an output :math:`n^2` vector. It puts the first row of the matrix (:math:`n` elements) in the first :math:`n` elements of the vector (from :math:`0` to :math:`n-1`), the second row of the matrix in the elements from :math:`n` to :math:`2n-1` of the vector and so on.

    **Members/Variables**
    
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

.. _N:

.. cpp:function:: extern_C_OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status)    
    
    Located in file: *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *OptimalFilterSIRENA* data structure.
    
    **Members/Variables**
    
    .. cpp:member:: int* const status

        Input/output status
        
        
.. cpp:function:: extern_C_PulsesCollection* newPulsesCollection(int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *PulsesCollection* data structure. 
    
    **Members/Variables**
    
    .. cpp:member:: int* const status
        
        Input/output status
        
        
.. cpp:function:: extern_C_ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status)
    
    Located in file *integraSIRENA.cpp*
    
    Constructor. It returns a pointer to an empty *ReconstructInitSIRENA* data structure.
    
    **Members/Variables**
    
    .. cpp:member:: int* const status
        
        Input/output status
    
.. _O:

.. _P:


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
    
    .. cpp:member:: gsl_vector *x_fit
    
        Input GSL with *x* vector
        
    .. cpp:member:: gsl_vector *y_fit
    
        Input GSL with *y* vector
        
    .. cpp:member:: double *a
    
        Fit coefficient of the quadratic term
        
    .. cpp:member:: double *b
    
        Fit coefficient of the linear term    
        
    .. cpp:member:: double *c
    
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
    
    .. cpp:member:: gsl_vector *x_fit
    
        Input GSL with *x* vector
        
    .. cpp:member:: gsl_vector *y_fit
    
        Input GSL with *y* vector
        
    .. cpp:member:: double *a
    
        Fit coefficient of the linear term
        
    .. cpp:member:: double *b
    
        Fit coefficient (independent term)    

.. cpp:function:: void print_error( const char* const func, string message, int status)

    Located in file: *genutils.cpp*

    This function prints out error messages.

    **Members/Variables**
    
    .. cpp:member:: const char* const func
    
        Function name whose error is printed 
        
    .. cpp:member:: string msg
    
        Error message to be printed 

    .. cpp:member::  int status
        
        Status

        
.. cpp:function:: int procRecord (ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, PulsesCollection *foundPulses)
        
    Located in file: *tasksSIRENA.cpp*
    
    This function processes the input record:

    1) Declare and initialize variables
    
    2) Allocate GSL vectors
    
    3) (Low-pass filtering and) differentiation
    
    4) Find the events (pulses) in the record
    
    5) Calculate the end time of the found pulses and check if the pulses are saturated
    
    6) Obtain the approximate rise and fall times of each pulse (to be done)
    
    7) Load the found pulses data in the input/output *foundPulses* structure
    
    8) Write test info (if *reconstruct_init->intermediate* = 1)
    
    9) Write pulses info in intermediate output FITS file (if *reconstruct_init->intermediate* = 1)
    
    10) Free allocated GSL vectors

    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  In particular, this function uses parameters 
        to filter (:option:`scaleFactor`), to find pulses (:option:`PulseLength`, :option:`samplesUp`, :option:`nSgms`, :option:`LrsT`, :option:`LbT`, 
        :option:`EventListSize`) and to write info (:option:`detectFile`) if :option:`intermediate` = 1

    .. cpp:member:: double tstartRecord 

        Starting time of the record (in order to calculate absolute times)

    .. cpp:member:: double samprate 

        Sampling rate (in order to low-pass filter)

    .. cpp:member:: fitsfile* dtcObject

        Object which contains information of the intermediate FITS file (to be written if :option:`intermediate` = 1)

    .. cpp:member:: gsl_vector* record

        GSL vector with signal values of input record

    .. cpp:member:: PulsesCollection* foundPulses 

        Input/output structure where the info about found pulses is stored
        
        
.. cpp:function:: int pulseGrading(ReconstructInitSIRENA *reconstruct_init, int grade1, int grade2_1, int OFlength_strategy, int *pulseGrade, long *OFlength)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function provides the pulse grade (Pileup=-2, Rejected=-1, HighRes=1, MidRes=2, LowRes=3) and the optimal filter length by taking into account the info read from the XML file and the :option:`OFStrategy` (**FREE**, **BASE2**, **BYGRADE** or **FIXED**).

    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses the  :option:`OFLength` and *grading*
     
    .. cpp:member:: int grade1
     
        Pulse duration (length of optimal filter applied)
     
    .. cpp:member:: int grade2_1
    
        Difference between the start time of the pulse and the start time of the previous pulse
        
    .. cpp:member:: int OFlength_strategy
     
        Same as :option:`OFStrategy` (input)
        
    .. cpp:member:: int* pulseGrade
    
        Pulse grade (output)
    
    .. cpp:member:: long* OFlength
    
        Optimal filter length (= :option:`OFLength` only if :option:`OFStrategy` = **FIXED**) (output)
        
        
.. _Q:

.. _R:

.. cpp:function:: int readAddOrderParams(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int eventcntLib, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight)
    
    Located in file: *tasksSIRENA.cpp*

    This function reads the library data, add new data and sort the data according to an energy-ascending order.

    - Declare variables
    - Load values already in the library
    - Add new values
    - Realign
    - Add intermeadiate values
    - Recalculate intermediate values of some new pairs
    - Write values in the library
    - Free allocated GSL vectors
    
    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses :option:`FilterMethod`, :option:`PulseLength`, *library_collection*, :option:`monoenergy`, :option:`mode` and *noise_spectrum*
     
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
        
        
.. cpp:function:: int readFitsComplex(IOData obj, gsl_matrix **result)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a complex column of a FITS file. After that, the function puts them into a GSL matrix for an easier processing.
    
    **Members/Variables**
    
    .. cpp:member:: IOData obj
    
        Input object for complex FITS column
    
    .. cpp:member:: gsl_matrix** result
    
        Output GSL matrix
        
        
.. cpp:function:: int readFitsSimple(IOData obj, gsl_vector **result)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a simple column of a FITS file. After that, the function puts them into a GSL vector for an easier processing.
    
    **Members/Variables**
    
    .. cpp:member:: IOData obj
    
        Input object for simple FITS column
    
    .. cpp:member:: gsl_vector** result
    
        Output GSL vector
        
    
.. cpp:function:: extern_C_void reconstructRecordSIRENA(TesRecord* record, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init,  int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status)
    
    Located in file: *integraSIRENA.cpp*
    
    This function is the main wrapper function to detect, grade and calculate energy of pulses in input records.
    
    - Inititalize structure *PulsesCollection*
    - Check consistency of some input parameters
    - Detect pulses in input record (:cpp:func:`runDetect`). Save information of detected pulses
    
        - If PCA, pulses energies are already written in the :cpp:member:`pulsesAll` structure

    - If in RECONSTRUCTION (:option:`mode` = 1):
    
        - Filter record and calculate energy of pulses (:cpp:func:`runEnergy`)
        
    - Populate output event list with pulses energies, arrival time and grading
    
    **Members/Variables**
    
    .. cpp:member:: TesRecord* record 

        Instance of *TesRecord* structure that contains the input record 
    
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
        
        
.. cpp:function:: void runDetect(TesRecord* record, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
    
    Located in file: *tasksSIRENA.cpp*

    This function is responsible for the **detection** in SIRENA, record by record. It is used both for library creation (:option:`mode` = 0) and energy reconstruction (:option:`mode` = 1) runnings.

    Conditions:   

        - If first record and :option:`mode` = 1  :math:`\Rightarrow`  Run :cpp:func:`filderLibrary`

        - If last record and :option:`mode` = 0 :math:`\Rightarrow` Run :cpp:func:`calculateTemplate` and :cpp:func:`writeLibrary`

        - If :option:`intermediate` = 1 :math:`\Rightarrow` :cpp:func:`writeTestInfo` and :cpp:func:`writePulses`

        - If :option:`mode` = 0 :math:`\Rightarrow` Find pulses by using :cpp:func:`findPulsesCAL`
        
        - If :option:`mode` = 1 :math:`\Rightarrow` Find pulses by :cpp:func:`InitialTriggering` and :cpp:func:`FindSecondaries`
        
    Steps:

        1) Create library file if it is necessary: last record and :option:`mode` = 0 (run :cpp:func:`createLibrary`)           

        2) Create intermediate output FITS file if required (:cpp:func:`createDetectFile`)     

        3) (Filter and) differentiate the *models* of the library (only for the first record in :option:`mode` = 1 run  (:cpp:func:`filderLibrary`)

        4) Store the input record in *invector* (:cpp:func:`loadRecord`)

        5) Convert *I* into *R* if :option:`EnergyMethod` = **I2R**, **I2RALL**, **I2RNOL** or **I2RFITTED** (:cpp:func:`convertI2R`)
        
        6) Process each record (:cpp:func:`proceRecord`): 

                - (Low-pass filter and) differentiate                          
                - Find pulses                                                           
                - Load the found pulses data in the input/output *foundPulses* structure
                - Write test info in intermediate output FITS file if *reconstruct_init->intermediate* = 1 (:cpp:func:`writeTestInfo`)                            
                - Write pulses info in intermediate output FITS file if *reconstruct_init->intermediate* = 1 (:cpp:func:`writePulses`)
                
           **From this point forward, I2R, I2RALL, I2RNOL and I2RFITTED are completely equivalent to OPTFILT**     
           
        7) If last record in :option:`mode` = 0 run:                                         

                * :cpp:func:`calculateTemplate` (and :cpp:func:`weightMatrix`)
                * :cpp:func:`writeLibrary`
                
        8) If last record and PCA
        
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

    .. cpp:member:: TesRecord* record 

        Member of *TesRecord* structure that contains the input record 

    .. cpp:member:: int lastRecord

        Integer to verify whether record is the last one (=1) to be read (and thus if library file will be created)
        
    .. cpp:member:: PulsesCollection* pulsesAll

        Member of *PulsesCollection* structure to successively store all the pulses used to create the library. Re-populated after each processed record
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)

    .. cpp:member:: PulsesCollection** pulsesInRecord

        Member of *PulsesCollection* structure to store all the pulses found in the input record

                
.. cpp:function:: void runEnergy(TesRecord* record, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function calculates the pulse energy applying different methods (from :option:`EnergyMethod`).
    It only runs in RECONSTRUCTION mode (:option:`mode` = *1*) (except to :option:`EnergyMethod` = **PCA**).

    - Declare variables
    - Store the :cpp:member:`record` in *invector* (:cpp:func:`loadRecord`)
    - Subtract the baseline if :option:`EnergyMethod` = **OPTFILT** and *runF0orB0val* = 1 (:option:`FilterMethod` = **B0**)
    - Subtract the baseline if :option:`EnergyMethod` = **WEIGHT**
    - Check Quality
    - For each pulse:

        - Establish the pulse grade (HighRes=0, MidRes=1, LowRes=-1) and the optimal filter length
        - Pulse: Load the proper piece of the record in *pulse*
        - If :option:`OFIter` = 1, in the first iteration ( *numiteration* = 0) the values of *maxDER* and *maxDERs* are used in 
          :cpp:func:`find_matchedfilter`, :cpp:func:`find_matchedfilterDAB`, :cpp:func:`find_optimalfilter`, :cpp:func:`find_optimalfilterDAB` or :cpp:func:`find_Esboundary` getting the values of the *energies* which straddle the *maxDER* (*Ealpha* and *Ebeta*). It will have more iterations if the calculated *energy* is out of *[Ealpha, Ebeta]*. If *energy* is in *[Ealpha, Ebeta]* the iterative process stops.
            
                - If :option:`EnergyMethod` = **OPTFILT** (or **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED**) and :option:`OFLib` = 0:
                
                    - Filter (find the matched filter and load it in *filter*) 
                    
		      - If *MF* :math:`\Rightarrow` :cpp:func:`find_matchedfilter`  
		      - If *DAB* :math:`\Rightarrow` :cpp:func:`find_matchedfilterDAB`
		      
                    - Calculate the optimal filter
                    
                - If :option:`EnergyMethod` = **OPTFILT** (or **I2R**, **I2RALL**, **I2RNOL**, **I2RFITTED**) and :option:`OFLib` = 1:
                
                    - Choose the base-2 system value closest (lower than or equal) to the pulse length
                    - Filter (find the optimal filter and load it in *filter*)
                      
                      - If *MF* :math:`\Rightarrow` :cpp:func:`find_optimalfilter`  
		      - If *DAB* :math:`\Rightarrow` :cpp:func:`find_optimalfilterDAB`
		      
                    - If :option:`FilterDomain` = **F** :math:`\Rightarrow` Bin0 = 0 and conjugate
                    
                - If :option:`EnergyMethod` = **WEIGHT** or **WEIGHTN**:
                
                    - Get the indexes of the two energies which straddle the pulse
                    - If :option:`EnergyMethod` = **WEIGHTN** and :option:`OFLib` = 1:
                    
                       - Choose the base-2 system value closest (lower than or equal) to the pulse length
                       - :cpp:func:`find_prclwn` to find the appropriate values of the PRECALWN HDU (**PRCLx** columns)
                    
                - Calculate the energy of each pulse
        - Subtract the pulse model from the record
        - Write info of the pulse in the output intemediate file if :option:`intermediate` = 1
    - Free allocated GSL vectors

    **Members/Variables**

    .. cpp:member:: TesRecord** record

        Structure that contains the input ADC record
        
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values)
    
    .. cpp:member:: PulsesCollection* pulsesInRecord

        Collection of pulses found in the current record
        
    .. cpp:member:: OptimalFilterSIRENA **optimalFilter
    
        Optimal filters used in reconstruction

.. _S:

.. cpp:function:: int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
    
    Located in file: *tasksSIRENA.cpp*

    This function returns as :cpp:member:`vectorout` the :cpp:member:`vectorin` delayed :cpp:member:`m` samples.

    **Members/Variables**

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
    
    .. cpp:member:: int m

        Advance in samples

    .. cpp:member:: gsl_vector* vectorin

        GSL vector with input vector

    .. cpp:member:: gsl_vector* vectorout 

        GSL with input vector (:cpp:member:`vectorin`) moved forward :cpp:member:`m` samples
            
.. _T:

.. cpp:function:: int toGslMatrix(void **buffer, gsl_matrix **matrix, long numCol, int numRow, int type, int eventini)
    
    Located in file: *inoututils.cpp*
    
    The function puts the values of the input buffer into an output GSL matrix. Columns and rows are input parameters.
    
    **Members/Variables**
    
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
    
    This function converts an input :math:`n^2` vector into an output :math:`[n \times n]` square matrix. It puts the first :math:`n` elements of the vector in the first row of the matrix, the second group of :math:`n` elements (from :math:`n` to :math:`2n-1`) of the vector in the second row and so on.
    
    **Members/Variables**
    
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
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  In particular, this function uses :option:`PulseLength` and :option:`EnergyMethod`
        
    .. cpp:member:: bool saturatedPulses

        If *true*, all the pulses ( :option:`mode` = 0 :math:`\Rightarrow` all the pulses have the same energy) are saturated
        
    .. cpp:member:: PulsesCollection* pulsesAll 

        Collection of pulses found in the previous records
        
    .. cpp:member:: PulsesCollection* pulsesInRecord

        Collection of pulses found in the current record
        
    .. cpp:member:: long nonpileupPulses

        Number of non piled-up pulses
        
    .. cpp:member: gsl_vector* nonpileup

        GSL vector containing info about all the pulses informing if they are piled-up or not
        
    .. cpp:member:: gsl_vector** pulseaverage

        GSL vector with the pulseaverage (= template = model) of the non piled-up pulses

    .. cpp:member:: gsl_matrix** covariance

        GSL matrix with covariance matrix
        
    .. cpp:member:: gsl_matrix** weight

        GSL matrix with weight matrix
    
    
.. cpp:function:: int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double normalizationFactor, double energy, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function runs in RECONSTRUCTION mode and writes the optimal filter info (in the FILTER HDU) for each pulse
    if :option:`intermediate` = 1.

    - Declare variables
    - Open intermediate FITS file
    - Create the FILTER HDU if it is the first pulse
    - Write data
    
      - **OPTIMALF** column
      
      - **OFLENGTH** column
      
      - **NFCTR** column
      
      - **FREQ** column
      
      - **OPTIMALFF** column
      
      - **ENERGY** column
      
    - Modify the Primary HDU
    - Close intermediate output FITS file if it is necessary
    - Free memory

    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses :option:`detectFile`, 
        :option:`clobber` and :option:`PulseLength`
        
    .. cpp:member:: int pulse_index
    
        Index of the pulse whose info is going to be written (to know if it is the first pulse)
        
    .. cpp:member:: double normalizationFactor
    
        Normalization factor in the expression of the optimal filtering (see :ref:`optimalFilter`)
        
    .. cpp:member:: double energy
    
        Estimated energy (eV)
        
    .. cpp:member:: gsl_vector* optimalfilter 
    
        Optimal filter in time domain
        
    .. cpp:member:: gsl_vector* optimalfilter_f
    
        Frequency axis of the optimal filter spectrum
        
    .. cpp:member:: gsl_vector* optimalfilter_FFT
    
        Optimal filter spectrum
        
    .. cpp:member:: fitsfile** dtcObject 
    
        Fitsfile object for intermeadiate file name
    
    
.. cpp:function:: int writeFitsComplex(IOData obj, gsl_matrix *matrix)
    
    Located in file: *inoututils.cpp*
    
    This function reads values of a GSL matrix. After that, the function puts them into a complex column of the output FITS file.
    
    **Members/Variables**
    
    .. cpp:member:: IOData obj
    
        Object for FITS column to be written
        
    .. cpp:member:: gsl_matrix* matrix
    
        Input GSL matrix with data
    
    
.. cpp:function:: int writeFitsSimple(IOData obj, gsl_vector *vector)
    
    Located  in file: *inoututils.cpp*
    
    This function reads values of a GSL vector.	After that, the function puts them into a column of the output FITS file.
    
    **Members/Variables** 
    
    .. cpp:member:: IOData obj
        
        Object for FITS column to be written
        
    .. cpp:member:: gsl_vector* vector
        
        Input GSL vector with data
        
    
.. cpp:function:: int writeLibrary(ReconstructInitSIRENA *reconstruct_init, double samprate, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject)
    
    Located in file: *tasksSIRENA.cpp*
    
    This function writes the library (reordering if it is necesary and calculating some intermediate parameters)

        - Adding a new row to the library
        - Write the first row of the library

    **Members/Variables**
    
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
    
    
.. cpp:function:: void writeLog(FILE *fileRef, string type, int verbosity, string message)

    Located  in file: *inoututils.cpp*
    
    This function includes the processing of the each level of message in the log file and the output screen:
 
      - Verbosity = 0 :math:`\Rightarrow` The log file and the output screen include Errors
      - Verbosity = 1 :math:`\Rightarrow` The log file and the output screen include Errors and Warnings
      - Verbosity = 2 :math:`\Rightarrow` The log file and the output screen include Errors, Warnings and Alerts
      - Verbosity = 3 :math:`\Rightarrow` The log file and the output screen include Errors, Warnings, Alerts and Log messages

    **Members/Variables**
    
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

    This function writes the data of the pulses found in the record in the intermediate FITS file (in the **PULSES** HDU). The pulses info given is: **TSTART**, **I0** (the pulse itself), **TEND**, **TAURISE**, **TAUFALL** and **QUALITY**.
    
    **Members/Variables**
    
    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values). In particular, this function uses :option:`detectFile`, 
        :option:`PulseLength` and :option:`clobber`.

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
        
        1 :math:`\Rightarrow` Truncated pulses
        
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

    This function writes the **TESTINFO** HDU in the intermediate FITS file. The written columns are **FILDER** (low-pass filtered and differentiated record) and **THRESHOLD**.

    **Members/Variables**

    .. cpp:member:: ReconstructInitSIRENA** reconstruct_init

        Member of *ReconstructInitSIRENA* structure to initialize the reconstruction parameters (pointer and values).  In particular, this function uses :option:`detectFile`.
        
    .. cpp:member:: gsl_vector* recordDERIVATIVE 

        GSL vector with input record (low-pass filtered and) differentiated
        
    .. cpp:member:: double threshold

        Threshold value used to find pulses
        
    .. cpp:member:: fitsfile dtcObject

        Object which contains information of the intermediate FITS file
   
.. _X:

.. _Y:

.. _Z:




    
    
    
.. |mem| replace:: **Members/Variables**    
    
    
    
    
    
    
