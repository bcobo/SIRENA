.. Description of SIRENA tools command line

.. role:: bred
.. role:: red
.. role:: blue

.. _SIRENAtools:

##########################
SIRENA Tools CLI
##########################

Here we describe the command line options for the tools run for SIRENA reconstruction (``gennoisespec`` and ``tesreconstruction``).
At the end, links to the documentation of other SIXTE tools required for the simulation of data files are provided (``tesconstpileup``) and XIFUSIM.

.. _gennoisespec: 

gennoisespec
=============

The goal of the ``gennoisespec`` tool is to calculate the current noise spectral density and the noise weight matrixes.
The input data from which it would be calculated should be a FITS file with the data splitted into records (see :ref:`noise-records <noise-records>`) with or without photon events (pulses).

The user must supply the following input parameters:

.. _gennoisePars:

.. option:: inFile=<str>

	Name of the input FITS file (stream splitted into records).

	Default: *a.fits*

.. option:: outFile=<str>

	Name of the output FITS file with the current noise spectral density.

	Default: *a_noisespec.fits*

.. option:: intervalMinSamples=<int>

	Minimum length (in samples) of a pulse-free interval to use. 
	
	It will be redefined as the base-2 system value closest-lower than or equal than :option:`intervalMinSamples`.
	
	Default: 8192

.. option:: nplPF=<real>

	Number of pulse lengths after the end of the pulse to start the pulse-free interval searching (only relevant if pulse detection in the stream has to be performed).

	Default: 0

.. option:: nintervals=<int>

	Number of pulse-free intervals to use for the noise average.

	Default: 1000

.. _scaleFactor_gennoisespec:

.. option:: scaleFactor=<real>
        
	Scale factor to apply to make possible a variable cut-off frequency of the low-pass filter. In fact, the cut-off frequency of the filter is :math:`1/(\pi \cdot sF)` and therefore, the box-car length is :math:`\pi \cdot sF \cdot samprate` (see :ref:`Low-Pass filtering <lpf>`).
	
	If the :option:`scaleFactor` makes the box-car length :math:`\leq 1` is equivalent to not filter (cut-off frequency of the low-pass filter is too high). If the :option:`scaleFactor` is too large, the low-pass filter band is too narrow, and not only noise is rejected during the filtering, but also the signal.
	
	Default: 0

.. _samplesUp_gennoisespec:

.. option:: samplesUp=<int>

	Consecutive samples that the signal must cross over the threshold to trigger a pulse detection (only relevant if pulse detection in the stream has to be performed).

	Default: 3

.. _nSgms_gennoisespec:

.. option:: nSgms=<real> 

	Number of quiescent-signal standard deviations to establish the threshold through the *kappa-clipping* algorithm (only relevant if pulse detection in the stream has to be performed).

	Default: 3.5

.. option:: pulse_length=<int> 

	Pulse length in samples. 

	Default: 8192

.. _LrsT_gennoisespec:

.. option:: LrsT=<secs> 

	Running sum (RS) length for the RS-filtering for raw energy estimation, in seconds. 

	Default: 3.E-5

.. _LbT_gennoisespec:

.. option:: LbT=<secs> 

	Baseline averaging length, in seconds. 

	Default: 1.E-3
	
.. option:: weightMS=<yes|no> 

	Calculate and write the weight matrixes if *yes*.

	Default: *no*
	
.. option:: namelog=<str>

	Output log file name. 

	Default: *noise_log.txt*

.. _clobber_gennoisespec:

.. option:: clobber=<yes|no> 
	
	Overwrite output files if they exist. 

	Default: *no*

.. option:: verbosity=<1|2|3> 

	Verbosity level of the output log file. 

	Default: 3

.. option:: matrixSize=<int> 

	Size of noise matrix if only one to be calculated, in samples. 

	Default: 0

.. option:: rmNoiseInterval=<yes|no> 

	Remove some noise intervals before calculating the noise spectrum if *yes*.

	Default: *no*

A typical command line run of this tool would be:

::

	> gennoisespec inFile=noise.fits outFile=noiseSpec.fits intervalMinSamples=pulseLength \
    		pulse_length=pulseLength nintervals=1000 

The sampling rate is calculated by using some keywords in the input FITS file. In case of ``tessim`` simulated data files, using the ``DELTAT`` keyword *samplingRate=1/deltat*. In case of ``xifusim`` simulated data files, every detector type defines a master clock-rate ``TCLOCK`` and the sampling rate is calculated either from a given decimation factor ``DEC_FAC`` (FDM and NOMUX) as *samplingRate=1/(tclock·dec_fac)*, or from the row period  ``P_ROW`` and the number of rows ``NUMROW`` (TDM) as *samplingRate=1/(tclock·numrow·p_row)*. In case of old simulated files, the sampling rate could be read from the ``HISTORY`` keyword in the *Primary* HDU. If the sampling frequency can not be get from the input file after all, a message will ask the user to include the ``DELTAT`` keyword (inverse of the sampling rate) in the input FITS file before running again.

.. _outNoise:

The output FITS file contains three HDUs, *NOISE*, *NOISEALL* and *WEIGHTMS*.
The *NOISE* HDU contains three columns:

* **FREQ**: Noise positive frequencies in Hz

* **CSD**: Current noise spectral density. Amount of current per unit of frequency (spectral density) in :math:`A/\sqrt(Hz)`

* **SIGMACSD**: CSD Standard error of the mean in :math:`A/\sqrt(Hz)` (not filled yet)

The *NOISE* HDU contains two keywords:

* ``BSLN0``: Noise baseline (it will be propagated to the library as ``BASELINE`` in the *Library* HDU when building the library FITS file)

* ``NOISESTD``: Noise standard deviation 

The *NOISEALL* HDU contains **FREQ** and **CSD** columns for positive and negative frequencies.

If :option:`weightMS` = *yes*, the *WEIGHTMS* HDU contains **Wx** columns. The lengths *x* will be base-2 values and will vary from the base-2 system value closest-lower than or equal-to the :option:`intervalMinSamples` decreasing until 2. If :option:`matrixSize` is different from 0, only the **Wx** column being *x* equals to :option:`matrixSize` is calculated (although the rest columns appear in the HDU, they are filled with 0's).


.. _tesreconstruction:


tesreconstruction
=================

The ``tesreconstruction`` tool is a wrapper to perform the energy reconstruction of the photon events by means of two different implementations: ``Rcmethod=PP`` runs the preliminary branch developed by Philippe Peille and ``Rcmethod=SIRENA`` runs the SIRENA code in this documentation.

SIRENA code takes a FITS input file of data, optionally performs the detection of the events, then grades them and finally reconstructs their energy following the algorithm selected by the user in the input command line of ``tesreconstruction``.

The :ref:`input data <inputFiles>` should be a FITS file with the data splitted into :ref:`records <records>`. 

To run SIRENA implementation, the user must supply the following input parameters (see :ref:`reconMethods` for a detailed description in the context of the reconstruction methods to which they apply):


.. _tesreconPars:


.. option:: Rcmethod=<SIRENA>

	SIRENA Reconstruction method.

.. option::  RecordFile=<str>

	Input record FITS file.
	
	Default: *record.fits*

.. option::  TesEventFile=<str>

	Output event list FITS file.
	
	Default: *event.fits*

.. option::  PulseLength=<int>

	Pulse length in samples.
	
	Default: 8192

.. option::  EventListSize=<str> 

	Default size of the event list. 
 
	Default: 1000

.. option::  LibraryFile=<str>

	FITS file with calibration library. 

	Default: *library.fits*

.. option::  scaleFactor=<real> 
	
	Scale factor to apply to make possible a variable cut-off frequency of the low-pass filter. In fact, the cut-off frequency of the filter is :math:`1/(\pi \cdot sF)` and therefore, the box-car length is :math:`\pi \cdot sF \cdot samprate` (see :ref:`Low-Pass filtering <lpf>`).
	
	If the :option:`scaleFactor` makes the box-car length :math:`\leq 1` is equivalent to not filter (cut-off frequency of the low-pass filter is too high). If the :option:`scaleFactor` is too large, the low-pass filter band is too narrow, and not only noise is rejected during the filtering, but also the signal.
	
	Default: 0

.. option::  samplesUp=<int> 

	Number of consecutive samples up for threshold trespassing (only used in calibration run, and in production run with STC detection mode).

	Default: 3
	
.. option::  samplesDown=<int> 

	Number of consecutive samples below the threshold to look for other pulse (only used in production run with STC detection mode).

	Default: 4

.. option::  nSgms=<real> 

	Number of quiescent-signal standard deviations to establish the threshold through the kappa-clipping algorithm.

	Default: 3.5

.. option::  detectSP=<0|1>

	Detect secondary pulses (1) or not (0).

	Default: 1
	
.. option::  LrsT=<secs>

	Running sum (RS) length for the RS raw energy estimation, in seconds (only used in calibration run).
	
	Default: 30E-6

.. option::  LbT=<secs>

	Baseline averaging length, in seconds.

	Default: 1.E-3

.. option::  monoenergy=<eV>

	Monochromatic energy of the pulses in the input FITS file in eV (only used in calibration run).
	
.. option::  hduPRECALWN=<yes|no>

	Add or not the *PRECALWN* HDU in the library file (only used in calibration run).

	Default: *no*	

.. option::  hduPRCLOFWM=<yes|no>

	Add or not the *PRCLOFWM* HDU in the library file (only used in calibration run).

	Default: *no*	
	
.. option::  largeFilter=<int>

	Length (in samples) of the longest fixed filter (only used in calibration run).
	
	Default: -999
	
.. option:: opmode=<0|1>

	Calibration run for library creation (0) or energy reconstruction run (1).

	Default: 1
	
.. option:: detectionMode=<AD | STC>

	Adjusted Derivative (AD) or Single Threshold Crossing (STC). Not used in library creation mode (:option:`opmode` = 0).

	Default: STC

.. option::  NoiseFile=<str>

	Noise FITS file with noise spectrum. 

	Default: *noise.fits*

.. option::  FilterDomain=<T | F> 

	Filtering Domain: Time(T) or Frequency(F). Not used in library creation mode (:option:`opmode` = 0).

	Default: *T*

.. option::  FilterMethod=<F0 | B0>
	
	Filtering Method: *F0* (deleting the zero frequency bin) or *B0* (deleting the baseline). 

	Default: *F0*

.. option::  EnergyMethod=<OPTFILT | WEIGHT | WEIGHTN | I2R | I2RALL | I2RNOL | IRFITTED | PCA>

	:ref:`reconMethods` Energy calculation Method: OPTFILT (Optimal filtering), WEIGHT (Covariance matrices), WEIGHTN (Covariance matrices, first order), I2R, I2RALL, I2RNOL and I2RFITTED (Linear Transformations), or PCA (Principal Component Analysis). Not used in library creation mode (:option:`opmode` = 0).
	
	If :option:`EnergyMethod` = OPTFILT and :option:`PulseLength` < :option:`OFLength`, 0-padding is applied (:option:`OFLength` length filters will be used but padding with 0's from :option:`PulseLength`).

	Default: *OPTFILT*
	
.. option::  filtEeV=<eV>

	Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R, I2RALL, I2RNOL and I2RFITTED).

	Default: 6000
	
.. option::  OFNoise=<NSD | WEIGHTM>

	It has only sense if :option:`EnergyMethod` = OPTFILT and it means to use the noise spectrum density (NSD) or the noise weight matrix (WEIGHTM).

	Default: *NSD*

.. option::  LagsOrNot=<0|1> 

	Use LAGS == 1 or NOLAGS == 0 to indicate whether subsampling pulse arrival time is required. Currently only implemented for :option:`EnergyMethod` =OPTFILT, and :option:`EnergyMethod` =WEIGHTN combined with :option:`OFLib` =yes.

	Default: 1

.. option::  nLags=<int> 

	Number of lags (samples) to be used if :option:`LagsOrNot` =1. It has to be a positive odd number.

	Default: 9

.. option::  Fitting35=<3|5> 

	Number of lags to analytically calculate a parabola (3) or to fit a parabola (5).

	Default: 3

.. option::  OFIter=<0|1>

	Iterate (1) or not iterate (0) to look for the closest energy interval. When iterations are activated, there will be more iterations if the calculated energy is out of the interval [Ealpha, Ebeta] straddling the predicted energy according the pulse shape.   

	Default: 0

.. option::  OFLib=<yes|no>  

	Work with a library with optimal filters (OFLib=yes) or instead do Optimal Filter calculation on-the-fly (OFLib=no).
	
	Default: yes 

.. option::  OFStrategy=<FREE | BYGRADE | FIXED> 

	Optimal Filter length Strategy: FREE (no length restriction), BYGRADE (length according to event grading) or FIXED (fixed length). These last 2 options are only for checking and development purposes; a normal run with *on-the-fly* calculations will be done with :option:`OFStrategy` = *FREE*.
	Only used if :option:`OFLib` =no. Not used in library creation mode (:option:`opmode` = 0). 

	Default: *BYGRADE*

.. option::  OFLength=<int> 

	Fixed Optimal Filter length (only if :option:`OFStrategy` = **FIXED**,  :option:`opmode` = 1 and :option:`OFLib` =no).

	Default: 8192
	
.. option::  preBuffer=<int> 

	Some samples added before the starting time of a pulse.

	Default: 0

.. option::  intermediate=<0|1>  

	Write intermediate files: Y(1), N(0)? 

	Default: 0

.. option::  detectFile=<str>

	Intermediate detections FITS file (if :option:`intermediate` = 1).

	Default: *detections.fits*
	
.. option::  errorT=<int> 

	Additional error (in samples) added to the detected time. Logically, it changes the reconstructed energies. For deveplopment purposes.

	Default: 0
	
.. option::  Sum0Filt=<0|1>  

	If 0-padding, subtract (1) or not subtract (0) the sum of the filter. For deveplopment purposes. 

	Default: 0

.. option::  tstartPulse1=<str> 
	
	Start time (in samples) of the first pulse (0  if detection should be performed by the system; greater than 0 if provided by the user) or file name containing the tstart (in seconds) of every pulse. For development purposes.

	Default: 0

.. option::  tstartPulse2=<int>  

	Start time (in samples) of the second pulse in the record (0  if detection should be performed by the system; greater than 0 if provided by the user). For development purposes.

	Default: 0

.. option::  tstartPulse3=<int> 
	
	Start time (in samples) of the third pulse in the record (0  if detection should be performed by the system; greater than 0 if provided by the user). For development purposes.

	Default: 0
	
.. option::  energyPCA1=<real>

	First energy (in eV) (only for PCA).
	
	Default: 500

.. option::  energyPCA2=<real>

	Second energy (in eV) (only for PCA).
	
	Default: 1000
	
.. option::  XMLFile=<str>

	XML input FITS file with instrument definition.

	Default: *xifu_pipeline.xml*
	
.. option::  clobber=<yes|no> 
	
	Overwrite output files if they exist.

	Default: *no*

.. option::  history=<yes|no> 

	Write program parameters into output FITS file.

	Default: *yes*

The sampling rate is calculated by using some keywords in the input FITS file. In case of ``tessim`` simulated data files, using the ``DELTAT`` keyword *samplingRate=1/deltat*. In case of ``xifusim`` simulated data files, every detector type defines a master clock-rate ``TCLOCK`` and the sampling rate is calculated either from a given decimation factor ``DEC_FAC`` (FDM and NOMUX) as *samplingRate=1/(tclock·dec_fac)*, or from the row period  ``P_ROW`` and the number of rows ``NUMROW`` (TDM) as *samplingRate=1/(tclock·numrow·p_row)*. In case of old simulated files, the sampling rate could be read from the ``HISTORY`` keyword in the *Primary* HDU or even from the input XML file. If the sampling frequency can not be get from the input files after all, a message will ask the user to include the ``DELTAT`` keyword (inverse of the sampling rate) in the input FITS file before running again.

The output file will also be a FITS file storing one event per row with the following information in the HDU named *EVENTS*:

* **TIME**: arrival time of the event (in s)

* **SIGNAL**: energy of the event (in keV)

* **AVG4SD**: average of the first 4 samples of the derivative of the pulse

* **ELOWRES**: energy provided by a low resolution energy estimator filtering with a 4-samples-length filter (in keV)

* **GRADE1**: length of the filter used, i.e., the distance to the following pulse (in samples) or the :option:`PulseLength` if the next event is further than this value or if there are no more events in the same record

* **GRADE2**: distance to the end of the preceding pulse (in samples). If pulse is the first event in the record, this is fixed to the :option:`PulseLength` value

* **PHI**: arrival phase (offset relative to the central point of the parabola) (in samples) 

* **LAGS**: number of samples shifted to find the maximum of the parabola

* **BSLN**: mean value of the baseline in general 'before' a pulse (according the value in samples of :option:`LbT`)

* **RMSBSLN**: standard deviation of the baseline in general 'before' a pulse (according the value in samples of :option:`LbT`)

* **PIX_ID**: pixel number

* **PH_ID**: photon number identification for cross matching with the impact list

* **RISETIME**: rise time of the event (in s)

* **FALLTIME**: fall time of the event (in s)

* **GRADING**: Pulse grade (HighRes=1, MidRes=2, LimRes=3, LowRes=4, Rejected=-1, Pileup=-2)

.. _xifusim:

xifusim
=======

http://www.sternwarte.uni-erlangen.de/research/sixte/ 

.. _tesconstpileup:

tesconstpileup
==============

http://www.sternwarte.uni-erlangen.de/research/sixte/


