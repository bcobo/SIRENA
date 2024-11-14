[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_bronze.png)](https://api.eu.badgr.io/public/assertions/5GEdTdkzR2KDlXkZTCgrLg "SQAaaS bronze badge achieved")
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md) 
# SIRENA
Package developed at the Instituto de Física de Cantabria ([IFCA](https://ifca.unican.es/en-us)), CSIC-UC to process events in Transition Edge Sensors (TES) detectors
like the [X-IFU](https://x-ifu.irap.omp.eu/) detector, to be onboard the future ESA mission [Athena](https://www.the-athena-x-ray-observatory.eu/en).

It reads as input simulated or laboratory data files with the pulses produced by the X-ray events on the TES detector.
Different algorithms can then be applied to process such events generating as output an event list with the energy, localtion, arrival time and grading of all the detected events.
A full description of the algorithms that SIRENA contains as well as as the guide to use it can be found in [SIRENA ReadTheDocs](https://sirena.readthedocs.io/en/latest/SIRENA.html)

The SIRENA installation procedure is as follows:

1. Required software and instrument files that must be installed and initialized following the instructions at their own sites:
   * [SIXTE v3.0](https://www.sternwarte.uni-erlangen.de/sixte/sixte-beta/) 
   * [XIFU instrument files](https://www.sternwarte.uni-erlangen.de/sixte/sixte-beta)
   * [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/) 
2. Clone SIRENA repository 
   ```
     > git clone https://github.com/bcobo/SIRENA.git
     > cd SIRENA
   ```
3. Configure project and create build directory:
   ```
     > cmake -S . -B build -DCMAKE_INSTALL_PREFIX=sirenadir -DSIMPUT_ROOT=simputdir -DSIXTE_ROOT=sixtedir
   ```
   where `sirenadir`, `simputdir` and `sixtedir` are the full paths to the SIRENA, SIMPUT and SIXTE installation folders respectively.
   
4. Compile and link
   ```
     > cmake --build build
   ```
5. Install SIRENA  
   ```
     > cmake --install build
   ``` 
6. Initialize SIRENA to use it (same as for SIXTE/SIMPUT):
   For bash shell:
   ```
     > export SIRENA=sirenadir
     > . ${SIRENA}/bin/sirena-install.sh
   ```
   For csh shell:
   ```
     > setenv SIRENA sirenadir
     > source ${SIRENA}/bin/sirena-install.csh
   ```

Help Desk: 
      ```
         cobo@ifca.unican.es
      ```
