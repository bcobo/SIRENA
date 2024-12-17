## Tests for SIRENA   

These tests are designed to be run with `pytest` after new code upload.

### Initial steps   
Before running the tests:
  - initialize SIMPUT, SIXTE and SIRENA
  - download [*input*](https://nextcloud.ifca.es/index.php/s/742eePYbDPatrso) and *reference*  [1](https://nextcloud.ifca.es/index.php/s/2Ssi9rDTZtKHooy) and [2](https://nextcloud.ifca.es/index.php/s/L6ZwppncyPKHCg5) data and untar them under *test* folder.

### Tests description
* `test_noise_spectrum.py`: test the creation of the noise spectrum (used by optimal filter)
* `test_noise_spectrum_matrices.py`: test the creation of the noise spectrum with the noise weight matrices (used by optimal filter) **-> run time > 1h**
* `test_library_optfilt.py`: creation of a library with an optimal filter at 6 keV
* `test_library_covar2E.py`: creation of a library prepared for covariance methods of reconstruction (2 monochromatic energies) **-> run time > 1h**
* `test_recons_mono_0pad_FULL.py`: test to check that reconstruction of monochromatic pulses with an optimal filter (OPTFILT) is equivalent to reconstruction with 0-padding method of the same length.
* `test_recons_pairs_optfilt_ns.py`: reconstruction of pairs of pulses with optimal filter and noise spectrum
* `test_recons_pairs_optfilt_wn.py`: reconstruction of pairs of pulses with optimal filter and noise weight matrix
* `test_recons_mono_intcovar.py`: reconstruction of monochromatic pulses with INTCOVAR method **-> run time > 1h**
* `test_recons_mono_covar.py`: reconstruction of monochromatic pulses with COVAR method

### Runing the tests
  Once the *input* and *reference* data files are located in their folders, go to `test` folder and run:
  ```
  test> pytest name_of_the_test.py
  ```
  
