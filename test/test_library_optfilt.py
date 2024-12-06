#############################################################################
# Runs the generation of a library with optimal filters  (NS & WM)  
# Then roughly compares this libraries with reference library files.
#############################################################################

import os
from subprocess import run
import glob

# Define environmental variables
sixte = os.environ["SIXTE"]
sirena = os.environ["SIRENA"]
headas = os.environ["HEADAS"]

# if sirena or sixte are not initiated send message and stop
assert sirena != "", "Please set the environmental variable SIRENA and initiate the software using (for sh-like shell): . $SIRENA/bin/sirena-install.sh"
assert sixte != "", "Please set the environmental variable SIXTE and initiate the software using (for sh-like shell): . $SIXTE/bin/sixte-install.sh"
assert headas != "", "Please set the environmental variable HEADAS and initiate the software using (for sh-like shell): . $SIXTE/headas-init.sh"

outdir = "./data/output"
indir = "./data/input"
reference = "./data/reference"
# read config file date form file config_version.txt
with open(f"{indir}/config_version.txt", "r") as file:
    date_config_file = file.read().strip()

# Find name of (unique) xml file in indir directory
xmlfile = glob.glob(f"{indir}/config*.xml")[0]

# Define input common variables
file_input = f"{indir}/mono6keV_1000p_{date_config_file}.fits"
template_eV = 6000

# Create a test for teslib tool to be used only for optimal filter(s) and noise spectrum (NSD)
def test_teslib_optfilt_nsd():
    """
    Test the teslib function for optimal filter  and NS (noise spectrum).
    This function performs the following steps:
    1. Defines input and output files.
    2. Sets the template energy to 6000 eV.
    3. Runs the teslib tool with the specified parameters (Noise Spectrum and no covar analysis).
    4. Compares the output file with a reference file.
    5. Removes the output files.
    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    lib_output = f"{outdir}/library_6keV_optfilt_ns.fits"
    noiseSpectrumFile = f"{reference}/noise_spectrum_reference_{date_config_file}.fits"
    
    # run teslib tool
    comm = (f"teslib Recordfile={file_input}"
            f" TesEventFile=myLibEvents.fits"
            f" monoenergy={template_eV}"
            f" LibraryFile={lib_output}"
            f" NoiseFile={noiseSpectrumFile}"
            f" addCOVAR=no"
            f" addINTCOVAR=no"
            f" addOFWN=no"
            f" XMLFile={xmlfile}"
    )
    
    output_teslib = run(comm, shell=True, capture_output=True)
    assert output_teslib.returncode == 0, f"teslib failed to run: {comm}"

    # compare output file with reference file
    ref_file = f"{reference}/library_6keV_optfilt_ns_reference_{date_config_file}.fits"
    assert os.path.exists(ref_file), f"Reference library file {ref_file} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {lib_output} {ref_file}'
    output_fdiff = run(comm, shell=True, capture_output=True)

    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), "Comparison failed: produced library does not match reference library"

    # remove output file
    os.remove(lib_output)
    os.remove("myLibEvents.fits")

# Create a test for teslib tool to be used only for optimal filter(s) and WEIGHTN for the noise
def test_teslib_optfilt_wn():
    """
    Test the teslib function for optimal filter  and WEIGHTN for the noise.
    This function performs the following steps:
    1. Defines input and output files.
    2. Sets the template energy to 6000 eV.
    3. Runs the teslib tool with the specified parameters (WEIGHTN and no covar analysis).
    4. Compares the output file with a reference file.
    5. Removes the output files.
    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    lib_output = f"{outdir}/library_6keV_optfilt_wn.fits"
    noiseSpectrumFile = f"{reference}/noise_spectrum_mat_reference_{date_config_file}.fits"

    # run teslib tool
    comm = (f"teslib Recordfile={file_input}"
            f" TesEventFile=myLibEvents.fits"
            f" monoenergy={template_eV}"
            f" LibraryFile={lib_output}"
            f" NoiseFile={noiseSpectrumFile}"
            f" addCOVAR=no"
            f" addINTCOVAR=no"
            f" addOFWN=yes"
            f" XMLFile={xmlfile}"
            )
    
    output_teslib = run(comm, shell=True, capture_output=True)
    assert output_teslib.returncode == 0, f"teslib failed to run: {comm}"

    # compare output file with reference file
    ref_file = f"{reference}/library_6keV_optfilt_wn_reference_{date_config_file}.fits"
    assert os.path.exists(ref_file), f"Reference library file {ref_file} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {lib_output} {ref_file}'
    output_fdiff = run(comm, shell=True, capture_output=True)

    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), "Comparison failed: produced library does not match reference library"

    # remove output file
    os.remove(lib_output)
    os.remove("myLibEvents.fits")

