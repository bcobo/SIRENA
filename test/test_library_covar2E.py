###################################################################################
# Runs the generation of libraries for covariance analysis (COVAR & INTCOVAR)
# adding only 2 energies (less memory and running time)                             
# Then roughly compares these libraries with  reference library files.            
# Generation of library with covar information is time-consuming: optional test   
###################################################################################

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

print("Warning: This test takes several hours to run")

outdir = "./data/output"
indir = "./data/input"
reference = "./data/reference"
# read config file date form file config_version.txt
with open(f"{indir}/config_version.txt", "r") as file:
    date_config_file = file.read().strip()

# Find name of (unique) xml file in indir directory
xmlfile = glob.glob(f"{indir}/config*.xml")[0]

# Create a test for teslib tool with all caovariance matrices
def test_teslib_covar2E():
    """
    Test the teslib function.
    This function performs the following steps:
    1. Defines input and output files.
    2. Sets the templates energy to 5000 eV and 7000 eV.
    3. Runs the teslib tool with the specified parameters.
    4. Compares the output file with a reference file.
    5. Removes the output files.
    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    file1_input = f"{indir}/mono5keV_1000p_{date_config_file}.fits"
    file2_input = f"{indir}/mono7keV_1000p_{date_config_file}.fits"
    energies = [5000., 7000.]
    files_input = [file1_input, file2_input]
    lib_output = f"{outdir}/library_5keV_7keV_all.fits"
    noiseSpectrumFile = f"{reference}/noise_spectrum_reference_{date_config_file}.fits"
    
    # run teslib tool in a loop for each of the input energies
    for ie in range(len(files_input)):
        file_input = files_input[ie]
        print(file_input)
        template_eV = energies[ie]
        comm = (f"teslib Recordfile={file_input}"
                f" TesEventFile=myLibEvents.fits"
                f" monoenergy={template_eV}"
                f" LibraryFile={lib_output}"
                f" NoiseFile={noiseSpectrumFile}"
                f" addCOVAR=yes"
                f" addINTCOVAR=yes"
                f" addOFWN=no"
                f" XMLFile={xmlfile}"
                )
        
        output_teslib = run(comm, shell=True, capture_output=True)
        assert output_teslib.returncode == 0, f"teslib failed to run: {comm}"
        
    # compare output file with reference file
    ref_file = f"{reference}/library_5keV_7keV_all_reference_{date_config_file}.fits"
    assert os.path.exists(ref_file), f"Reference library file {ref_file} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {lib_output} {ref_file}'
    output_fdiff = run(comm, shell=True, capture_output=True)
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), "Comparison failed: produced library does not match reference library"

    # remove output file
    os.remove(lib_output)
    os.remove("myLibEvents.fits")