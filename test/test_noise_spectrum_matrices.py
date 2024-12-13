################################################################################
# Runs the generation of a noise spectrum with 'gennoisespec' (noise matrix)
# Then roughly compares this noise spectrum with a reference file.  
# WARNING: It takes several hours to run                                          
################################################################################

import os
from subprocess import run

# Define environmental variables
sixte = os.environ["SIXTE"]
sirena = os.environ["SIRENA"]
headas = os.environ["HEADAS"]

# if sirena or sixte are not initiated send message and stopa
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


# Create a test for gennoisespec tool to be tested with pytest
def test_gennoisespec_matrices():
    """
    Test the gennoisespec function.

    This function tests the gennoisespec function with noise matrix by performing the following steps:
    1. Define input and output files.
    2. Run the gennoisespec tool from the shell.
    3. Compare the output file with the reference file.
    4. Remove the output file.

    Parameters:
        None

    Returns:
        None
    """
    # Define input and output files
    ns_input = f"{indir}/noise_20240917.fits"
    ns_output = f"{outdir}/noise_spectrum_mat.fits"
    comm = (f" gennoisespec inFile={ns_input}"
            f" outFile={ns_output}"
            f" intervalMinSamples=8192"
            f" nintervals=1000"
            f" EnergyMethod=OPTFILT"
            f" rmNoiseIntervals=no"
            f" clobber=yes"
            f" weightMS=yes"
    )
    # run tool from shell
    output_gennoisespec = run(comm, shell=True, capture_output=True)
    assert output_gennoisespec.returncode == 0, f"gennoisespec failed to run"
    assert os.path.exists(ns_output), f"gennoisespec did not produce an output file"
    
    # compare output file with reference file
    ref_file = f"{reference}/noise_spectrum_mat_reference_{date_config_file}.fits"
    assert os.path.exists(ref_file), f"Reference noise spectrum file {ref_file} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {ns_output} {ref_file}'
    output_fdiff = run(comm, shell=True, capture_output=True)
    
    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), f"Comparison failed: produced noise spectrum does not match reference noise spectrum"
    
    # remove output file
    os.remove(ns_output)
