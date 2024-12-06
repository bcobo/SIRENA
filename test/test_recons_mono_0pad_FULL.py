####################################################################################
# Runs the reconstruction of a file of monochromatic events with 'tesrecons' tool. 
# using the '0PAD' energy method (but no truncation) & FULL optimal filter to check    
# that they are consistent                                                              
####################################################################################

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


# non truncated 0PAD test
def test_recons_0pad8192_full():
    """
    Test the 'tesrecons' tool with the '0PAD (FULL length)' energy method.
    This function runs the 'tesrecons' tool with the '0PAD' energy method. 
    It takes a FITS file as input, performs reconstruction using the specified library file, 
    and generates an output FITS file with reconstructed events. 
    The output file (0PAD not truncated) is then compared with a reference file 
    (FULL optimal filter) to ensure correctness. 
    Finally, the output file is removed.

    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    file_input = f"{indir}/mono6keV_1000p_{date_config_file}.fits"
    file_output_0pad8192 = f"{outdir}/events_mono_6keV_0pad8192.fits"
    file_output_FULL = f"{outdir}/events_mono_6keV_optfilt8192.fits"
    library = f"{reference}/library_6keV_optfilt_ns_reference_{date_config_file}.fits"
        
    # run tesrecons tool with 0pad=FULL
    comm = (f"tesrecons Recordfile={file_input}"
            f" TesEventFile={file_output_0pad8192}"
            f" LibraryFile={library}"
            f" XMLFile={xmlfile}"
            f" clobber=yes"
            f" EnergyMethod=0PAD"
            f" prebuff_0pad=1000"
            f" flength_0pad=8192"
            f" OFStrategy=FIXED"
            f" filtEeV=6000"
            )
    
    output_tesrecons = run(comm, shell=True, capture_output=True)
    assert output_tesrecons.returncode == 0, f"tesrecons failed to run with command: {comm}"

    # run tesrecons tool with FULL
    comm = (f"tesrecons Recordfile={file_input}"
            f" TesEventFile={file_output_FULL}"
            f" LibraryFile={library}"
            f" XMLFile={xmlfile}"
            f" clobber=yes"
            f" EnergyMethod=OPTFILT"
            f" OFStrategy=FIXED"
            f" OFLength=8192"
            f" filtEeV=6000"
            )
    
    output_tesrecons = run(comm, shell=True, capture_output=True)
    assert output_tesrecons.returncode == 0, f"tesrecons failed to run with command: {comm}"

    # compare output files
    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {file_output_0pad8192} {file_output_FULL}'
    output_fdiff = run(comm, shell=True, capture_output=True)

    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), f"Comparison failed: reconstruction with 0PAD does not match reconstruction with FULL"
    
    # remove output file
    os.remove(file_output_0pad8192)
    os.remove(file_output_FULL)
