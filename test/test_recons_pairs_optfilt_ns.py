####################################################################################
# Runs the reconstruction of a file of pairs of events with 'tesrecons' tool
# using optimal filter and Noise spectrum                                      
# Check results are compatible with reference file.                                
####################################################################################

import os
from subprocess import run, PIPE
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


def test_recons_pairs_optfilt_ns():
    """
    Test the 'tesrecons' tool with the OPTFILT energy method and NS noise spectrum.
    This function tests the 'tesrecons' tool with the OPTFILT energy method.
    It takes a FITS file as input, performs reconstruction using the specified library file, 
    and generates an output FITS file with reconstructed events. 
    The output file is then compared with a reference file to ensure correctness. 
    Finally, the output file is removed.

    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    file_input = f"{indir}/xifusim_pairs_3keV_6keV_sep3500samples_{date_config_file}.fits"
    file_output = f"{outdir}/events_pairs_3keV_6keV_sep3500samples_optfilt_ns.fits"
    library = f"{reference}/library_6keV_optfilt_ns_reference_{date_config_file}.fits"
        
    # run tesrecons tool with OPTFILT (+ NSD)
    comm = (f"tesrecons Recordfile={file_input}"
            f" TesEventFile={file_output}"
            f" LibraryFile={library}"
            f" XMLFile={xmlfile}"
            f" clobber=yes"
            f" EnergyMethod=OPTFILT"
            f" OFStrategy=BYGRADE"
            f" filtEeV=6000"
            f" OFNoise=NSD"
            )
    
    output_tesrecons = run(comm, shell=True, capture_output=True, text=True)
    assert output_tesrecons.returncode == 0, f"tesrecons failed to run with command: {comm}"

    # compare output files
    file_reference = f"{reference}/events_pairs_3keV_6keV_sep3500samples_optfilt_ns_reference_{date_config_file}.fits"
    assert os.path.exists(file_reference), f"Reference file {file_reference} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {file_output} {file_reference}'
    output_fdiff = run(comm, shell=True, capture_output=True)
        
    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), f"Comparison failed: reconstructed file does not match reference file"
    
    # remove output file
    os.remove(file_output)
