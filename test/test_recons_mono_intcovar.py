####################################################################################
# Runs the reconstruction of file of monochromatic events with 'tesrecons' tool.   
# Use INTCOVAR energy method. 
# Check that results are compatible with reference file.                                
####################################################################################

import os
from subprocess import run, PIPE

# Define environmental variables
sixte = os.environ["SIXTE"]
sirena = os.environ["SIRENA"]
headas = os.environ["HEADAS"]

# if sirena or sixte are not initiated send message and stop
assert sirena != "", "Please set the environmental variable SIRENA and initiate the software using (for sh-like shell): . $SIRENA/bin/sirena-install.sh"
assert sixte != "", "Please set the environmental variable SIXTE and initiate the software using (for sh-like shell): . $SIXTE/bin/sixte-install.sh"
assert headas != "", "Please set the environmental variable HEADAS and initiate the software using (for sh-like shell): . $SIXTE/headas-init.sh"

outdir = "/data/output"
indir = "/data/input"
reference = "/data/reference"
xmlfile = f"./{indir}/1pix_lpa2.5a_fll_SIRENAintegration.xml"

def test_recons_mono_intcovar():
    """
    Test the 'tesrecons' tool with the INTCOVAR reconstruction method.
    It takes a FITS file as input, performs reconstruction using the specified library file, 
    and generates an output FITS file with reconstructed events. 
    The output file is then compared with a reference file 
    Finally, the output file is removed.

    Parameters:
    None
    Returns:
    None
    """
   
    # Define input and output files
    file_input = f".{indir}/mono6.0keV_1000p.fits"
    file_output = f".{outdir}/events_6.0keV_intcovar.fits"
    library = f".{reference}/library_5keV_7keV_all_reference.fits"
        
    # run tesrecons tool with FULL
    comm = (f"tesrecons Recordfile={file_input}"
            f" TesEventFile={file_output}"
            f" LibraryFile={library}"
            f" XMLFile={xmlfile}"
            f" clobber=yes"
            f" EnergyMethod=INTCOVAR"
            f" OFStrategy=BYGRADE"
            )
    
    output_tesrecons = run(comm, shell=True, capture_output=True, text=True)
    assert output_tesrecons.returncode == 0, f"tesrecons failed to run with command '{comm}' with problem: {output_message_tesrecons.stdout}"

    # compare output files
    file_reference = f".{reference}/events_mono_6.0keV_intcovar_reference.fits"
    assert os.path.exists(file_reference), f"Reference file {file_reference} does not exist"

    comm = f'fdiff caldsum=yes cmpdata=no exclude="CREADATE, SIRENAV" {file_output} {file_reference}'
    output_fdiff = run(comm, shell=True, capture_output=True)
    
    # check if the comparison was successful
    # Look for the string "End of file comparison:  0 differences were found" in the output of fdiff
    success_message = "End of file comparison:  0 differences were found"
    assert success_message in str(output_fdiff.stdout), f"Comparison failed: reconstructed file does not match reference file"
    
    # remove output file
    os.remove(file_output)
