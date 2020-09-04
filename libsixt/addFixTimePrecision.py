#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday May 16 2019

@author: cobo
"""

from astropy.io import fits
from random import uniform
import shlex
from subprocess import check_call, STDOUT
import sys

def changeTimeColumn(piximpactFile,eventsFile,timePrecision):
#    print(piximpactFile)
    imp = fits.open(piximpactFile, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']
        
    random_numbers = [uniform(-float(timePrecision),float(timePrecision)) for _ in range(len(impTimes))]
#    print(random_numbers)
#    print(random_numbers[0],random_numbers[1],random_numbers[2])
#    print(impTimes[0])
    impTimes = impTimes + random_numbers
#    print(impTimes[0])
    impTab['TIME'] = impTimes
    
    imp.writeto(piximpactFile + '.' + timePrecision)
    
    imp.close()
    
    # Remove rows because number of pulses in events file is differents from number of photon in impacts file
    #eventsFile = piximpactFile.replace('piximpact','fits')
    evt = fits.open(eventsFile, memmap=True)
    num_evts = evt[1].header["NAXIS2"]
    if (len(impTimes) != num_evts):
        num_rowstodelete_end = len(impTimes)-num_evts-1
    else:
        num_rowstodelete_end = 0
    evt.close()
#    print(len(impTimes))
#    print(num_evts)
#    print(num_rowstodelete_end)
    if (num_rowstodelete_end != 0):
            # Last rows
        comm = "fdelrow infile=" + piximpactFile + '.' + timePrecision + "[1] firstrow=" + str(len(impTimes)-num_rowstodelete_end+1) + " nrows=" + str(num_rowstodelete_end) + " confirm=no proceed=yes"
    #    print(comm)
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
    #        print("     Last rows removed in ", piximpactFile)
        except RuntimeError:
            print("Error running ", comm, " to remove final rows in ", piximpactFile)
            raise
            # First row
        comm = "fdelrow infile=" + piximpactFile + '.' + timePrecision + "[1] firstrow=1 nrows=1 confirm=no proceed=yes"
    #    print(comm)
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
    #        print("     First row removed in ", piximpactFile)
        except RuntimeError:
            print("Error running ", comm, " to remove first row in ", piximpactFile)
            raise
        
def addFixTimePrecision(piximpactFile,eventsFile,timePrecision):
    # It creates a new .piximpact.timePrecision file with the TIME column changed and the proper rows deleted
    # timePrecision (sec) => photon_time(PIXIMPACT)+-timePrecision
    # Applied only to a file or to several files if name file starts with '@'
    #
    # python addFixTimePrecision.py --piximpactFile="@listadePIXIMPACT.txt" --eventsFile="@listadeEVENTS.txt" --timePrecision=0.05e-6
    # python addFixTimePrecision.py --piximpactFile=mono-1000.piximpact --eventsFile=mono-1000.fits --timePrecision=0.05e-6
    
#    print("piximpactFile: ",piximpactFile)
#    print("eventsFile: ",eventsFile)
#    print("timePrecision: ",timePrecision)
    
    if ((((piximpactFile[0] == '@') and (eventsFile[0] != '@'))) or (((piximpactFile[0] != '@') and (eventsFile[0] == '@')))):
        sys.exit("Error: Both input files must be ASCII files (starting with @) or single FITS files")
    
    # Check if piximpact file starts with '@'
    if (piximpactFile[0] == '@'):
        # To know the number of 
        piximpactFile=piximpactFile[1:len(piximpactFile)]  # Delete the '@'
        listPIXIMPACTFile=open(piximpactFile,"r")
        file_lines = listPIXIMPACTFile.readlines()
        numfiles = len(file_lines)
        last_line = file_lines[len(file_lines)-1]
        if (last_line == '\n'):
            numfiles = numfiles-1
#        print("last line: ",last_line)
#        print("Number of files: ",numfiles)
        listPIXIMPACTFile.close()
        
        eventsFile=eventsFile[1:len(eventsFile)]  # Delete the '@'
        listEVENTSFile=open(eventsFile,"r")
        eventsfile_lines = listEVENTSFile.readlines()
        events_numfiles = len(eventsfile_lines)
        eventslast_line = eventsfile_lines[len(eventsfile_lines)-1]
        if (eventslast_line == '\n'):
            events_numfiles = events_numfiles-1
        if (events_numfiles != numfiles):
            sys.exit("Error: Both input ASCII files must have the same number of lines")
        listEVENTSFile.close()
            
    else:
        numfiles = 1
        #print("Number of files: ",numfiles)
        
    for i in range(0,numfiles):
        if numfiles == 1:
            changeTimeColumn(piximpactFile,eventsFile,timePrecision)
        else: 
            changeTimeColumn(file_lines[i][0:len(file_lines[i])-1],eventsfile_lines[i][0:len(eventsfile_lines[i])-1],timePrecision)
            
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Add a fix time precision to the photons time (in a new PIXIMPACT file)',
                                     prog='addFixTimePrecision ')
    parser.add_argument('--piximpactFile', help='file with impact times or ASCII file with several impact files (if it starts with @)', required=True)
    parser.add_argument('--eventsFile', help='file with events or ASCII file with several events files (if it starts with @)', required=True)
    parser.add_argument('--timePrecision', help='time precision', required=True)

inargs = parser.parse_args()

addFixTimePrecision(piximpactFile=inargs.piximpactFile, eventsFile=inargs.eventsFile,timePrecision=inargs.timePrecision)
