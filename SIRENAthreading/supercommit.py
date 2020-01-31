"""
This script changes adds one to the last number in versionSIRENA 
(previous to commit) to GitHub

"""

import os
from datetime import datetime
filein = open('./versionSIRENA.h','rt')
fileout = open('./versionSIRENA2.h','wt')


for line in filein:
    if 'define SIRENA_VERSION "' in line:
        lp=line.rfind('.')  # last dot
        lc=line.rfind('"')  # last "
        nv=int(line[lp+1:lc])+1
        newver=line[:lp+1]+str(nv)+'"\n'
        fileout.write(newver)
    elif "DATE:" in line:
        newdate = datetime.now().strftime("%Y/%m/%d, %H:%M:%S")
        fileout.write("//   DATE: " + newdate + "\n")
    else:
        fileout.write(line)
filein.close()
fileout.close()
os.rename('versionSIRENA2.h','versionSIRENA.h')
print("SIRENA version updated to ",newver)

