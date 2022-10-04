import os
import sys
import re
import fileinput

icPath = "$WORK2/mergedInitialConditions/"
icFiles = os.popen(" ".join(["ls", icPath + "And*.i.*.nc"])).read().splitlines()


#PICK A FREQUENCY!
everyWhich = int(sys.argv[1])
for i in range(len(icFiles)):
    if (i+1)%everyWhich==0:
        os.system(" ".join(["python", "singlerun.py", icFiles[i]]))

