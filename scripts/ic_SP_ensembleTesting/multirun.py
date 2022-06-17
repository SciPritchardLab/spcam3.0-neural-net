import os
import sys
import re
import fileinput
#pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al/"
icPath = "$WORK2/mergedInitialConditions/"
icFiles = os.popen(" ".join(["ls", icPath + "And*.i.*.nc"])).read().splitlines()



# THIS SCRIPT REQUIRES TWO INPUTS. 
# 1.) THE MODEL NUMBER YOU WISH TO USE 
# 2.) THE FREQUENCY OF INITIAL CONDITION FILES USED. (ie every 7th, 8th, etc)

# PICK A MODEL!
#modelNum = int(sys.argv[1])
#modelNum = str(modelNum)

#PICK A FREQUENCY!
everyWhich = int(sys.argv[1])
for i in range(len(icFiles)):
    if (i+1)%everyWhich==0:
        os.system(" ".join(["python", "singlerun.py", icFiles[i]]))

