import os
import sys
pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al/"
icPath = "$SCRATCH/icSpawn32/run/"
icFiles = os.popen(" ".join(["ls", icPath + "And*.i.*.nc"])).read().splitlines()


# THIS FILE REQUIRES YOU TO INPUT TWO NUMBERS. 
# 1.) THE MODEL NUMBER YOU WISH TO USE 
# 2.) THE FREQUENCY OF INITIAL CONDITION FILES USED. (ie every 7th, 8th, etc)

# PICK A MODEL!
modelNum = int(sys.argv[1])
modLabel = "%05d" %modelNum
rundir = "$SCRATCH/biasTrainingData/run_" + modLabel + "/"
os.system(" ".join(["mkdir", "-p", rundir]))
os.system(" ".join(["cp", "-r", "baseline/*", rundir]))
modelNum = str(modelNum)
os.system(" ".join(["cp", pathOttmodels + modLabel + ".txt", rundir + "/keras_matrices/model.txt"]))
#Making a custom run file
os.system(" ".join(["cp", "run.template", "run.csh"]))
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("RUNDIR", rundir), end='')
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("XXX", modLabel), end='')
os.system(" ".join(["mv", "run.csh", rundir]))

#PICK A FREQUENCY
everyWhich = int(sys.argv[2])
for i in range(len(icFiles)):
    if (i+1)%everyWhich==0: 
        os.system(" ".join(["python", "singlerun.py", modelNum, icFiles[i]]))
        os.system(" ".join(["rm", rundir, "atm_in"]))
