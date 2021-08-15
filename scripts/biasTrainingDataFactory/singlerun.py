import os
import fileinput
import sys
pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al/"
pathICs = "/scratch/08110/tg874091/icSpawn32/run/"

# PICK A MODEL NUMBER!
modelNum = int(sys.argv[1])
modLabel = "%05d" %modelNum
rundir = "$SCRATCH/biasTrainingData/run_" + modLabel

#This should be commented out if using multi run as this is done in that python script.
#os.system(" ".join(["mkdir", "-p", rundir]))
#os.system(" ".join(["cp", "-r", "baseline/*", rundir]))


# PICK AN INITIAL CONDITION FILE!
initFile = sys.argv[2]

#Copying in Jordan's NN to overwrite keras_matrices/model.txt
#os.system(" ".join(["cp", pathOttmodels + modLabel + ".txt", rundir + "/keras_matrices/model.txt"]))
#Copying in initial condition file
os.system(" ".join(["cp", pathICs + initFile, rundir]))
#Making a custom atm_in file
os.system(" ".join(["cp", "atm_in.template", "atm_in"]))
with fileinput.FileInput("atm_in", inplace=True) as file:
    for line in file:
        print(line.replace("XXX", modLabel), end='')
with fileinput.FileInput("atm_in", inplace=True) as file:
    for line in file:
        print(line.replace("initialConditionFileName", initFile), end='')
os.system(" ".join(["mv", "atm_in", rundir]))
#Making a custom run file
#This code is not necessary when using the multirun.py script
#os.system(" ".join(["cp", "run.template", "run.csh"]))
#with fileinput.FileInput("run.csh", inplace=True) as file:
#    for line in file:
#        print(line.replace("RUNDIR", rundir), end='')
#with fileinput.FileInput("run.csh", inplace=True) as file:
#    for line in file:
#        print(line.replace("XXX", modLabel), end='')
#os.system(" ".join(["mv", "run.csh", rundir]))

#RUNNING THE MODEL
os.system(" ".join(["sbatch", rundir + "/run.csh"]))
