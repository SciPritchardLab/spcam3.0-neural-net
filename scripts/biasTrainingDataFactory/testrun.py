import os
import fileinput
import sys
pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al"
pathInitialConditions = "/scratch/08110/tg874091/icSpawn32/run"
# PICK A MODEL NUMBER!
modelNum = int(sys.argv[1])
modLabel = "%05d" %modelNum
rundir = "$SCRATCH/biasTrainingData/run_" + modLabel
os.system(" ".join(["mkdir", "-p", rundir]))
os.system(" ".join(["cp", "-r", "baseline/*", rundir]))
# PICK A MONTH!
initMonth = int(sys.argv[2])
initMonth = "%02d" %initMonth
# PICK A DAY!
initDay = int(sys.argv[3])
initDay = "%02d" %initDay
#Copying in Jordan's NN to overwrite keras_matrices/model.txt
os.system(" ".join(["cp", pathOttmodels + "/" + modLabel + ".txt", rundir + "/keras_matrices/model.txt"]))
#Making a custom atm_in file
os.system(" ".join(["cp", "atm_in.template", "atm_in"]))
with fileinput.FileInput("atm_in", inplace=True) as file:
    for line in file:
        print(line.replace("XXX", modLabel), end='')
with fileinput.FileInput("atm_in", inplace=True) as file:
    for line in file:
        print(line.replace("ZZ", initMonth), end='')
with fileinput.FileInput("atm_in", inplace=True) as file:
    for line in file:
        print(line.replace("YY", initDay), end='')
os.system(" ".join(["mv", "atm_in", rundir]))
#Making a custom run file
os.system(" ".join(["cp", "run.template", "run.csh"]))
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("RUNDIR", rundir), end='')
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("XXX", modLabel), end='')
os.system(" ".join(["mv", "run.csh", rundir]))
#Running the model
os.system(" ".join(["sbatch", rundir + "/run.csh"]))
