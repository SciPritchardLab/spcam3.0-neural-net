import os
import fileinput
import sys
import re
pathOttmodels = "/work/08110/tg874091/stampede2/all_models_fkb_ott_et_al/"
pathICs = "/scratch/08110/tg874091/icSpawn32/run/"

#THIS FILE REQUIRES TWO INPUTS. 
# 1.) THE MODEL NUMBER YOU WISH TO USE 
# 2.) THE NAME OF THE INITIAL CONDITION FILE YOU'D LIKE TO USE

# PICK A MODEL NUMBER!
modelNum = int(sys.argv[1])
modLabel = "%05d" %modelNum

# PICK AN INITIAL CONDITION FILE!
initFile = sys.argv[2]
# making nickname for job
nickname = re.sub("^[^-]+", "", initFile)
nickname = re.sub("[^-]+$", "", nickname)
nickname = re.sub("-$", "", nickname)
nickname = re.sub("-", "_", nickname)
nickname = modLabel + nickname

rundir = "$SCRATCH/biasTrainingData/run_" + modLabel + "/" + nickname 

os.system(" ".join(["mkdir", "-p", rundir]))
os.system(" ".join(["cp", "-r", "baseline/*", rundir]))

#Copying in Jordan's NN to overwrite keras_matrices/model.txt
os.system(" ".join(["cp", pathOttmodels + modLabel + ".txt", rundir + "/keras_matrices/model.txt"]))

#Copying in initial condition file
os.system(" ".join(["cp", initFile, rundir]))

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
os.system(" ".join(["cp", "run.template", "run.csh"]))
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("RUNDIR", rundir), end='')
with fileinput.FileInput("run.csh", inplace=True) as file:
    for line in file:
        print(line.replace("INITCOND", nickname), end='')
os.system(" ".join(["mv", "run.csh", rundir]))

#RUNNING THE MODEL
os.system(" ".join(["sbatch", rundir + "/run.csh"]))
