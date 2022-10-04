import os
import fileinput
import sys
import re


# PICK AN INITIAL CONDITION FILE!
initFile = sys.argv[1]

# making nickname for job
nickname = re.sub("^[^-]+", "", initFile)
nickname = re.sub("[^-]+$", "", nickname)
nickname = re.sub("-$", "", nickname)
nickname = "training" + re.sub("-", "_", nickname)

rundir = "$SCRATCH/dataProcessingFactory/cleanDataFolder/" + nickname 
copydir = "$SCRATCH/spreadTrainingData/" + nickname + "/"

icFolders = os.popen(" ".join(["ls", copydir])).read().splitlines()
slimFolders = [x for x in icFolders if "h1" in x][10:]


os.system(" ".join(["mkdir", "-p", rundir]))
os.system(" ".join(["cp", "-r", "baseline/*", rundir]))

for h1File in slimFolders:
    os.system(" ".join(["cp", copydir + h1File, rundir]))

os.system(" ".join(["cp", "minibatch.template", "custombatch.template"]))
with fileinput.FileInput("custombatch.template", inplace=True) as file:
    for line in file:
        print(line.replace("TRAINPATH", nickname), end='')
os.system(" ".join(["mv", "custombatch.template", rundir]))

#RUNNING THE MODEL

os.system(" ".join(["sbatch", rundir + "/custombatch.template"]))
