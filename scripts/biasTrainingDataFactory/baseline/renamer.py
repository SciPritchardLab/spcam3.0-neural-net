import re
import os

#reading in file names
initial = os.popen(" ".join(["ls", "And*.i.*.nc"])).read()
h1Files = os.popen(" ".join(["ls", "And*.h1.*.nc"])).read().splitlines()

#extracting relevant dates
initial = re.sub("^[^-]+-", "", initial)
initial = re.sub("-[^-]+$", "", initial)
month = re.sub("-[0-9]+", "", initial)
day = re.sub("[0-9]+-", "", initial)
monthNum = int(re.sub("0", "", month))
dayNum = int(re.sub("0", "", day))

#renaming h1 files
for i in range(len(h1Files)):
    dayNum = dayNum + 1
    if dayNum > 30:
        dayNum = dayNum%30
        monthNum = monthNum + 1
        month = "%02d" %monthNum
    day = "%02d" %dayNum
    temp = re.sub("01-[0-9][0-9]-0", "01-" + day + "-0", h1Files[i])
    temp = re.sub("0000-[0-9][0-9]-", "0000-" + month + "-", temp)
    os.system(" ".join(["mv", h1Files[i], temp]))

