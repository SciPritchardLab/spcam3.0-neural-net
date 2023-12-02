# SPCAM3 with Neural Fortran Integration

### Building Neural Fortran
```
git clone https://github.com/jordanott/neural-fortran.git
cd neural-fortran
sh build_steps.sh
cd ..
```

### Compile SPCAM3
To compile SPCAM3 run:
`sh compile_run.sh`

To submit a single model to the queue run:
`sh comile_run.sh 00001 /scratch/06528/tg858273/SCAM3-NeuralFortran-Workflow-Copy/saved_models/output_weights/batch_1/00001.txt`

```
# compile_run.sh
# $1 : Case ID (number of model)
# $2 : Path/to/model.txt or ensemble file if $3 is not blank
# $3 : Blank if not using ensembles

if [ "$1" != "" -a  "$2" != "" ]; then
    echo "Your command line contains $# arguments, using first 2 as:"
    echo 'Output Directory = ' $1
    echo 'Path to txt file = ' $2
    caseid=$1
    txt_path=$2
else
    # if both arguments are blank: COMPILE

    $HOME/repositories/neural-net/scripts/build-spsld-stampede2-intel.csh
    ls -l $HOME/jordan-nn-build/run/

    exit 1
fi

# change to scratch
cds

mkdir $caseid
cd $caseid

# move cam file
cp $HOME/jordan-nn-build/run/cam ./

# moving weights and bias files
mkdir keras_matrices

# move normalization files
cp /scratch/06528/tg858273/SCAM3-NeuralFortran-Workflow-Copy/saved_models/007_32col_pnas_exact/inp* ./keras_matrices
cp /scratch/06528/tg858273/SCAM3-NeuralFortran-Workflow-Copy/saved_models/007_32col_pnas_exact/out* ./keras_matrices

# move model.txt file
cp $txt_path ./keras_matrices/model.txt

# copy initial conditions file
cp /work/05488/tg847872/stampede2/ic_bc/spinup_AndKua_aqua_SPCAM3.0.cam2.i.0000-12-02-00000.nc ./

# moving name list file - replace caseid for current run
cat $HOME/run_files/atm_in | sed "s@CASEID@$caseid@g" > atm_in

if [ "$3" != "" ]; then
    mkdir Models
    cp $txt_path Models

    num_threads=32
    num_nodes=8
else
    num_threads=1
    num_nodes=2
fi

cat /scratch/06528/tg858273/SCAM3-NeuralFortran-Workflow-Copy/run_files/run.csh | sed "s@CASEID@$caseid@g" | sed "s@NUM_OF_NODES@$num_nodes@g" | sed "s@NUM_OF_THREADS@$num_threads@g" > run.csh
sbatch run.csh
```

### Submitting multiple jobs at once

```
#!/bin/bash
# directory of models
model_dir=$1

echo $model_dir

for file_path in $model_dir*.txt; do
  file_name=$(basename $file_path)
  model_id="${file_name%.*}"

  sh compile_run.sh $model_id $file_path
done
```
