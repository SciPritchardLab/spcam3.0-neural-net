#!/bin/bash
#SBATCH -J spcam3.0
#SBATCH -o myjob.o%j       
#SBATCH -e myjob.e%j       
#SBATCH -p normal
#SBATCH -N 4
#SBATCH -n 64
#SBATCH -t 00:10:00        
##SBATCH --mail-user=fernando.iglesias-suarez@dlr.de
##SBATCH --mail-type=all
#SBATCH -A TG-ATM190002
export I_MPI_COMPATIBILITY=4
export OMP_NUM_THREADS=16
ibrun -n 64 ./cam < atm_in > logfile

