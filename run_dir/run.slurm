#!/bin/bash
#SBATCH -J cbrain-cam
#SBATCH -o myjob.o%j       
#SBATCH -e myjob.e%j       
#SBATCH -p normal          
#SBATCH -N 4              
#SBATCH -n 64              
#SBATCH -t 23:59:00        
#SBATCH --mail-user=s.rasp@lmu.de
#SBATCH --mail-type=all    
#SBATCH -A TG-ATM160016       
export I_MPI_COMPATIBILITY=4
export OMP_NUM_THREADS=16
cd /scratch/05488/tg847872/cloudbrain_ctrl_aquaplanet_03
ibrun -n 64 ./cam < atm_in

