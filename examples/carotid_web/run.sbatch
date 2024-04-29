#!/bin/bash

#Do the following before running
#change initcond.F90 layers variable
#change the #cpus and #cores here
#change log file name to reflect proc and layer count
#change num processes in srun to reflect the total execution threads

#run this in the case directory of RBC3D

#SBATCH --account=gts-coc-coda20
#SBATCH -N8 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4GB
#SBATCH -t8:00:00
#SBATCH -q inferno
#SBATCH --mail-user=smullangi3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./D/carotid_vertweb.log"
#SBATCH --open-mode=append

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/storage/home/hcoda1/6/smullangi3/lib/

cd $SLURM_SUBMIT_DIR

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

timeout 7h srun tube 
EX_CODE=$?

# move files to scratch, except for log file and restart.LATEST
cd D
cp * /storage/home/hcoda1/6/smullangi3/scratch/carotid_vertweb/
rm r* x* w*
cp /storage/home/hcoda1/6/smullangi3/scratch/carotid_vertweb/restart.LATEST.dat ./restart.LATEST.dat
cd ../

if [[ ${EX_CODE} == 124 ]]; then 
  sbatch run.sbatch
fi