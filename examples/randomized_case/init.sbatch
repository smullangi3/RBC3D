#!/bin/bash

#Do the following before running
#change initcond.F90 layers variable
#change the #cpus and #cores here
#change log file name to reflect proc and layer count
#change num processes in srun to reflect the total execution threads

#run this in the case directory of RBC3D

#SBATCH --account=gts-coc-coda20
#SBATCH -N1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G
#SBATCH -t8:00:00
#SBATCH -q inferno
#SBATCH --mail-user=smullangi3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o "./D/init.log"

cd $SLURM_SUBMIT_DIR

cd D
rm x* r* w* 1* 2* 3*
cd ../

ml gcc mvapich2 netcdf-c netcdf-cxx netcdf-fortran fftw

make
srun ./initcond
