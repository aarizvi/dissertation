#!/bin/bash

#SBBATCH --time=24:00:00

#SBATCH --nodes=14

#SBATCH --ntasks-per-node=4

#SBATCH --mem-per-cpu=5000

#SBATCH --job-name=myjobname

#SBATCH --output=myjob.out

#SBATCH --error=myjoberr.err

#SBATCH --partition=general-compute

#SBATCH --mail-user=abbasriz@buffalo.edu

#SBATCH --mail-type=ALL

#SBATCH --mail-type=END

#Get date and time

tstart=`date`

echo "###### start time:"$tstart

cd /projects/rpci/lsuchest/abbasriz/candidate_gene/result_files_cg/res.directories/
          
echo "Run program"
          
module load R 

R CMD BATCH independent_results.R

echo "program finished"

echo "All Done!"

tend=`date`

echo "###### end time: $tend"  
