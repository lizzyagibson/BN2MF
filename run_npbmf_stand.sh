#!/bin/bash
#$ -cwd -S /bin/bash
#$ -l mem=5G
#$ -l time=:10000:
#$ -M eag2186@cumc.columbia.edu

$MODULESHOME/init/bash
module load matlab/2018a

clear

matlab -nojvm -nodisplay -nosplash  -singleCompThread -nodesktop  -logfile "/ifs/scratch/msph/ehs/eag2186/log-mat{$SGE_TASK_ID}"  -r "run /ifs/scratch/msph/ehs/eag2186/nmf_sim_hpc.m"


