#!/bin/bash

#SBATCH --mail-user=sterling.baird@icloud.com
#SBATCH --mail-type=TIME_LIMIT,ARRAY_TASKS

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#---------------------Call to MATLAB Code--------------
mkdir ~/.matlab/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID
unset TZ
module load matlab/r2018b

tid=$SLURM_ARRAY_TASK_ID

if [ "$for_type" = "parallel" ]; then
    pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '~/.matlab/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID'; parpool(pc,$inf_cores); parfor_Q = 1;"
    echo "pc_opts: $pc_opts"
else
    pc_opts="parfor_Q = 0;"
    echo "pc_opts: $pc_opts"
fi

echo "par_filepath: $par_filepath
jid: $jid
tid: $tid
dim: $dim
for_type: $for_type
"
walltime_name="walltime$jid"
walltime=${!walltime_name}

echo "walltime: $walltime"

matlab -nodisplay -nosplash -r "clear all; $pc_opts; jspace_driver($par_filepath,'$dim',$tid,parfor_Q)"

rm -rf ~/.matlab/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID

# extra code
#export PATH="$PATH:usr/sbaird9/1DOF_functional"
# export SLURM_ARRAY_JOB_ID=1000
# export SLURM_ARRAY_TASK_ID=1000
