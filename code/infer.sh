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
pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '~/.matlab/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID'; parpool(pc,$inf_cores);"
echo "pc_opts: $pc_opts"
fi

echo "load_filepath: $load_filepath
jid: $jid
tid: $tid
for_type: $for_type
"
walltime_name="walltime$jid"
walltime=${!walltime_name}

echo "walltime: $walltime"

matlab -nodisplay -nosplash -r "clear all; $pc_opts; sse_inference($load_filepath, $jid, $tid)"

rm -rf ~/.matlab/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID


# extra code
#export PATH="$PATH:usr/sbaird9/1DOF_functional"
# export SLURM_ARRAY_JOB_ID=1000
# export SLURM_ARRAY_TASK_ID=1000
