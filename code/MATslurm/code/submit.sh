#!/bin/bash

TMPDIR=/tmp/$SLURM_JOB_ID

# prepare function to clean up directory if job is killed, etc
# NOTE: THIS IS NOT EXECUTED UNTIL THE SIGNAL IS CALLED.
cleanup_scratch() {
    echo "Cleaning up temporary directory inside signal handler, meaning I either hit the walltime, or deliberately deleted this job using scancel"
    rm -rfv $TMPDIR
    echo "Signal handler ended at:"
    date
    exit 1
}

#Now, associate this function with the signal that is called when the job is killed, or hits its walltime
trap 'cleanup_scratch' TERM

#Now, initial setup of temporary scratch directory
echo "Creating scratch directory at $TMPDIR"
mkdir -pv $TMPDIR 2>&1

#PUT CODE TO COPY YOUR DATA INTO $TMPDIR HERE IF NECESSARY
#DO WHATEVER YOU NEED TO DO TO GET YOUR SOFTWARE TO USE $TMPDIR. THIS WILL DEPEND ON THE SOFTWARE BEING USED
#PUT CODE TO RUN YOUR JOB HERE

tid=$SLURM_ARRAY_TASK_ID

echo "parpath: $parpath
jid: $jid
tid: $tid
"
walltime_name="walltime$jid"
walltime=${!walltime_name}

echo "walltime: $walltime"


module load matlab
pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '$TMPDIR'; parpool(pc,$cores,'SpmdEnabled',logical($spmdQ));"
echo $pc_opts
#matlab -nodisplay -nosplash -r "$pc_opts; run; exit"
#matlab -nodisplay -nosplash -r "$pc_opts; randOctParityData; randOctParityPlot; exit"

matlab -nodisplay -nosplash -r "$pc_opts; addpath(genpath('.')); exec_combs($parpath, $jid, $tid)"
echo "Cleaning up temporary directory at end of script, meaning that the job exited cleanly"
rm -rfv $TMPDIR

git add figures/


#extra code

#pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '$TMPDIR'; parpool(pc,$cores,'SpmdEnabled',false)" #spmd parameter set to false so loop continues even if worker aborts
#pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '$TMPDIR'; parpool(pc,12,'SpmdEnabled',false)" #spmd parameter set to false so loop continues even if worker aborts
