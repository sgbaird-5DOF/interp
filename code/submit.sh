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


module load matlab/r2019b
#pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '$TMPDIR'; parpool(pc,$cores);"
pc_opts="pc=parcluster('local'); pc.JobStorageLocation = '$TMPDIR'; parpool(pc,$cores,'SpmdEnabled',false)" #spmd parameter set to false so loop continues even if worker aborts
echo $pc_opts
matlab -nodisplay -nosplash -r "$pc_opts; run; exit"

echo "Cleaning up temporary directory at end of script, meaning that the job exited cleanly"
rm -rfv $TMPDIR
