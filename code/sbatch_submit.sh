#!/bin/bash
module load matlab/r2019b
matlab -nodisplay -nosplash -r "$pc_opts; get_cmd_test; exit"

