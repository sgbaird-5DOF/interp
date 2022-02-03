#!/bin/bash
module load matlab
#fn='get_cmd_test'
fn='randOctParityData'
matlab -nodisplay -nosplash -r "cd ../../../; addpath(genpath('.')); $fn; exit"
