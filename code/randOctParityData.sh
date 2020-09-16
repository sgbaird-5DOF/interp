#!/bin/bash
module load matlab/r2019b
fn='randOctParityData'
matlab -nodisplay -nosplash -r "$fn; exit"
