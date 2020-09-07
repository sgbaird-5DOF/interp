# octonion-mesh
 code related to meshing and interpolation of grain boundary octonions

## Dependencies
### MATLAB Version
MATLAB 2019b or higher (mainly for the "arguments" syntax checking at
the beginning of functions, which is used extensively throughout)
### Toolboxes
-Statistics and Machine Learning Toolbox (for Gaussian Process Regression: fitrgp(), fitrgp.predict())
-Symbolic Math Toolbox (optional, for numStabBary.m)

## Usage
### Step0: download the code
`git clone --recurse-submodules https://github.com/sgbaird/octonion-inference.git`
See [cloning a repository](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) and[git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for more information or other options.

### Step1: navigate to [octonion-inference/code/](code/)
`cd octonion-inference/code/`

### Step3: open MATLAB and call [run.m](code/run.m) or [interp5DOF_test.m](code/interp5DOF_test.m)
`matlab; run`

## Accessing functions via addpathdir()
dir() and addpath() commands are used to locate functions in subfolders of the current working directory via a custom function [code/addpathdir.m](https://github.com/sgbaird/octonion-inference/blob/master/code/addpathdir.m). This could give anomalous behavior if the directory structure is changed such that filenames are non-unique in sub-folders of the parent folder where addpathdir() gets called, or if files with the same name are present elsewhere on the user's MATLAB path. This is also the only function that is shadowed (to my knowledge) within this repository ([by octonion-mesh](code/octonion-mesh/)); however, the functionality is fairly basic, and I don't anticipate any changes to the functionality. In other words, it shouldn't matter which one gets called.

## Getting started
Look at [code/interp5DOF.m](code/interp5DOF.m), which is a top-level function for creating a mesh, importing/generating data, triangulating a mesh, identifying the intersecting facet for datapoints, and finally, computing an interpolation. Also consider looking at [code/run.m](code/run.m) which is a script that contains more options, but may be less portable than the interp5DOF() function. interp5DOF.m can be called in other functions/scripts to produce interpolation results using 5DOF misorientation/boundary plane normal pairs and grain boundary property values. It was written with loosely similar input/output structure to the MATLAB built-in function interpn() involving input points, input values, query points, and output values.

## Test functions
Most functions have a corresponding "test" function (e.g. hsphext_subdiv.m -> hsphext_subdiv_test.m) which gives simple usage example(s). These are useful for debugging, visualizations, and understanding the functions without having to do a full run which could be time-consuming. The various test functions generally run to completion within a few seconds, and the parameters may be changed freely (e.g. dimension, number of points, etc.) where applicable.

## parfor loops
Parfor loops are used by default where there is potential for significant speed-up. A parfor-compatible text progress bar is encoded into each of these. Using disp() or fprintf() inside the parfor loop (aside from what's already inside the nested text progress bar function) may cause odd behavior on the command line, but should not affect the integrity of the code execution. Because these contain nested functions, in order to add new variables to the "static workspace", you can either assign variables to "ans" (a special variable that is still accessible), output statements directly to the command line terminal, or comment the nested function, nUpdateProgress(). If the parallel computing toolbox is not installed, the parfor loops will execute as regular for loops. If the parallel computing toolbox is installed, and you want to use only a single core, start a parallel pool with only one core before running any of the functions. The loop will still do all of its internal parfor executions. A parfor loop with a single core and parallel computing toolbox should not run any slower than a regular for loop when they are contained within functions. A parfor loop executed within a script, however, is likely to result in significant slow-down. To debug within a parfor loop, simply change it to a for loop temporarily. MATLAB "find files" (Alt,E,F,F) searching for the keywords "parfor compatible" should help in keeping track of which ones have been changed.
