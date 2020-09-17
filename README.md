# Five Degree-of-Freedom (5DOF) Interpolation
 code related to meshing and interpolation of grain boundaries by representing 5DOF of grain boundaries as octonions and forming a closed mesh.
 
 See
 * [GB_octonion_code](https://github.com/ichesser/GB_octonion_code)
 * Chesser, I., Francis, T., De Graef, M., & Holm, E. A. (2020). Learning the Grain Boundary Manifold: Tools for Visualizing and Fitting Grain Boundary Properties. Acta Materialia. https://doi.org/10.2139/ssrn.3460311
 * Francis, T., Chesser, I., Singh, S., Holm, E. A., & De Graef, M. (2019). A geodesic octonion metric for grain boundaries. Acta Materialia, 166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034

## Dependencies
### MATLAB Version
MATLAB R2019b or higher (mainly for the [arguments ... end syntax checking](https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html) at
the beginning of functions, which is used extensively throughout). For users of R2007a - R2019a, I suggest removing the arguments ... end syntax for any functions that use this and replacing it with corresponding [inputParser()](https://www.mathworks.com/help/matlab/ref/inputparser.html) and [varargin](https://www.mathworks.com/help/matlab/ref/varargin.html) code to deal with variable input arguments, default parameter values, and repeating arguments. Alternatively, you could remove the arguments ... end syntax lines for each function and update every place that the function is called so that all input arguments are specified. Open up an issue if you need more details on this.

### Toolboxes
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html) (for Gaussian Process Regression: [fitrgp()](https://www.mathworks.com/help/stats/fitrgp.html), [fitrgp.predict()](https://www.mathworks.com/help/stats/compactregressiongp.predict.html))
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) (optional, but for fitrgp() may need to change `hyperopts = struct('UseParallel',true,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` to `hyperopts = struct('UseParallel',false,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` in [interp5DOF.m](code/interp5DOF.m) under "method-specific interpolation" section --> 'gpr' case.)
- [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html) (optional, for [numStabBary.m](code/numStabBary.m))
- [Signal Processing Toolbox](https://www.mathworks.com/products/symbolic.html)

### Files
See [File dependencies](https://github.com/sgbaird/octonion-mesh/blob/master/README.md#file-dependencies)

## Usage
Linux-specific commands given.

### Step0: download the code
See [cloning a repository](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) and [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for more information or other options such as using GitHub Desktop (Windows, Linux, etc.) or downloading a .zip file.

`git clone --recurse-submodules https://github.com/sgbaird/octonion-inference.git`

### Step1: navigate to [interp-5DOF/code/](code/)
`cd interp-5DOF/code/`

### Step3: open MATLAB and call [interp5DOF_test.m](code/interp5DOF_test.m) or [run.m](code/run.m)
`matlab`

\>\> `interp5DOF_test`

## Accessing functions via addpathdir()
dir() and addpath() commands are used to locate functions in subfolders of the current working directory via a custom function [addpathdir.m](code/addpathdir.m). This could give anomalous behavior if the directory structure is changed such that filenames are non-unique in sub-folders of the parent folder where addpathdir() gets called, or if files with the same name are present elsewhere on the user's MATLAB path. This is also the only function that is shadowed (to my knowledge) within this repository (it's shadowed by [octonion-mesh](code/octonion-mesh/)); however, the functionality is fairly basic, and I don't anticipate any changes to the functionality. In other words, it shouldn't matter which one gets called.

## Getting started
Look at [interp5DOF.m](code/interp5DOF.m), which is a top-level function for creating a mesh, importing/generating data, triangulating a mesh, identifying the intersecting facet for datapoints, and finally, computing an interpolation. Also consider looking at [run.m](code/run.m) which is a script that contains more options, but may be less portable than the interp5DOF() function.

interp5DOF.m can be called in other functions/scripts to produce interpolation results using 5DOF misorientation/boundary plane normal pairs (qm/nA) and grain boundary property values. It was written with loosely similar input/output structure to the MATLAB built-in function [interpn()](https://www.mathworks.com/help/matlab/ref/interpn.html) involving input points, input values, query points, and query values.

## Test functions
Most functions have a corresponding "test" function (e.g. hsphext_subdiv.m --> hsphext_subdiv_test.m) which gives simple usage example(s). These are useful for debugging, visualizations, and understanding the functions without having to do a full run which could be time-consuming. This also allows for the non-test function code to be more succinct, as different plotting routines can be moved to the test function instead. The various test functions generally run to completion within a few seconds, and the parameters may be changed freely (e.g. dimension, number of points, etc.) where applicable. Some test functions have specific plotting routines for 1-sphere (2D) and 2-sphere (3D) cases since a 7-sphere is difficult to visualize and interpret ([n-sphere](https://en.wikipedia.org/wiki/N-sphere)).

## parfor loops
Parfor loops are used by default where there is potential for significant speed-up. A parfor-compatible text progress bar is encoded into many of these. Using disp() or fprintf() inside the parfor loop (aside from what's already inside the nested text progress bar function) may cause odd behavior on the command line, but should not affect the integrity of the code execution. Because these contain nested functions, in order to [deal with the inability to add variables to static workspaces](https://www.mathworks.com/help/matlab/matlab_prog/variables-in-nested-and-anonymous-functions.html), you can either assign variables to "ans" (a special variable that is still accessible), output statements directly to the command line terminal (no variable assignment), or comment the nested function, nUpdateProgress().

If the parallel computing toolbox is not installed, the parfor loops will execute as regular for loops. If the parallel computing toolbox is installed and only want to use a single core, start a parallel pool with only one core before running any of the functions via `parpool(1)`. The loop will still run as a parfor loop, however. A parfor loop with a single core and parallel computing toolbox should not run any slower than a regular for loop as long as they are contained within functions. A parfor loop executed within a script, however, is likely to result in significant slow-down.

To debug within a parfor loop, simply change it to a for loop while debugging. I added "parfor compatible" as a comment next to the parfor statements. Thus, you can use MATLAB [find files](https://www.mathworks.com/help/matlab/matlab_env/finding-files-and-folders.html#:~:text=To%20open%20the%20Find%20Files,on%20the%20MATLAB%20search%20path.) (Ctrl+Shift+F) to search for the keyword "parfor compatible" (including quotes) in order to keep track of which parfor loops have been changed to for loops.

## File dependencies
Take a look at [parseReqFiles_test.m](code/parseReqFiles_test.m) for generating a list of file dependencies for [interp5DOF.m](code/inter5DOF.m) (below) or other files.

1. [GB5DOF.m](code/GB5DOF.m)
1. [GB5DOF_setup.m](code/GB5DOF_setup.m)
1. [GBdist4.m](code/GBdist4.m)
1. [GBoct2five.m](code/GBoct2five.m)
1. [addpathdir.m](code/addpathdir.m)
1. [allcomb.m](code/allcomb.m)
1. [constructGBMatrices.m](code/constructGBMatrices.m)
1. [disorientation.m](code/disorientation.m)
1. [findgeometry.m](code/findgeometry.m)
1. [get_interp.m](code/get_interp.m)
1. [get_octpairs.m](code/get_octpairs.m)
1. [get_omega.m](code/get_omega.m)
1. [get_sympairs.m](code/get_sympairs.m)
1. [inmisFZ.m](code/inmisFZ.m)
1. [interp5DOF.m](code/interp5DOF.m)
1. [intersect_facet.m](code/intersect_facet.m)
1. [misFZcon.m](code/misFZcon.m)
1. [mustBeSqrt2Norm.m](code/mustBeSqrt2Norm.m)
1. [mustContainFields.m](code/mustContainFields.m)
1. [normr.m](code/normr.m)
1. [numStabBary.m](code/numStabBary.m)
1. [osymset.m](code/osymset.m)
1. [osymsets.m](code/osymsets.m)
1. [plotFZrodriguez.m](code/plotFZrodriguez.m)
1. [plotFZrodriguez_vtx.m](code/plotFZrodriguez_vtx.m)
1. [proj_down.m](code/proj_down.m)
1. [proj_up.m](code/proj_up.m)
1. [projfacet2hyperplane.m](code/projfacet2hyperplane.m)
1. [projray2hyperplane.m](code/projray2hyperplane.m)
1. [projray2hypersphere.m](code/projray2hypersphere.m)
1. [sphbary.m](code/sphbary.m)
1. [sphbary_setup.m](code/sphbary_setup.m)
1. [sphconvhulln.m](code/sphconvhulln.m)
1. [sqrt2norm.m](code/sqrt2norm.m)
1. [symaxis.m](code/symaxis.m)
1. [var_names.m](code/var_names.m)
1. [vert2con.m](code/vert2con.m)
1. [zeta_min2.m](code/zeta_min2.m)
