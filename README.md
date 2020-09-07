# octonion-mesh
 code related to meshing and interpolation of grain boundary octonions

## Dependencies
### MATLAB Version
MATLAB R2019b or higher (mainly for the "arguments" syntax checking at
the beginning of functions, which is used extensively throughout). For users of R2007a - R2019a, I suggest removing the arguments ... end syntax for any functions that use this and replacing it with corresponding [inputParser()](https://www.mathworks.com/help/matlab/ref/inputparser.html) and [varargin](https://www.mathworks.com/help/matlab/ref/varargin.html) code to deal with variable input arguments, default parameter values, and repeating arguments. Alternatively, you could remove the arguments ... end syntax lines for each function and update every place that the function is called so that all input arguments are specified. Open up an issue if you need more details on this.

### Toolboxes
- Statistics and Machine Learning Toolbox (for Gaussian Process Regression: [fitrgp()](https://www.mathworks.com/help/stats/fitrgp.html), [fitrgp.predict()](https://www.mathworks.com/help/stats/compactregressiongp.predict.html))
- Parallel Computing Toolbox (optional, but for fitrgp() may need to change

`hyperopts = struct('UseParallel',true,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` to `hyperopts = struct('UseParallel',false,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` in [interp5DOF.m](code/interp5DOF.m) under "method-specific interpolation" section --> 'gpr' case.)
- Symbolic Math Toolbox (optional, for numStabBary.m)

## Usage
### Step0: download the code
`git clone --recurse-submodules https://github.com/sgbaird/octonion-inference.git`
See [cloning a repository](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) and [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for more information or other options.

### Step1: navigate to [octonion-inference/code/](code/)
`cd octonion-inference/code/`

### Step3: open MATLAB and call [interp5DOF_test.m](code/interp5DOF_test.m) or [run.m](code/run.m)
`matlab`

\>\> `interp5DOF_test`

## Accessing functions via addpathdir()
dir() and addpath() commands are used to locate functions in subfolders of the current working directory via a custom function [code/addpathdir.m](https://github.com/sgbaird/octonion-inference/blob/master/code/addpathdir.m). This could give anomalous behavior if the directory structure is changed such that filenames are non-unique in sub-folders of the parent folder where addpathdir() gets called, or if files with the same name are present elsewhere on the user's MATLAB path. This is also the only function that is shadowed (to my knowledge) within this repository ([by octonion-mesh](code/octonion-mesh/)); however, the functionality is fairly basic, and I don't anticipate any changes to the functionality. In other words, it shouldn't matter which one gets called.

## Getting started
Look at [code/interp5DOF.m](code/interp5DOF.m), which is a top-level function for creating a mesh, importing/generating data, triangulating a mesh, identifying the intersecting facet for datapoints, and finally, computing an interpolation. Also consider looking at [code/run.m](code/run.m) which is a script that contains more options, but may be less portable than the interp5DOF() function. interp5DOF.m can be called in other functions/scripts to produce interpolation results using 5DOF misorientation/boundary plane normal pairs and grain boundary property values. It was written with loosely similar input/output structure to the MATLAB built-in function interpn() involving input points, input values, query points, and output values.

## Test functions
Most functions have a corresponding "test" function (e.g. hsphext_subdiv.m --> hsphext_subdiv_test.m) which gives simple usage example(s). These are useful for debugging, visualizations, and understanding the functions without having to do a full run which could be time-consuming. This also allows for the non-test function code to be more succinct, as different plotting routines can be moved to the test function instead. The various test functions generally run to completion within a few seconds, and the parameters may be changed freely (e.g. dimension, number of points, etc.) where applicable. Some test functions have specific plotting routines for 1-sphere (2D) and 2-sphere (3D) cases since a 7-sphere is difficult to visualize and interpret ([n-sphere](https://en.wikipedia.org/wiki/N-sphere)).

## parfor loops
Parfor loops are used by default where there is potential for significant speed-up. A parfor-compatible text progress bar is encoded into many of these. Using disp() or fprintf() inside the parfor loop (aside from what's already inside the nested text progress bar function) may cause odd behavior on the command line, but should not affect the integrity of the code execution. Because these contain nested functions, in order to add new variables to the "static workspace", you can either assign variables to "ans" (a special variable that is still accessible), output statements directly to the command line terminal, or comment the nested function, nUpdateProgress(). If the parallel computing toolbox is not installed, the parfor loops will execute as regular for loops. If the parallel computing toolbox is installed, and you want to use only a single core, start a parallel pool with only one core before running any of the functions. The loop will still do all of its internal parfor executions. A parfor loop with a single core and parallel computing toolbox should not run any slower than a regular for loop when they are contained within functions. A parfor loop executed within a script, however, is likely to result in significant slow-down. To debug within a parfor loop, simply change it to a for loop while debugging. MATLAB [find files](https://www.mathworks.com/help/matlab/matlab_env/finding-files-and-folders.html#:~:text=To%20open%20the%20Find%20Files,on%20the%20MATLAB%20search%20path.) (Ctrl+Shift+F) searching for the keywords "parfor compatible" (including quotes) should help in keeping track of which ones have been changed. If you do not have the Parallel Computing Toolbox

## File dependencies
Take a look at [parseMyFiles_test.m](code/parseMyFiles_test.m) for generating a list of function dependencies for [interp5DOF.m](code/inter5DOF.m) or other files.

1. [addpathdir.m](octonion-inference\codeaddpathdir.m)
1. [interp5DOF.m](octonion-inference\codeinterp5DOF.m)
1. [GB5DOF.m](octonion-inference\code\octonion-mesh\codeGB5DOF.m)
1. [GB5DOF_setup.m](octonion-inference\code\octonion-mesh\codeGB5DOF_setup.m)
1. [PGnames.mat](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\crystal_symmetry_opsPGnames.mat)
1. [PGsymops.mat](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\crystal_symmetry_opsPGsymops.mat)
1. [GBfive2oct.m](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\octonion_functionsGBfive2oct.m)
1. [qinv.m](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\octonion_functionsqinv.m)
1. [qmult.m](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\octonion_functionsqmult.m)
1. [ax2qu.m](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\rotation_conversionsax2qu.m)
1. [qu2ax.m](octonion-inference\code\octonion-mesh\code\GB_octonion_code-master_CMU\TutorialCode\rotation_conversionsqu2ax.m)
1. [GBdist4.m](octonion-inference\code\octonion-mesh\codeGBdist4.m)
1. [GBoct2five.m](octonion-inference\code\octonion-mesh\codeGBoct2five.m)
1. [allcomb.m](octonion-inference\code\octonion-mesh\codeallcomb.m)
1. [constructGBMatrices.m](octonion-inference\code\octonion-mesh\codeconstructGBMatrices.m)
1. [misFZfeatures.mat](octonion-inference\code\octonion-mesh\code\datamisFZfeatures.mat)
1. [disorientation.m](octonion-inference\code\octonion-mesh\codedisorientation.m)
1. [findgeometry.m](octonion-inference\code\octonion-mesh\codefindgeometry.m)
1. [get_interp.m](octonion-inference\code\octonion-mesh\codeget_interp.m)
1. [get_octpairs.m](octonion-inference\code\octonion-mesh\codeget_octpairs.m)
1. [get_omega.m](octonion-inference\code\octonion-mesh\codeget_omega.m)
1. [get_sympairs.m](octonion-inference\code\octonion-mesh\codeget_sympairs.m)
1. [inmisFZ.m](octonion-inference\code\octonion-mesh\codeinmisFZ.m)
1. [intersect_facet.m](octonion-inference\code\octonion-mesh\codeintersect_facet.m)
1. [misFZcon.m](octonion-inference\code\octonion-mesh\codemisFZcon.m)
1. [mustBeSqrt2Norm.m](octonion-inference\code\octonion-mesh\codemustBeSqrt2Norm.m)
1. [mustContainFields.m](octonion-inference\code\octonion-mesh\codemustContainFields.m)
1. [normr.m](octonion-inference\code\octonion-mesh\codenormr.m)
1. [numStabBary.m](octonion-inference\code\octonion-mesh\codenumStabBary.m)
1. [osymset.m](octonion-inference\code\octonion-mesh\codeosymset.m)
1. [osymsets.m](octonion-inference\code\octonion-mesh\codeosymsets.m)
1. [plotFZrodriguez.m](octonion-inference\code\octonion-mesh\codeplotFZrodriguez.m)
1. [plotFZrodriguez_vtx.m](octonion-inference\code\octonion-mesh\codeplotFZrodriguez_vtx.m)
1. [proj_down.m](octonion-inference\code\octonion-mesh\codeproj_down.m)
1. [proj_up.m](octonion-inference\code\octonion-mesh\codeproj_up.m)
1. [projfacet2hyperplane.m](octonion-inference\code\octonion-mesh\codeprojfacet2hyperplane.m)
1. [projray2hyperplane.m](octonion-inference\code\octonion-mesh\codeprojray2hyperplane.m)
1. [projray2hypersphere.m](octonion-inference\code\octonion-mesh\codeprojray2hypersphere.m)
1. [gmat2q.m](octonion-inference\code\octonion-mesh\code\quaternionsgmat2q.m)
1. [q2gmat.m](octonion-inference\code\octonion-mesh\code\quaternionsq2gmat.m)
1. [q2rod.m](octonion-inference\code\octonion-mesh\code\quaternionsq2rod.m)
1. [qconj.m](octonion-inference\code\octonion-mesh\code\quaternionsqconj.m)
1. [qinv_johnson.m](octonion-inference\code\octonion-mesh\code\quaternionsqinv_johnson.m)
1. [qmultiply.m](octonion-inference\code\octonion-mesh\code\quaternionsqmultiply.m)
1. [qnorm.m](octonion-inference\code\octonion-mesh\code\quaternionsqnorm.m)
1. [sphbary.m](octonion-inference\code\octonion-mesh\codesphbary.m)
1. [sphbary_setup.m](octonion-inference\code\octonion-mesh\codesphbary_setup.m)
1. [sphconvhulln.m](octonion-inference\code\octonion-mesh\codesphconvhulln.m)
1. [sqrt2norm.m](octonion-inference\code\octonion-mesh\codesqrt2norm.m)
1. [symaxis.m](octonion-inference\code\octonion-mesh\codesymaxis.m)
1. [var_names.m](octonion-inference\code\octonion-mesh\codevar_names.m)
1. [vert2con.m](octonion-inference\code\octonion-mesh\codevert2con.m)
1. [zeta_min2.m](octonion-inference\code\octonion-mesh\codezeta_min2.m)
