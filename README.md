# Five Degree-of-Freedom Grain Boundary Interpolation
 Code related to meshing and interpolation of grain boundaries by representing 5DOF of grain boundaries as grain boundary octonions and mapping them into a Voronoi Fundamental Zone. ([https://github.com/sgbaird-5DOF/interp](https://github.com/sgbaird-5DOF/interp))
 
See
> 1. Baird, S. G.; Homer, E. R.; Fullwood, D. T.; Johnson, O. K. Five Degree-of-Freedom Property Interpolation of Arbitrary Grain Boundaries via Voronoi Fundamental Zone Octonion Framework. 2021. https://arxiv.org/abs/2104.06575 (Submitted to Computational Materials Science)
> 1. [GB_octonion_code](https://github.com/ichesser/GB_octonion_code)
> 1. Chesser, I., Francis, T., De Graef, M., & Holm, E. A. (2020). Learning the Grain Boundary Manifold: Tools for Visualizing and Fitting Grain Boundary Properties. Acta Materialia. https://doi.org/10.2139/ssrn.3460311
> 1. Francis, T., Chesser, I., Singh, S., Holm, E. A., & De Graef, M. (2019). A geodesic octonion metric for grain boundaries. Acta Materialia, 166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034

## Table of Contents
Use the <img src=https://user-images.githubusercontent.com/45469701/116359125-a68e5200-a7bb-11eb-8eba-dbea2acf653c.png width=30> button at the top-left of the README window

## Dependencies
### MATLAB Version
MATLAB R2019b or higher (mainly for the [arguments ... end syntax checking](https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html) at
the beginning of functions, which is used extensively throughout). For users of R2007a - R2019a, I suggest removing the arguments ... end syntax for any functions that use this and replacing it with corresponding [inputParser()](https://www.mathworks.com/help/matlab/ref/inputparser.html) and [varargin](https://www.mathworks.com/help/matlab/ref/varargin.html) code to deal with variable input arguments, default parameter values, and repeating arguments. Alternatively, you could remove the arguments ... end syntax lines for each function and update every place that the function is called so that all input arguments are specified. Open up an issue if you need more details on this. Other functions may need to be replaced if they aren't available in early MATLAB versions.

### Toolboxes
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html) (for Gaussian Process Regression: [fitrgp()](https://www.mathworks.com/help/stats/fitrgp.html), [fitrgp.predict()](https://www.mathworks.com/help/stats/compactregressiongp.predict.html))
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) (optional)
- [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html) (optional, for [numStabBary.m](code/numStabBary.m))

<!---, but for fitrgp() may need to change `hyperopts = struct('UseParallel',true,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` to `hyperopts = struct('UseParallel',false,'Optimizer','bayesopt','MaxObjectiveEvaluations',maxhyperobj);` in [interp5DOF.m](code/interp5DOF.m) under "method-specific interpolation" section 'gpr' case.) --->
### Files
By cloning the GitHub repository according to the instructions here, all file dependencies should be included. See [File dependencies](https://github.com/sgbaird/octonion-mesh/blob/master/README.md#file-dependencies) for a list of files that `interp5DOF.m` depends on.

## Installation

### Basic steps:
* Step 0: download the [code](https://github.com/sgbaird-5DOF/interp.git)
* Step 1: set [interp-5DOF/code/](code/) as working directory
* Step 2: add subfolders to path (`addpath(genpath('.'))`) and run [interp5DOF_test.m](code/interp5DOF_test.m)

### Platform-specific directions

#### Windows
##### Step 0: download the code
[Download GitHub Desktop](https://desktop.github.com/) and [Git Bash](https://git-scm.com/downloads). For Git Bash, the default installation options should be fine. I prefer to use [Atom](https://atom.io/) as the text editor which has some slick integrations with git. Login to GitHub Desktop and make a dummy repository via `Ctrl+N` so that you can open Git Bash via GitHub Desktop. Then clone and/or fork `https://github.com/sgbaird-5DOF/interp.git` by opening the Git Bash command line (i.e. Menubar --> Repository --> "Open in Git Bash" or 
```
Ctrl+`
```
) then run the following command in the directory where you want the `interp` directory to appear:
```bash
git -c submodule.interp5DOF-paper.update=none clone --recurse-submodules https://github.com/sgbaird-5DOF/interp.git
```

The `-c submodule.xxxx.update=none` flag indicates that a particular (private) submodule be ignored. The public submodules should be downloaded automatically. In order to update these submodules, add these to GitHub desktop as well ("add from existing", navigate within interp folder to the submodule, click on submodule folder, and "add").

Alternatively, you can try cloning directly in GitHub Desktop or via the "Open in GitHub Desktop" button under ![image](https://user-images.githubusercontent.com/45469701/116357284-907f9200-a7b9-11eb-81a3-3f55d27b8017.png), but it will likely throw an error, and it may not correctly clone the `MATslurm` submodule by the time it reaches that error. `MATslurm` is a bare minimum requirement for running `interp5DOF.m`. Attempt at your own risk.

##### Step 1: open MATLAB and navigate to navigate to [interp-5DOF/code/](code/)
Set [interp-5DOF/code/](code/) as working directory via `cd` or GUI

##### Step 2: Add subfolders to path and run [interp5DOF_test.m](code/interp5DOF_test.m)

\>\> `addpath(genpath('.'))`

\>\> `interp5DOF_test`

<!---
## Accessing functions via addpathdir()
dir() and addpath() commands are used to locate functions in subfolders of the current working directory via a custom function [addpathdir.m](code/addpathdir.m). This could give anomalous behavior if the directory structure is changed such that filenames are non-unique in sub-folders of the parent folder where addpathdir() gets called, or if files with the same name are present elsewhere on the user's MATLAB path.
--->

#### Linux
##### Step 0: download the code
```bash
git -c submodule.interp5DOF-paper.update=none clone --recurse-submodules https://github.com/sgbaird-5DOF/interp.git
```

Verify that [MATslurm](https://github.com/sgbaird-5DOF/MATslurm) is not an empty directory.

##### Step 1: open MATLAB and navigate to navigate to [interp-5DOF/code/](code/)
`matlab`

\>\> `cd interp-5DOF/code/`

##### Step 2: Add subfolders to path and run [interp5DOF_test.m](code/interp5DOF_test.m)

\>\> `addpath(genpath('.'))`

\>\> `interp5DOF_test`

#### Troubleshooting
See [cloning a repository](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) and [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for more information or other options such as using GitHub Desktop (Windows, Linux, etc.) or downloading a .zip file. [Forking](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo), pull requests, and opening of issues are welcome/encouraged. Note that the .zip files will not contain submodules, which means that you'll need to download any submodules individually if you go that route (bare minimum would be the [MATslurm](https://github.com/sgbaird-5DOF/MATslurm) repository). Instead of downloading a .zip file, I suggest instead downloading and using GitHub Desktop. If you run into issues with the repository, you can download a [lightweight version of interp5DOF.m and its dependencies](code/interp5DOF-lightweight.zip), unzip it, and read README.txt for instructions, but please open up an issue in the GitHub repo if you do run into trouble with the repo. That way, you can get the most recent updates and others can benefit from it.

## Getting started
Look at [interp5DOF.m](code/interp5DOF.m), which is a high-level function for Gaussian Process Regression (GPR), barycentric, nearest neighbor (NN), and inverse-distance weighting (IDW) interpolation. This involves importing/generating data and computing an interpolation.

`interp5DOF.m` can be called in other functions/scripts to produce interpolation results using 5DOF misorientation/boundary plane normal pairs (qm/nA) and grain boundary property values. It was written with loosely similar input/output structure to the MATLAB built-in function [interpn()](https://www.mathworks.com/help/matlab/ref/interpn.html) involving input points and values, query points and values, and options.

### Simple Example Data
Separate from [interp5DOF_test.m](code/interp5DOF_test.m), the following is a fast, bare-bones example to show the basic input/output format of interp5DOF.m. See also [get_cubo.m](code/get_cubo.m)
```matlab
npts = 100;
qm = get_cubo(npts); nA = normr(rand(npts,3)); %random (qm,nA) pairs
propList = 1:npts; %property values
qm2 = get_cubo(npts); nA2 = normr(rand(npts,3)); %random (qm,nA) pairs
method = 'gpr'; %interpolation method
[propOut,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,propList,qm2,nA2,method)
```

## Test functions
Most functions have a corresponding "test" function (e.g. `hsphext_subdiv.m` --> `hsphext_subdiv_test.m`, `interp5DOF.m` --> `interp5DOF_test.m`) which gives simple usage example(s). These are useful for debugging, visualizations, and understanding the functions without having to do a full run which could be time-consuming. This also allows for the non-test function code to be more succinct, as certain plotting routines can be moved to the test function instead. The various test functions generally run to completion within a few seconds, and the parameters can generally be changed freely (e.g. dimension, number of points). Some test functions have specific plotting routines for 1-sphere (2D) and 2-sphere (3D) cases since a 7-sphere is difficult to visualize and interpret ([n-sphere](https://en.wikipedia.org/wiki/N-sphere)). For example, see [sphbary_test.m](code/sphbary_test.m) and [toBPFZ_test.m](code/toBPFZ_test.m).

## Distance Calculations
If you only want to (manually) compute distances in the VFZ sense, first you need to map all GBs into a VFZ.
```matlab
npts = 100;
o = get_ocubo(npts); %generate some random data
o = get_octpairs(o); %symmetrize (using default reference GBO)
```
At this point, you can get the VFZ pairwise distance matrix via:
```matlab
pd = pdist(o);
```
Alternatively, `pdist2()` may be of interest if you want pairwise distances between two sets of points, or `vecnorm()` if you want to calculate distances between two lists of GBOs:
```matlab
npts2 = 100;
o2 = get_ocubo(npts2);
o2 = get_octpairs(o2);
d = vecnorm(o1-o2,2,2); %o1 and o2 need to be the same size
```
If you want the "true" minimum distances (i.e. essentially the same implementation as [GB_octonion_code](https://github.com/ichesser/GB_octonion_code), but vectorized and parallelized), you may use `GBdist4.m` directly with two sets of GBOs.
```matlab
d = GBdist4(o1,o2);
```
It depends on the application, but if you want to compute large pairwise distance matrices that are nearly identical to the traditional GBO distances, I recommend using the ensembled VFZ distance via `ensembleGBdist.m` with `K >= 10`. This will be much faster than using `GBdist4.m`. For reference, this corresponds to (from the main paper when `K==10`):  
<img src=https://user-images.githubusercontent.com/45469701/116044929-bcbad780-a62e-11eb-8c59-58a4354badbb.png width=300>  
This is distinct from `ensembleVFZO.m`, which takes the average interpolated property from `K` different VFZs.

Drop me a note in "Issues" if you have something you'd like to do or something you'd like to clarify, but can't figure out among the (many) options and functions in the `interp` repo. With a few details, there's a good chance I can offer some suggestions that will save a lot of time.

## Plots
Most plots in the paper are produced in the script: [plotting.m](code/plotting.m). First, you need to create a dummy folder:
```bash
mkdir interp/code/interp5DOF-paper/figures/
```

The larger file dependencies, `gitID-0055bee_uuID-475a2dfd_paper-data6.mat` and `gpr46883_gitID-b473165_puuID-50ffdcf6_kim-rng11.mat` can be [downloaded at figshare](https://doi.org/10.6084/m9.figshare.14405924.v3) and have the following citation:
> @misc{baird_homer_fullwood_johnson_2021, title={Five Degree-of-Freedom Grain Boundary Interpolation}, url={https://figshare.com/articles/dataset/gitID-0055bee_uuID-475a2dfd_paper-data6_mat/14405924/3}, DOI={10.6084/m9.figshare.14405924.v3}, abstractNote={These are larger MATLAB .mat files required for reproducing plots from the sgbaird-5DOF/interp repository for grain boundary property interpolation. gitID-0055bee_uuID-475a2dfd_paper-data6.mat contains multiple trials of five degree-of-freedom interpolation model runs for various interpolation schemes. gpr46883_gitID-b473165_puuID-50ffdcf6_kim-rng11.mat contains a Gaussian Process Regression model trained on 46883 Fe simulation GBs.}, publisher={figshare}, author={Baird, Sterling and Homer, Eric R and Fullwood, David and Johnson, Oliver K.}, year={2021}, month={Apr} }

## Data Preparation
Input GBs can take on the following forms:
- misorientation (`qm`) / boundary plane normal (`nA`) pairs
- octonions (`o`)

Prediction GBs have a `2` appended, as in `qm2`. The inputs are stacked vertically (i.e. 1st row corresponds to 1st GB, 2nd row to 2nd GB, etc.).

Data can be converted between forms/conventions using the various rotation functions available in the repository (modified versions of [GB Octonion Code](https://github.com/ichesser/GB_octonion_code/tree/master/TutorialCode/rotation_conversions)).

### Active vs. Passive Rotation Convention
Active rotation is specified with an input parameter `'epsijk'` == `1` and is the default used throughout the work. Passive rotation, in theory, can be specified via `epsijk` == `-1` which should be propagated throughout the entire codebase via the top-level input, but this has not undergone the same extensive testing as the active rotation convention. To convert from the passive to the active rotation convention, I suggest converting the GB into octonion form and applying `qinv.m` to each grain boundary, that is:
```matlab
o = [qinv(o(1:4,:) qinv(o(5:8,:))];
```
or simply
```matlab
o = oflip(o);
```
## parfor loops
Parfor loops are used by default where there is potential for significant speed-up. A parfor-compatible text progress bar is encoded into many of these. Adding disp() or fprintf() inside the parfor loop (aside from what's already inside the nested text progress bar function) may cause odd behavior on the command line output, but should not affect the integrity of the code execution. Because the parfor-compatible text progress bars need to be nested functions, in order to [deal with the inability to add variables to static workspaces](https://www.mathworks.com/help/matlab/matlab_prog/variables-in-nested-and-anonymous-functions.html) while debugging, you can either assign variables to "ans" (a special variable that is still accessible), output statements directly to the command line terminal (no variable assignment). Alternatively, you can comment the nested function, `nUpdateProgress()`.

If the parallel computing toolbox is not installed, the `parfor` loops will execute as regular `for` loops. If the parallel computing toolbox is installed and you only want to use a single core, start a parallel pool with only one core before running any of the functions via `parpool(1)`. The loop will still run internally as a `parfor` loop, however. A `parfor` loop with a single core and parallel computing toolbox should not run any slower than a regular `for` loop as long as they are contained within functions. A `parfor` loop executed within a script, however, is likely to result in significant slow-down.

To debug within a `parfor` loop, simply change it to a `for` loop while debugging and change it back afterwards. I added "parfor compatible" as a comment next to the parfor statements. Thus, you can use MATLAB [find files](https://www.mathworks.com/help/matlab/matlab_env/finding-files-and-folders.html#:~:text=To%20open%20the%20Find%20Files,on%20the%20MATLAB%20search%20path.) (`Ctrl+Shift+F`) to search for the keyword "parfor compatible" (including quotes) in order to keep track of which parfor loops have been changed to for loops. If you make changes and an error arises, it is possible it will only give useful information at the top-level where the `parfor` started, hence the debugging suggestion above.

## File dependencies
Take a look at [parseReqFiles_test.m](code/parseReqFiles_test.m) for generating a list of file dependencies for [interp5DOF.m](code/inter5DOF.m) (below) or other files.

1. [GB5DOF.m](code/GB5DOF.m)
1. [GB5DOF_setup.m](code/GB5DOF_setup.m)
1. [PGnames.mat](code/GB_octonion_code-master_CMU/TutorialCode/crystal_symmetry_ops/PGnames.mat)
1. [PGsymops.mat](code/GB_octonion_code-master_CMU/TutorialCode/crystal_symmetry_ops/PGsymops.mat)
1. [GBfive2oct.m](code/GB_octonion_code-master_CMU/TutorialCode/octonion_functions/GBfive2oct.m)
1. [qinv_francis.m](code/GB_octonion_code-master_CMU/TutorialCode/octonion_functions/qinv_francis.m)
1. [qmult.m](code/GB_octonion_code-master_CMU/TutorialCode/octonion_functions/qmult.m)
1. [GetPyramid.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/GetPyramid.m)
1. [ax2qu.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/ax2qu.m)
1. [cu2ho.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/cu2ho.m)
1. [cu2qu.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/cu2qu.m)
1. [ho2ax.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/ho2ax.m)
1. [ho2qu.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/ho2qu.m)
1. [qu2om.m](code/GB_octonion_code-master_CMU/TutorialCode/rotation_conversions/qu2om.m)
1. [GBdist4.m](code/GBdist4.m)
1. [GBlab2oct.m](code/GBlab2oct.m)
1. [get_gitcommit.m](code/MATslurm/code/get_gitcommit.m)
1. [get_uuid.m](code/MATslurm/code/get_uuid.m)
1. [structhorzcat.m](code/MATslurm/code/structhorzcat.m)
1. [var_names.m](code/MATslurm/code/var_names.m)
1. [addpathdir.m](code/addpathdir.m)
1. [allcomb.m](code/allcomb.m)
1. [constructGBMatrices.m](code/constructGBMatrices.m)
1. [get_cubo.m](code/get_cubo.m)
1. [get_errmetrics.m](code/get_errmetrics.m)
1. [get_interp.m](code/get_interp.m)
1. [get_knn.m](code/get_knn.m)
1. [get_octpairs.m](code/get_octpairs.m)
1. [get_ocubo.m](code/get_ocubo.m)
1. [get_omega.m](code/get_omega.m)
1. [get_sympairs.m](code/get_sympairs.m)
1. [idw.m](code/idw.m)
1. [get_ppts.m](code/interpfns/get_ppts.m)
1. [get_pts.m](code/interpfns/get_pts.m)
1. [interp_avg.m](code/interpfns/interp_avg.m)
1. [interp_bary.m](code/interpfns/interp_bary.m)
1. [interp_bary_fast.m](code/interpfns/interp_bary_fast.m)
1. [interp_gpr.m](code/interpfns/interp_gpr.m)
1. [interp_idw.m](code/interpfns/interp_idw.m)
1. [interp_nn.m](code/interpfns/interp_nn.m)
1. [intersect_facet.m](code/intersect_facet.m)
1. [mustBeSqrt2Norm.m](code/mustBeSqrt2Norm.m)
1. [mustContainFields.m](code/mustContainFields.m)
1. [normr.m](code/normr.m)
1. [numStabBary.m](code/numStabBary.m)
1. [osymset.m](code/osymset.m)
1. [osymsets.m](code/osymsets.m)
1. [proj_down.m](code/proj_down.m)
1. [proj_up.m](code/proj_up.m)
1. [projfacet2hyperplane.m](code/projfacet2hyperplane.m)
1. [projray2hyperplane.m](code/projray2hyperplane.m)
1. [projray2hypersphere.m](code/projray2hypersphere.m)
1. [gmat2q.m](code/quaternions/gmat2q.m)
1. [q2gmat.m](code/quaternions/q2gmat.m)
1. [qconj.m](code/quaternions/qconj.m)
1. [qinv_johnson.m](code/quaternions/qinv_johnson.m)
1. [qmultiply.m](code/quaternions/qmultiply.m)
1. [qnorm.m](code/quaternions/qnorm.m)
1. [sphbary.m](code/sphbary.m)
1. [sphbary_setup.m](code/sphbary_setup.m)
1. [sphconvhulln.m](code/sphconvhulln.m)
1. [sqrt2norm.m](code/sqrt2norm.m)
1. [zeta_min2.m](code/zeta_min2.m)
