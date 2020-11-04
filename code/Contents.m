% CODE
%
% Files
%   addpathdir                        - add folders of filenames using addpath() and dir() 
%   allcomb                           - All combinations (Cartesian Product)
%   axpolytope                        - Indices of a polytope with vertices on every axis and every intersection of every plane formed by axes.
%   axpolytope_test                   - generate vertices on the axes of a polytope in n-D
%   chebycenter                       - Compute Chebyshev center of polytope Ax <= b.
%   constant_colorbar                 - script to generate a colorbar with constant limits
%   constructGBMatrices               - Make Olmsted GB matrices for GB5DOF using labframe quaternion/normal inputs
%   correctdis                        - convert to disorientation (for plotting)
%   cprnd                             - Draw from the uniform distribution over a convex polytope.
%   datagen                           - generate octonion data for various literature datasets (or
%   datagen_setup                     - setup for generating octonion data (random or from literature)
%   datagen_test                      - test generation of octonion data
%   degdelaunayn                      - compute delaunay triangulation of d-1 hyperplane in d-dimensions
%   disorientation                    - determines the unique misorientation between two adjacent
%   DisplayRequiredFunctions          - List required files and toolboxes.  Displays them in the command window or console window (if deployed).
%   dynamicCellExample                - Dynamically creating a nested cell structure (unfinished/deprecated)
%   facet_subdiv                      - Project a facet from n-dimensional space to a simplex in n-1
%   facet_subdiv_test                 - SUBDIV_TEST  e.g. subdivide a d-1 or d-2 simplex (e.g. triangle or line in 3D Cartesian coordinates) and plot
%   findgeometry                      - FINDGEOMETRY: Output the misorientation FZ geometry of a point given a quaternion (e.g. 'OAB').
%   findgeometry_test                 - test findgeometry for 'OAB', 'OBCE', 'OADE', etc. of misorientation FZ
%   fz_inserter                       - Compute misorientations in the fundamental zone for two lists of orientations
%   GB5DOF                            - computes the energy of an arbitrary boundary in FCC metals (BRK energy function)
%   GB5DOF_setup                      - Compute 5DOF GB energy from BRK function
%   GB5DOF_setup_test                 - GB5DOF_setup test
%   GBdist2                           - troubleshooting for GBdist() (deprecated)
%   GBdist3                           - troubleshooting for GBdist() (deprecated)
%   GBdist4                           - modified version of GBdist function by CMU group. Keeps o1 constant.
%   GBdist4_r2018a                    - arguments
%   GBdistdis                         - modified version of GBdist() (deprecated)
%   GBdistEucl                        - modified version of GBdist() for Euclidean distances (deprecated)
%   GBlab2oct                         - convert lab coordinate grain boundareis to octonions
%   GBlab2oct_test                    - Test from GB Octonion Tutorial script
%   GBoct2five                        - (incorrect, derivation error, 2020-11-03) Inverse operation of GBfive2oct.m (GB_octonion_code), output "five" struct
%   GBoct2five_r2018a                 - arguments
%   GBoct2five_test                   - 
%   GBoct2five_test2                  - 
%   GBpair                            - (deprecated) Method 1: Find o3 that has the minimum summed distances, o1-->o3, and o2-->o3.
%   GBpair2                           - deprecated
%   GBpair_r2018a                     - arguments
%   GBpair_test                       - test
%   gcfpos                            - get current figure position and copy to clipboard
%   get_cubo                          - get n quaternions from randomly or uniformly sampled cubochoric points
%   get_cubo_r2018a                   - arguments
%   get_errmetrics                    - get various error metrics for measured data relative to
%   get_failed                        - find failed jobs (under construction)
%   get_failed_test                   - test getfailed.m
%   get_fname                         - Get a filename to load or save (deprecated, specific to old naming scheme)
%   get_interp                        - Interpolate query point values based on spherical or planar barycentric coords in a mesh
%   get_octpairs                      - Get a set of octonions that are symmetrized with respect to a fixed reference GB (rng seed == 10)
%   get_octpairs2                     - deprecated version of get_octpairs()
%   get_octpairs_r2018a               - arguments
%   get_octpairs_test                 - get_octpairs test
%   get_ocubo                         - get octonions formed by pairs of quaternions from randomly or uniformly sampled cubochoric points.
%   get_ocubo_r2018a                  - arguments
%   get_ocubo_test                    - get_ocubo test
%   get_omega                         - calculate octonion distance
%   get_omega_r2018a                  - arguments
%   get_repsets                       - Find sets of points with at least a degeneracy of n
%   get_repsets_test                  - nonunique test
%   get_sympairs                      - get all combinations (pairs) of operators for a point group
%   get_sympairs_r2018a               - arguments
%   hsphext_subdiv                    - Hypersphere exterior hull subdivision.
%   hsphext_subdiv_test               - hypersphere exterior hull subdivision test
%   hypercube                         - Calculate vertices and convex hull triangulation of a hypercube
%   hypercube_test                    - hypercube test
%   hyperquadrant                     - Generate n-dimensional Ellipsoid or Sphere
%   hypersphere                       - Generate n-dimensional Ellipsoid or Sphere
%   hypersphere_subdiv                - Subdivide a "spherical" convex hull and collapse the triangulation.
%   hypersphere_subdiv_test           - 
%   hypersphereSetup                  - wrapper function for hypersphere.m, includes 'orthant' option
%   inBPFZ                            - return logical of which boundary plane normals are in BP fundamental zone
%   InCubicFZ                         - Test whether or not a set of rodrigues vectors is within the cubic misorientation fundamental zone.
%   inhull                            - tests if a set of points are inside a convex hull
%   inhull_setup                      - wrapper function for inhull (deprecated)
%   inmisFZ                           - check which rodrigues vectors fall inside the misorientation fundamental zone via vert2con (FEX)
%   inmisFZ_test                      - inmisFZ test
%   insertrows                        - Insert rows into a matrix at specific locations
%   interp5DOF                        - Convert misorientation and boundary plane normal 5DOF input
%   interp5DOF_setup                  - setup for interpolating five-degree-of-freedom property
%   interp5DOF_test                   - interp5DOF test
%   interpplot                        - Create parity plot with 5DOF plots on side (deprecated)
%   intersect_facet                   - Find intersection of ray with facet using barycentric coordinates.
%   intersect_facet_test              - intersect facet test
%   ismembc_test                      - ismembc vs. ismember test
%   knninterp                         - compute linear interpolation using hyperplane fitted to k-nearest neighbors
%   knninterp_test                    - knninterp test
%   mesh5DOF                          - generate "five" and "o" for a five degree-of-freedom fundamental zone (deprecated)
%   meshBP                            - generate boundary plane (BP) fundamental zone mesh for a
%   meshBP_test                       - test meshBP
%   meshFZ                            - mesh the misorientation fundamental zone (issues with tetgen)
%   meshgen                           - generate octonions randomly or from literature (deprecated)
%   meshgen_test                      - 
%   misFZcon                          - get constraints for misorientation fundamental zone (misFZ)
%   misFZcon_test                     - misFZcon test
%   mustBeSqrt2Norm                   - check that first octonion in list has norm == sqrt(2) and each quaternion has norm == 1
%   mustContainFields                 - check fieldnames(S) and make sure every checkname exists
%   myismember                        - do an ismembertol by rows with set precision and tolerance
%   n2c                               - t = num2cell(pts,1)
%   n2cplot                           - convert rows of points to cells of points, convenient for plotting.
%   normr                             - normalizes vectors row-by-row. Outputs zero vector if zero vector input (shadowed by Computer Vision Toolbox)
%   numStabBary                       - a numerically stable barycentric approach in high dimensions
%   numStabBary_example               - 
%   numSubplots                       - Calculate how many rows and columns of sub-plots are needed to neatly display n subplots
%   Oh_pg                             - Oh point group load/testing function
%   optimize_zeta                     - Minimize the pairwise distance matrix of a set of GBs w.r.t. zeta (twist angle via U(1) symmetry)
%   optimize_zeta_r2018a              - arguments
%   optimize_zeta_test                - optimize_zeta test
%   optimize_zeta_test_r2018a         - optimize_zeta test
%   orthoplex                         - Find indices and convex hull triangulation of an orthoplex
%   orthoplex_test                    - orthoplex test
%   OSLERP_setup                      - script to call OSLERP (deprecated)
%   OSLERP_setup2                     - another script to call OSLERP (deprecated)
%   osymset                           - get symmetrically equivalent octonions
%   osymset_r2018a                    - arguments
%   osymsets                          - Get symmetrically equivalent octonions (osymsets) for each octonion in a list of octonions
%   osymsets_r2018a                   - arguments
%   parityplot                        - Create a parity plot and pass options to (scatter() or hexscatter()) and refline().
%   parseReqFiles                     - parse required files & products using fileseparator sep
%   parseReqFiles_test                - parseMyFiles_test
%   pd_sse                            - get the error of a pairwise distance matrix relative to the true pairwise distance matrix
%   pd_sse_r2018a                     - arguments
%   plot5DOF                          - plotting five degree-of-freedom parameters in misorientation and boundary plane spaces
%   plotFZrodriguez                   - define the FZ vertices
%   plotFZrodriguez_test              - PLOTFZRODRIGUEZ_TEST
%   plotFZrodriguez_vtx               - plotFZrodriguez with 'A','B','C','D','E','O' vertices
%   proj_down                         - project down by removing null dimensions (i.e. a rotation and translation) via singular value decomposition
%   proj_down_r2018a                  - arguments
%   proj_down_test                    - 
%   proj_up                           - project up (restore null dimensions) using "USV" struct from proj_down.m
%   projfacet2hyperplane              - project facet vertices onto a hyperplane defined by nvec
%   projfacet2hyperplane_test         - project facet to hyperplane test
%   projray2hyperplane                - Project ray (pt) from unit hypersphere to tangent hyperplane at another point (nvec)
%   projray2hyperplane_test           - project ray to hyperplane test
%   projray2hypersphere               - project ray to hypersphere, compute barycentric coordinates, compute intersecting facet
%   projray2hypersphere_test          - project ray to hypersphere test
%   randOctParityData                 - submit sets of jobs to supercomputer or run locally for interp5DOF()
%   readNODE                          - Read vertex positions from .node file, .node files are used by Stellar and Triangle
%   rescale_test                      - rescale test (test of built-in function)
%   rotationmat3D                     - creates a rotation matrix based on axis and angle
%   RotMatrix                         - N-dimensional Rotation matrix
%   roundp                            - rounds values in x to n digits past the decimal sign.
%   run                               - run - Run MATLAB script
%   run2                              - run test
%   sbatch_setup                      - setup for SLURM sbatch submissions (deprecated)
%   simplex_subdiv                    - Compute the subdivision cartesian coordinates of a simplex
%   slightRot                         - rotate points slightly in every dimension
%   slightRot_test                    - slight rotation matrix test
%   sphbary                           - compute spherical barycentric coordinates of a facet or simplex [1]
%   sphbary_setup                     - sphbary coords and intersections for points
%   sphbary_test                      - spherical barycentric coordinates test
%   sphconvhulln                      - Compute the "convex hull" on a spherical surface
%   sphconvhulln_test                 - spherical convhulln test
%   sphere_stereograph                - compute the stereographic image of points on a sphere.
%   sphere_stereograph_inverse        - compute stereographic preimages of points.
%   sphere_stereograph_test           - sphere stereograph test
%   sphere_stereographic_inverse_test - sphere stereographic inverse test
%   sphtri_subdiv                     - Subdivision scheme that segments the triangle formed by vertices of a spherical triangle
%   sphtri_subdiv_test                - 
%   sqrt2norm                         - take a set of octonions and give each row norm == sqrt(2) if (norm == 1) || (norm == sqrt(2))
%   sqrt2norm_r2018a                  - arguments
%   symaxis                           - Return the symmetry axes (if multiple) and the geometry given a quaternion
%   symaxis_test                      - {
%   tblvertcat                        - vertically catenate tables with different variables, filling in dummy values as needed
%   tblvertcat_test                   - tblvertcat_test
%   testForGBFZSymmetryPlot           - 
%   toBPFZ                            - Rotate an arbitrary boundary plane normal for quaternion into standard boundary plane fundamental zone
%   tofiveFZ                          - Take 5DOF data and rotate it into misorientation and boundary plane fundamental zones.
%   tricollapse                       - collapse triangulation of points (i.e. a triangulation involving repeat points).
%   vecpair2rmat                      - Compute a (non-unique) rotation matrix to go from v1 to v2.
%   vert2con                          - convert a set of points to the set of inequality constraints
%   vert2lcon                         - An extension of Michael Kleder's vert2con function, handles degeneracy
%   write_video                       - write a set of images to a video named movname with some defaults
%   write_video_test                  - 
%   zeta_min2                         - ZETA_MIN  Alternative version of CMU group function zeta_min(), vectorized by Sterling Baird
%   zeta_min2_r2018a                  - arguments
%   GBoct2mat                         - convertion grain boundary octonions to orientation matrices (**unfinished**)
%   Kim2oct                           - load Kim2011 dataset, convert to norm-symmetrized octonions
%   avgrepeats                        - average values for duplicate points, remove all but one from degenerate set
%   avgrepeats_test                   - avgrepeats test
%   bugtest                           - checking a "bug" (turns out actually intended) with arguments..end syntax
%   el2po                             - elevation angles to polar angles
%   el2po_test                        - po2el_test
%   get_alen                          - get arclength of sphere (acos of dot product)
%   get_charlbl                       - get character labels, e.g. '(a)', '(b)', ... for figures.
%   get_five                          - generate random cubochoric misorientation and random boundary plane normal pairs
%   get_knn                           - k-nearest neighbor points, distances, mean, and std
%   get_walltimefn                    - get a function handle for computing interp5DOF SLURM walltimes (e.g. 'gpr')
%   hexscatter                        - A scatter-plot substitute - generate a density plot using hexagonal patches.
%   idw                               - inverse-distance weighting interpolation
%   idw_tovar                         - Original inverse distance weight function by Andrew Tovar
%   multiparity                       - create tiled parity plots using cell parity data
%   multixyplots                      - create tiled xy-plots from variables in a table
%   nnhist                            - nearest neighbor distance histogram
%   nnhist_test                       - nnhist test
%   padcat                            - concatenate vectors with different lengths by padding with NaN
%   paperfigure                       - call figure in centimeters and with appropriate size
%   papertext                         - do a text() command for the figure labels '(a)', '(b)', etc. for num == 1, 2, etc.
%   plotting                          - plotting script for interp5DOF paper
%   po2el                             - convert polar angles to elevation angles
%   qlab2qm                           - Convert lab/sample frame quaternions of grain A and grain B and
%   qmA2nA                            - QLAB2FIVE Convert lab/sample frame quaternions of grain A and grain B and
%   quivplot                          - plot three quivers in the x-hat, y-hat, and z-hat directions
%   savefigpng                        - save a figure and print the figure as 300 DPI .png
%   sphplot                           - plot a simple sphere for visualization purposes
%   tblfilt                           - filter a table based on a struct of parameters (deprecated, see built-in findgroups() and splitapply())
%   tblfilt_test                      - parfilter_test
%   test_voronoisphere                - Script to test voronoisphere (modified from original (Bruno Luong) by Sterling Baird)
%   toBPFZ_test                       - 
%   vcell_solidangle                  - Compute the solid angles of All voronoi cell
%   voronoisphere                     - Compute the voronoi's diagram of points on the spheres S^2
%   xyplot                            - make an errorbar xyplot using a table (tbl) and output (G) from findgroups()
%   xyplots                           - plot multiple datasets on the same axes using a "master" table and findgroups()






