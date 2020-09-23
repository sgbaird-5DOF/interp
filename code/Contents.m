% CODE
%
% Files
%   addpathdir                        - add folders of filenames to path by searching through subfolders and
%   allcomb                           - - All combinations
%   axpolytope                        - --------------------------------------------------------------------------
%   axpolytope_test                   - axpolytope test
%   chebycenter                       - Compute Chebyshev center of polytope Ax <= b.
%   constant_colorbar                 - ----colorbar------
%   constructGBMatrices               - -------------------------------------------------------------------------%
%   correctdis                        - convert to disorientation (for plotting)
%   cprnd                             - Draw from the uniform distribution over a convex polytope.
%   datagen                           - --------------------------------------------------------------------------
%   datagen_setup                     - --------------------------------------------------------------------------
%   datagen_test                      - datagen test
%   degdelaunayn                      - {
%   disorientation                    - -------------------------------------------------------------------------%
%   DisplayRequiredFunctions          - List required files and toolboxes.  Displays them in the command window or console window (if deployed).
%   dynamicCellExample                - 
%   facet_subdiv                      - --------------------------------------------------------------------------
%   facet_subdiv_test                 - 
%   findgeometry                      - -------------------------------------------------------------------------
%   findgeometry_test                 - 
%   fz_inserter                       - Compute misorientations in the fundamental zone for two lists of
%   GB5DOF                            - energy    computes the energy of an arbitrary boundary in FCC metals. 
%   GB5DOF_setup                      - --------------------------------------------------------------------------
%   GB5DOF_setup_test                 - GB5DOF_setup test
%   GBdist2                           - INPUT DATA
%   GBdist3                           - INPUT DATA 
%   GBdist4                           - --------------------------------------------------------------------------
%   GBdist4_r2018a                    - arguments
%   GBdistdis                         - INPUT DATA 
%   GBdistEucl                        - INPUT DATA
%   GBlab2oct                         - -------------------------------------------------------------------------%
%   GBlab2oct_test                    - Test from GB Octonion Tutorial script
%   GBoct2five                        - --------------------------------------------------------------------------
%   GBoct2five_r2018a                 - arguments
%   GBoct2five_test                   - 
%   GBoct2five_test2                  - 
%   GBpair                            - --------------------------------------------------------------------------
%   GBpair2                           - 
%   GBpair_r2018a                     - arguments
%   GBpair_test                       - GBpair test
%   gcfpos                            - gcfpos
%   get_cubo                          - GET_OCUBO  get n quaternions from randomly or uniformly sampled cubochoric points
%   get_cubo_r2018a                   - arguments
%   get_errmetrics                    - get various error metrics for measured data relative to
%   get_failed                        - find failed jobs
%   get_failed_test                   - 
%   get_fname                         - --------------------------------------------------------------------------
%   get_interp                        - --------------------------------------------------------------------------
%   get_octpairs                      - Author: Sterling Baird
%   get_octpairs2                     - --------------------------------------------------------------------------
%   get_octpairs_r2018a               - arguments
%   get_octpairs_test                 - get_octpairs test
%   get_ocubo                         - --------------------------------------------------------------------------
%   get_ocubo_r2018a                  - arguments
%   get_ocubo_test                    - get_ocubo test
%   get_omega                         - --------------------------------------------------------------------------
%   get_omega_r2018a                  - arguments
%   get_repsets                       - --------------------------------------------------------------------------
%   get_repsets_test                  - nonunique test
%   get_sympairs                      - --------------------------------------------------------------------------
%   get_sympairs_r2018a               - arguments
%   hsphext_subdiv                    - --------------------------------------------------------------------------
%   hsphext_subdiv_test               - hypersphere exterior hull subdivision test
%   hypercube                         - --------------------------------------------------------------------------
%   hypercube_test                    - hypercube test
%   hyperquadrant                     - HYPERSPHERE - Generate n-dimensional Ellipsoid or Sphere
%   hypersphere                       - - Generate n-dimensional Ellipsoid or Sphere
%   hypersphere_subdiv                - --------------------------------------------------------------------------
%   hypersphere_subdiv_test           - 
%   hypersphereSetup                  - --------------------------------------------------------------------------
%   inBPFZ                            - --------------------------------------------------------------------------
%   InCubicFZ                         - --------------------------------------------------------------------------
%   inhull                            - inhull: tests if a set of points are inside a convex hull
%   inhull_setup                      - --------------------------------------------------------------------------
%   inmisFZ                           - --------------------------------------------------------------------------
%   inmisFZ_test                      - inmisFZ test
%   insertrows                        - - Insert rows into a matrix at specific locations
%   interp5DOF                        - Convert misorientation and boundary plane normal 5DOF input
%   interp5DOF_setup                  - setup for interpolating five-degree-of-freedom property
%   interp5DOF_test                   - interp5DOF test
%   interpplot                        - --------------------------------------------------------------------------
%   intersect_facet                   - --------------------------------------------------------------------------
%   intersect_facet_test              - intersect facet test
%   ismembc_test                      - ismembc vs. ismember test
%   knninterp                         - --------------------------------------------------------------------------
%   knninterp_test                    - knninterp test
%   mesh5DOF                          - --------------------------------------------------------------------------
%   meshBP                            - --------------------------------------------------------------------------
%   meshBP_test                       - test meshBP
%   meshFZ                            - --------------------------------------------------------------------------
%   meshgen                           - 
%   meshgen_test                      - 
%   misFZcon                          - --------------------------------------------------------------------------
%   misFZcon_test                     - misFZcon test
%   mustBeSqrt2Norm                   - NOTE: do not use "arguments" syntax here since this is a validation fn
%   mustContainFields                 - 
%   myismember                        - 
%   n2c                               - 
%   n2cplot                           - convert rows of points to cells of points, convenient for plotting.
%   normr                             - --------------------------------------------------------------------------
%   numStabBary                       - {
%   numStabBary_example               - 
%   numSubplots                       - function [p,n]=numSubplots(n)
%   Oh_pg                             - 
%   optimize_zeta                     - --------------------------------------------------------------------------
%   optimize_zeta_r2018a              - arguments
%   optimize_zeta_test                - optimize_zeta test
%   optimize_zeta_test_r2018a         - optimize_zeta test
%   orthoplex                         - --------------------------------------------------------------------------
%   orthoplex_test                    - orthoplex test
%   OSLERP_setup                      - OSLERP_setup
%   OSLERP_setup2                     - OSLERP_setup2
%   osymset                           - --------------------------------------------------------------------------
%   osymset_r2018a                    - arguments
%   osymsets                          - --------------------------------------------------------------------------
%   osymsets_r2018a                   - arguments
%   parityplot                        - --------------------------------------------------------------------------
%   parseReqFiles                     - parse required files & products
%   parseReqFiles_test                - parseMyFiles_test
%   pd_sse                            - --------------------------------------------------------------------------
%   pd_sse_r2018a                     - arguments
%   plot5DOF                          - plotting
%   plotFZrodriguez                   - define the FZ vertices
%   plotFZrodriguez_test              - 
%   plotFZrodriguez_vtx               - 
%   proj_down                         - --------------------------------------------------------------------------
%   proj_down_r2018a                  - arguments
%   proj_down_test                    - proj_down test
%   proj_up                           - --------------------------------------------------------------------------
%   projfacet2hyperplane              - --------------------------------------------------------------------------
%   projfacet2hyperplane_test         - project facet to hyperplane test
%   projray2hyperplane                - --------------------------------------------------------------------------
%   projray2hyperplane_test           - project ray to hyperplane test
%   projray2hypersphere               - --------------------------------------------------------------------------
%   projray2hypersphere_test          - project ray to hypersphere test
%   randOctParityData                 - 
%   readNODE                          - Read vertex positions from .node file, .node files are used by
%   rescale_test                      - rescale test (test of built-in function)
%   rotationmat3D                     - function R= rotationmat3D(radians,Axis)
%   RotMatrix                         - RotMatrix - N-dimensional Rotation matrix
%   roundp                            - -------------------------------------------------------------------------%
%   run                               - run - Run MATLAB script
%   run2                              - run test
%   sbatch_setup                      - function job_setup
%   simplex_subdiv                    - ---------------------------------------------------------------------------
%   slightRot                         - --------------------------------------------------------------------------
%   slightRot_test                    - slight rotation matrix test
%   sphbary                           - --------------------------------------------------------------------------
%   sphbary_setup                     - --------------------------------------------------------------------------
%   sphbary_test                      - spherical barycentric coordinates test
%   sphconvhulln                      - --------------------------------------------------------------------------
%   sphconvhulln_test                 - spherical convhulln test
%   sphere_stereograph                - *****************************************************************************80
%   sphere_stereograph_inverse        - *****************************************************************************80
%   sphere_stereograph_test           - sphere stereograph test
%   sphere_stereographic_inverse_test - sphere stereographic inverse test
%   sphtri_subdiv                     - --------------------------------------------------------------------------
%   sphtri_subdiv_test                - 
%   sqrt2norm                         - --------------------------------------------------------------------------
%   sqrt2norm_r2018a                  - arguments
%   symaxis                           - --------------------------------------------------------------------------
%   symaxis_test                      - {
%   tblvertcat                        - --------------------------------------------------------------------------
%   tblvertcat_test                   - tblvertcat_test
%   testForGBFZSymmetryPlot           - 
%   toBPFZ                            - function nA_out = toBPFZ(qlist,nAlist)
%   tofiveFZ                          - --------------------------------------------------------------------------
%   tricollapse                       - --------------------------------------------------------------------------
%   vecpair2rmat                      - --------------------------------------------------------------------------
%   vert2con                          - - convert a set of points to the set of inequality constraints
%   vert2lcon                         - An extension of Michael Kleder's vert2con function, used for finding the 
%   write_video                       - close all %not sure if you need this
%   write_video_test                  - 
%   zeta_min2                         - --------------------------------------------------------------------------
%   zeta_min2_r2018a                  - arguments
