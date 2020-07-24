# Overview

This repository contains a set of tutorials for grain boundary octonion computations to accompany the papers: 

1. [Submitted] Chesser, I., Francis, T., De Graef, M., & Holm, E. A. (2019). Learning the grain boundary manifold: tools for grain boundary data representation and visualization. *Acta Materialia*. 

2. Francis, Toby, et al. "[A geodesic octonion metric for grain boundaries](https://www.sciencedirect.com/science/article/abs/pii/S1359645418309844)." *Acta Materialia* 166 (2019): 135-147.

High performance octonion computations have been implemented in [EMsoft](https://github.com/EMsoft-org/EMsoft)


# GB Data

The primary dataset for paper 1 is the canonical Olmsted dataset consisting of GB crystallography data for 388 GBs in metals with cubic point group symmetry. Grain boundary energy and mobility values were computed for each structure using the Foiles-Hoyt FCC Ni EAM  potential. 

3. Olmsted, David L., Stephen M. Foiles, and Elizabeth A. Holm. "Survey of computed grain boundary properties in face-centered cubic metals: I. Grain boundary energy." *Acta Materialia* 57.13 (2009): 3694-3703.

4. Olmsted, David L., Elizabeth A. Holm, and Stephen M. Foiles. "Survey of computed grain boundary properties in face-centered cubic metals—II: Grain boundary mobility." *Acta Materialia* 57.13 (2009): 3704-3713.


The following relevant files are found in the Data directory: 

**olmsted_xstal_info.csv**: consolidates crystallographic info, including grain orientations, CSL/DSC lattice vectors

**olm_properties.txt**: GB energy Ni ($J/m^2$), energy Al , mobility Ni $m/(s GPa)$, dissipation energy Ni (eV/atom)

**olm_octonion_list.txt**: octonions corresponding to GBs in Olmsted dataset with BP rotated to lie along z direction

**olm_pairwise_distances_cubic.txt**: pairwise distance matrix for Olmsted dataset, cubic (432) symmetry applied

# Tutorial Code

Written in MATLAB and/or Python: 

**Example 0**: <br/>
Convert traditional grain boundary representations to the grain boundary octonion (GBO) representation. 
(M,n) --> octonion
(O1,O2) --> octonion 

**Example 1**: <br/>
Compute the symmetrized GB octonion distance for a pair of grain boundaries with arbitrary point group symmetry

**Example 2**: <br/>
Compute a pairwise distance matrix for a grain boundary dataset

