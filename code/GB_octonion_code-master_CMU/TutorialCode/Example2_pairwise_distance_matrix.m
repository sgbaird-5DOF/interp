
%% 1. LOAD OCTONION UTILITY FUNCTIONS 

addpath('crystal_symmetry_ops')
addpath('octonion_functions/')

%% 2. PROBLEM STATEMENT 

% Given a list of GB's, construct a pairwise distance matrix representing the connectivity of the dataset. Each 
% entry is the geodesic angle between a pair of grain boundaries o_i and o_j.  

% We will use the GBpd function to compute the pairwise distance matrix for a subset of the Olmsted dataset
% with cubic symmetry

% INPUT: 

% 1/2. oldpd, oldoct
% These will be set to empty lists in this example and are only relevant if
% we want to add to an already existing (partially computed) pairwise
% distance matrix 

% 3. newoct
% List of octonions to use for pairwise distance matrix construction (or
% expansion of a previously computed matrix). 

% 4. pgnum 
% point group number, names 

% pgnum: number of point group symmetry operators (from 1 to 32, cubic is 30, names given in below and crystal_symmetry_ops/PGnames.mat)

%     '1: 1'    '2: -1'    '3: 2'    '4: m'    '5: 2/m'    '6: 222'    '7: mm2'    '8: mmm'    '9: 4'    '10: -4'    '11: 4/m'    '12: 422'
%     '13: 4mm'    '14: -42m'    '15:4/mmm'    '16: 3'    '17: -3'    '18: 32'    '19: 3m'    '20: -3m'    '21: 6'    '22: -6'
%     '23:  6/m'    '24:  622'    '25: 6mm'    '26: -6m2'    '27: 6/mmm'    '28: 23'    '29: m3'    '30: 432'    '31: -43m'    '32: m-3m'

% 5. printbool, fname
% true if you want to output pairwise distance matrix to a text file named
% fname (include .txt in name)


%OUTPUT: 
% 1.pdtest
% pairwise distance matrix 


oldpd = []; oldoct = []; 
printbool = true; fname = 'little_test.txt';
pgnum = 30; 

test = importdata('../Data/olm_octonion_list.txt',' ',1); %list of GB octonions with number of octonions in file at top
data0 = test.data;
newoct = data0(1:3,:);


tic 
pdtest = GBpd(oldpd,oldoct,newoct,pgnum,printbool,fname);
toc 

% A faster implementation of this function can be found in the EMSOFT program EMGBOdm
