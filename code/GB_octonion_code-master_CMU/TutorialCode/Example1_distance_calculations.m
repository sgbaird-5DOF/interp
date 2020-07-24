
%% 1. LOAD OCTONION UTILITY FUNCTIONS 

addpath('crystal_symmetry_ops')
addpath('octonion_functions/')

%% 2. PROBLEM STATEMENT 

% Given a list of GB pairs, can we compute the geodesic angle between each pair? 

% We will use the GBdist function to compute the distance between ~200 GB's
% with cubic symmetry

%INPUT: 
% data: an N x 16 matrix of GB octonion pairs for the distance calculation
% a single octonion pair should have the form: (o1,o2) = (qA,qB,qC,qD), each of the four quaternions should be normalized

% pgnum: number of point group symmetry operators (from 1 to 32, cubic is 30, names given in below and crystal_symmetry_ops/PGnames.mat)

%     '1: 1'    '2: -1'    '3: 2'    '4: m'    '5: 2/m'    '6: 222'    '7: mm2'    '8: mmm'    '9: 4'    '10: -4'    '11: 4/m'    '12: 422'
%     '13: 4mm'    '14: -42m'    '15:4/mmm'    '16: 3'    '17: -3'    '18: 32'    '19: 3m'    '20: -3m'    '21: 6'    '22: -6'
%     '23:  6/m'    '24:  622'    '25: 6mm'    '26: -6m2'    '27: 6/mmm'    '28: 23'    '29: m3'    '30: 432'    '31: -43m'    '32: m-3m'

% genplot: plotting boolean to generate histogram of geodesic distances / U(1) angles 

%OUTPUT: 
% omega_new: Nx1 vector, minimum distance (geodesic distance) computed for each input GB pair
% oct_new: Nx16 matrix, symmetrized octonions that give minimum geodesic distance
% zeta_new: Nx1 vector, minimizing U(1) angle


test = importdata('../Data/olm_octonion_list.txt',' ',1); %list of GB octonions with number of octonions as first line in file
data0 = test.data;

% as a simple example, we will fold the Olmsted dataset in half to give a 194x16 matrix of GB pairs
data = zeros(388/2,16);
data(:,1:8) = data0(1:194,:);
data(:,9:16) = data0(195:end,:);

pgnum = 30; %cubic symmetry
genplot = false;

%this takes 48.2 seconds on a single core (which is quite slow). A much faster implementation of a
%similar routine can be found in the EMSOFT programs EMGBO / EMGBOdm

tic 
[omega_test, oct_test, zeta_test] = GBdist(data, pgnum, genplot);
toc 

