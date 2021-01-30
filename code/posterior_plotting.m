% posterior_plotting
%% setup
% files = dir(fullfile('**','interp5DOF-paper','figures'));
% figfolder = files(1).folder;

addpath(genpath('.'))

set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'defaultAxesFontSize',12)

figfolder = 'C:\Users\sterg\Documents\GitHub\posterior-sampling\figures';

%% test function truncated multi-variate normal distribution sampling
testPosteriorCovarianceGP
fname = 'testfn-mvn-tmvn';
savefigpng(figfolder,fname);

