% posterior_plotting
%% setup
% files = dir(fullfile('**','interp5DOF-paper','figures'));
% figfolder = files(1).folder;

addpath(genpath('.'))

setlatex()

figfolder = 'C:\Users\sterg\Documents\GitHub\posterior-sampling\figures';

%% test function truncated multi-variate normal distribution sampling
testPosteriorCovarianceGP()
fname = 'testfn-mvn-tmvn';
savefigpng(figfolder,fname);