function [ypred,ysd,ytrue,ci,covmat,kfntmp,kfntmp2,mdl,gprMdl2] = gprmix_test(fname)
arguments
   fname = 'gpr46883_gitID-5e6bd14_puuID-2e34b447_kim-paper-data3.mat' %'gpr46883_gitID-b473165_puuID-50ffdcf6_kim-rng11.mat' %gpr46883_gitID-9818307_puuID-e5d12412_kim-rng11.mat
end
% GPRMIX_TEST  Perform sigmoid-based mixture with model trained on subset of low property values.
%% setup
% load data
addpathdir(fname)
S = load(fname);

% unpack model
mdl = S.mdl;

% unpack model components
% gprMdl = mdl.gprMdl;
X = mdl.mesh.ppts;
y = mdl.mesh.props;
X2 = mdl.data.ppts;
ytrue = mdl.data.props;

% other gprmix inputs
thr = 1.1;
plotQ = true;
dispQ = true;

%% GPR mixture model
[ypred,ysd,ci,covmat,kfntmp,kfntmp2,gprMdl2] = gprmix(mdl,X,y,X2,ytrue,thr,'plotQ',plotQ,'dispQ',dispQ);
end
