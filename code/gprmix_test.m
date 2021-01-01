function [ypred,ysd,ytrue] = gprmix_test(fname)
arguments
   fname = 'gpr46883_gitID-233ae6e_puuID-5690bda0_kim-exponential-rng11.mat'
end
% GPRMIX_TEST  Perform sigmoid-based mixture with model trained on subset of low property values.
%% setup
% load data
addpathdir(fname)
S = load(fname);

% unpack model
mdl = S.mdl;

% unpack model components
gprMdl = mdl.gprMdl;
X = mdl.mesh.ppts;
y = mdl.mesh.props;
X2 = mdl.data.ppts;
ytrue = mdl.data.props;

% other gprmix inputs
thr = 1.1;
plotQ = true;
dispQ = true;

%% GPR mixture model
[ypred,ysd] = gprmix(gprMdl,X,y,X2,ytrue,thr,'plotQ',plotQ,'dispQ',dispQ);
end
