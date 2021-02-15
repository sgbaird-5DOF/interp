function [mdl,ntunnelpts,extend,A,B] = egprm_test(n,n2,K,thr,scl,epsijk,nv)
arguments
    n(1,1) double = 1000 %or 10000
    n2(1,1) double = 1000
    K(1,1) double = 10 %number of ensemble components
    thr(1,1) double = 1.1 %threshold for GPR mixture
    scl(1,1) double = 30
    epsijk(1,1) double = 1
    nv.mixQ(1,1) logical = true
    nv.loadQ(1,1) logical = true
    nv.seed(1,1) double = 10
    nv.mat char = 'Ni'
    nv.lgdloc char = 'best'
    nv.mdl = []
end
% EGPRM_TEST  perform a self-contained ensemble Gaussian process mixture run
loadQ = nv.loadQ;
seed = nv.seed;
mat = nv.mat;
mixQ = nv.mixQ;
lgdloc = nv.lgdloc;
mdl = nv.mdl;

%% setup
rng(seed);
[~,qm,nA] = get_five(n);
[~,qm2,nA2] = get_five(n2);
o = five2oct(qm,nA,epsijk);
o2 = five2oct(qm2,nA2,epsijk);
y = GB5DOF_setup(o(:,1:4),o(:,5:8),[0 0 1],mat,epsijk);
ytrue = GB5DOF_setup(o2(:,1:4),o2(:,5:8),[0 0 1],mat,epsijk);

%% load ensemble components
if loadQ
    load(['C:\Users\sterg\Documents\GitHub\posterior-sampling\figures\ensemble-interp-' int2str(n2)],...
        'mdllist');
end
%% perform EGPRM or EGPR
[ypost,ypred,ysd,ytrue,ci,covmat,kfn,mdl,gprMdl2list,X2] = ...
    egprm(qm,nA,y,qm2,nA2,K,'gpr','ytrue',ytrue,'mdls',mdllist,'thr',thr,...
    'scl',scl,'mixQ',mixQ,'egprmMdl',mdl);

%% filenaming
% method = 'egprm';
% fname = [method '-' int2str(n2)];

%% save
% save(fname,'mdl','-v7.3')

%% load
% load(fname,'mdl')

%% tunnel parameters
ntunnelpts = 300;

%% initial tunnelplot to get A,B
% paperfigure();
% [tpredlist,tsdlist,propList,methodlist,A,B] = tunnelplot(mdl,[],[],ntunnelpts);
pts = mdl.mdls(1).mesh.pts;
[A,B] = getAB(pts);

%% final parityplot and tunnelplot
% paperfigure(1,2,7.63916666666667);
% nexttile
% parityplot(ytrue,mdl.ypred);
% nexttile
extend = 0.5;
tunnelplot(mdl,A,B,ntunnelpts,'extend',extend,'nnQ',false,'nnQ2',false,'lgdloc',lgdloc);
% axis normal
