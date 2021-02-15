function [ypred,ytrue,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO_test(n,n2,K,method,epsijk)
arguments
    n(1,1) double = 1000
    n2(1,1) double = 1000
    K(1,1) double = 10
    method char = 'gpr'
    epsijk(1,1) double = 1
end
%ENSEMBLEVFZO_TEST  run a K-component ensemble for n input and n2 prediction GBs using the specified method
% addpathdir({'cu2qu.m','get_uuid.m','gmat2q.m'}) %use addpath(genpath('.')) in a top-level folder instead
[five,qm,nA] = get_five(n); %input
[five2,qm2,nA2] = get_five(n2); %prediction
y = GB5DOF_setup([],qm,nA,'Ni',epsijk);
ytrue = GB5DOF_setup([],qm2,nA2,'Ni',epsijk);
[ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,method);

errmetrics = get_errmetrics(ypred,ytrue,'dispQ',true);

%% CODE GRAVEYARD
%{
% y = GB5DOF_setup(five);
% ytrue = GB5DOF_setup(five2);
% n = 1000;
% n2 = 1000;
%}