%ENSEMBLEVFZO_TEST  
addpathdir({'cu2qu.m','get_uuid.m','gmat2q.m'})
n = 1000;
n2 = 1000;
[five,qm,nA] = get_five(n); %input
[five2,qm2,nA2] = get_five(n2); %prediction
epsijk = 1;
y = GB5DOF_setup([],qm,nA,'Ni',epsijk);
ytrue = GB5DOF_setup([],qm2,nA2,'Ni',epsijk);
K = 1;
[ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,'gpr');

errmetrics = get_errmetrics(ypred,ytrue,'dispQ',true);

paperfigure();
parityplot(ytrue,ypred);

% save('ensemble-interp.mat')

%% CODE GRAVEYARD
%{
% y = GB5DOF_setup(five);
% ytrue = GB5DOF_setup(five2);
%}