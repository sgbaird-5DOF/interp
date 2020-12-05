%ENSEMBLEVFZO_TEST  
addpathdir({'cu2qu.m','get_uuid.m','gmat2q.m'})
n = 10000;
n2 = 10000;
[five,qm,nA] = get_five(n); %input
[five2,qm2,nA2] = get_five(n2); %prediction
y = GB5DOF_setup(five);
ytrue = GB5DOF_setup(five2);
K = 10;
[ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,'gpr');

errmetrics = get_errmetrics(ypred,ytrue,'dispQ',true);

paperfigure()

save('ensemble-interp.mat')