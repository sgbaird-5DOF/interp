%sample_egprm_test
epsijk = 1;
K = 2;
loadQ = false;
if loadQ
    load('data\egprm\pcombs\egprm1000_gitID-79cf903_puuID-ae3ad8e6_brk.mat','mdl','mdlpars')
else
    ninputpts = 5000;
    npredpts = 1000;
    [~,~,mdl] = egprm_setup(ninputpts,npredpts,'K',K);
end

npts2 = 5000;
npostpts = 3;
l = 0.895255;
u = 1.2444;
[ypost,postMdls,o2] = sample_egprm(mdl,npts2,npostpts,l,u,K);

% evaluate the models that are fitted to the sampled posterior points
% npredpts = 1000;
% [~,qm2,nA2] = get_five(npredpts);
o = sqrt2norm(mdl.mesh.pts,'quat');
[qm2,nA2] = oct2five(o,epsijk);
ypred = nan(npts2,npostpts);
for i = 1:1 %npostpts
    postMdl = postMdls{i};
    [~,ypred(:,i)] = egprm([],[],[],qm2,nA2,'egprmMdl',postMdl);
end

paperfigure();
parityplot(mdl.mesh.props,ypred(:,i))

paperfigure();
tunnelplot(postMdls,'nnQ',false,'nnQ2',false,'extend',0.5)

save('posterior-sampling\sample_egprm_test.mat')