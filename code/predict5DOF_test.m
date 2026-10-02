%PREDICT5DOF_TEST  evaluate a GPR model trained by interp5DOF with predict5DOF
% Checks that predict5DOF returns one prediction per query GB and
% reproduces the predictions interp5DOF made for the same query GBs. Also
% prints how long each stage of predict5DOF takes (predict5DOF_timing) and
% how exact GPR prediction time grows with the number of training points.
%
% There are no %% sections on purpose: run-tests runs each section as a
% separate test in its own workspace.

rng(11)
ninputpts = 388;
npredpts = 500;

% random training and query GBs, BRK energies as the property
[~,qm,nA] = get_five(ninputpts);
[~,qm2,nA2] = get_five(npredpts);
y = GB5DOF_setup([],qm,nA);

% train, with exact predictions
[ypred,~,mdl] = interp5DOF(qm,nA,y,qm2,nA2,'gpr',1,...
    'mygpropts',{'PredictMethod','exact'},'dispQ',false);

% predict5DOF reproduces interp5DOF's predictions, one per query GB
ypred2 = predict5DOF(qm2,nA2,mdl,'noboundaryQ',false);
assert(isequal(size(ypred2),[npredpts 1]),...
    'predict5DOF returned %d predictions for %d query GBs',numel(ypred2),npredpts)
assert(max(abs(ypred2-ypred)) <= 1e-8*max(abs(ypred)),...
    'predict5DOF and interp5DOF predictions differ by up to %g',max(abs(ypred2-ypred)))

% uncertainty and the no-boundary constraint (predict5DOF's defaults)
[ypred3,ysd3,yint3] = predict5DOF(qm2,nA2,mdl,'qm',qm);
assert(isequal(size(ypred3),[npredpts 1]) && isequal(size(ysd3),[npredpts 1]) ...
    && isequal(size(yint3),[npredpts 2]),'unexpected predict5DOF output sizes')

% BP normals over a hemisphere at a fixed misorientation (Sigma3: 60 deg
% about [111]), plus low-index normals that sit on symmetry elements
[az,el] = meshgrid(linspace(0,2*pi,24),linspace(0,pi/2,7));
[n1,n2,n3] = sph2cart(az(:),el(:),1);
nA4 = [n1 n2 n3; normr([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 -1 0; 1 1 1; 1 1 -2; 1 2 3])];
qm4 = repmat(ax2qu([[1 1 1]/sqrt(3) pi/3]),size(nA4,1),1);
o4 = get_octpairs(five2oct(qm4,nA4),1,'oref',mdl.oref,'dispQ',false);
fprintf('Including ties in the VFZ mapping gives %d rows for %d hemisphere query GBs\n',...
    size(o4,1),size(qm4,1))
ypred4 = predict5DOF(qm4,nA4,mdl,'noboundaryQ',false);
assert(isequal(size(ypred4),[size(qm4,1) 1]),...
    'predict5DOF returned %d predictions for %d query GBs',numel(ypred4),size(qm4,1))

% time spent in each stage of predict5DOF
predict5DOF_timing(qm2,nA2,mdl);

% exact GPR prediction time vs. number of training points (synthetic 7D
% points, fixed hyperparameters, no fitting). The last two columns predict
% 10 points with sd, from the compact and the full model.
nquery = 1000;
Xq = normr(randn(nquery,7));
fprintf('exact GPR predict for %d query points (seconds)\n',nquery)
fprintf('%8s %10s %12s %16s %16s\n','N','mean','mean+sd','10 pts compact','10 pts full')
for N = [1000 2000 4000 8000]
    X = normr(randn(N,7));
    gpr = fitrgp(X,sum(X,2),'FitMethod','none','PredictMethod','exact',...
        'KernelParameters',[0.3;1],'Sigma',0.1);
    cgpr = compact(gpr);
    tic; ytmp = predict(cgpr,Xq); tmean = toc; %#ok<NASGU>
    tic; [ytmp,sdtmp] = predict(cgpr,Xq); tsd = toc; %#ok<ASGLU>
    tic; [ytmp,sdtmp] = predict(cgpr,Xq(1:10,:)); tsd10c = toc; %#ok<ASGLU>
    tic; [ytmp,sdtmp] = predict(gpr,Xq(1:10,:)); tsd10f = toc; %#ok<ASGLU>
    fprintf('%8d %10.3f %12.3f %16.3f %16.3f\n',N,tmean,tsd,tsd10c,tsd10f)
end
