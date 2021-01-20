function [ypred,ysd,ci,covmat,kfntmp,kfntmp2,gprMdl2] = gprmix(mdl,X,y,X2,ytrue,thr,nv)
arguments
    mdl
    X double = []
    y(:,1) double = []
    X2 double = []
    ytrue = []
    thr(1,1) double = 1.1
    nv.scl(1,1) double = 30
    nv.plotQ(1,1) logical = false
    nv.dispQ(1,1) logical = false
    nv.gprMdl2 = []
end
% GPRMIX  Gaussian process regression mixture model to better predict lower property values.
%--------------------------------------------------------------------------
% Inputs:
%  gprMdl - a Guassian process regression object
%  X - input points
%  y - input properties
%  X2 - output points
%  ytrue - true output properties
%  thr - property threshold for splitting the input data
%  nv.scl - higher values make the sigmoid function steeper
%  nv.plotQ - whether to do a parity plot
%  nv.dispQ - whether to display RMSE and MAE values
%
% Outputs:
%  ypred - predicted output properties of GPR mixture model
%  ysd - predicted output standard deviations of GPR mixture model
%
% Usage:
%  [ypred,ysd] = gprmix(gprMdl,X,y,X2);
%  [ypred,ysd] = gprmix(gprMdl,X,y,X2,[],thr);
%  [ypred,ysd] = gprmix(gprMdl,X,y,X2,ytrue,thr);
%
% Dependencies:
%  paperfigure.m (optional if nv.plotQ == false)
%  parityplot.m (optional if nv.plotQ == false)
%  get_errmetrics.m (optional if nv.dispQ == false)
%
% Notes:
%  *
%
% see also FITRGP
%
% Author(s): Sterling Baird
%
% Date: 2020-09-07
%--------------------------------------------------------------------------

%% setup
%unpack name-value pairs
scl = nv.scl;
plotQ = nv.plotQ;
dispQ = nv.dispQ;
gprMdl = mdl.gprMdl;

% account for case where no true values supplied
if isempty(ytrue)
    ytrue = nan(size(X2,1),1);
end

%% split input data with property values below threshold
ids = find(y < thr+0.1);
X = X(ids,:);
y = y(ids);

%% original model
%predict on all test points using original model
% orefs = get_orefs(8);

% if mdl.projQ
%     o2 = proj_up(X2,mdl.usv);
% else
%     o2 = X2;
% end
% for i = 1:K
% if i == 1
%     oref = get_ocubo(1,'random',[],10);
% else
%     oref = get_ocubo();
% end
% o2 = get_octpairs(o2,'oref',oref);
% end
[ypredtmp,ysdtmp,citmp] = predict(gprMdl,X2);

kfntmp = gprMdl.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl.Impl.ThetaHat);
covmattmp = kfntmp(X2,X2);
% [ypredtmpundoc,covmattmp,ciundoc] = predictExactWithCov(gprMdl.Impl,X2,alpha);

%find predictions below threshold and populate ypred
ids2 = ypredtmp <= thr;

ypred = zeros(size(X2,1),1);
ci = zeros(size(citmp));
ysd = zeros(size(ysdtmp));
covmat = zeros(size(covmattmp));

ypred(~ids2) = ypredtmp(~ids2);
ysd(~ids2) = ysdtmp(~ids2);
ci(~ids2,:) = citmp(~ids2,:);
covmat(~ids2) = covmattmp(~ids2);

%% lower subset model
%train on input points with properties less than threshold
if ~isempty(X)
    gprMdl2 = fitrgp(X,y,'PredictMethod','exact','KernelFunction','exponential');
    [ypredtmp2,ysdtmp2,citmp2] = predict(gprMdl2,X2);
    kfntmp2 = gprMdl2.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl2.Impl.ThetaHat);
    covmattmp2 = kfntmp2(X2,X2);
    % [ypredtmpundoc2,covmattmp2,ciundoc2] = predictExactWithCov(gprMdl2.Impl,X2,alpha);
elseif ~isempty(nv.gprMdl2)
    gprMdl2 = nv.gprMdl2;
    [ypredtmp2,ysdtmp2,citmp2] = predict(gprMdl2,X2);
    kfntmp2 = gprMdl2.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl2.Impl.ThetaHat);
    covmattmp2 = kfntmp2(X2,X2);
else
    [ypredtmp2,ysdtmp2] = deal([]);
end

% populate ypred with predictions above threshold
ypred(ids2) = ypredtmp2(ids2);
ysd(ids2) = ysdtmp2(ids2);
ci(ids2,:) = citmp2(ids2,:);
covmat(ids2) = covmattmp2(ids2);

ypredtmp3 = ypred;

%% GPR mixture model
% sigmoid function values
% sigfn = @(x,scl,xshift) 1./(1+exp(-scl*(x-xshift))); %changed to sigfn.m file
A = sigfn(ypred,scl,thr);
% x = [0,1.5];
% vq = [0,1];
% A = interp1(x,vq,ypred);
% A = 0.25;
B = 1-A;
% [A,B] = deal(0.5*ones(size(ypred)));

% sigmoid mixing
ypred = A.*ypredtmp+B.*ypredtmp2;
% if isempty(nv.gprMdl2)
ysd = A.*ysdtmp+B.*ysdtmp2;
ci = A.*citmp+B.*citmp2;
covmat = A.*covmattmp+B.*covmattmp2;

% % sigmoid mixing
% ypred = A.*ypredtmp+B.*ypredtmp2;
% % if isempty(nv.gprMdl2)
% ysd = A.*ysdtmp+B.*ysdtmp2;
% ci = A.*citmp+B.*citmp2;
% covmat = A.*covmattmp+B.*covmattmp2;
% end
% kfn = kmix(kfntmp,kfntmp2);

%% plotting & error metrics
if plotQ
%     paperfigure(2,2);
    % new parity plot
    multiparity({...
        struct('ytrue',ytrue,'ypred',ypredtmp),...
        struct('ytrue',ytrue,'ypred',ypredtmp2),...
        struct('ytrue',ytrue,'ypred',ypredtmp3),...
        struct('ytrue',ytrue,'ypred',ypred)})
    % sigmoid mixing function
%     nexttile(2)
%     cla reset
%     x = 0.5:0.01:1.5;
%     plot(x,sigfn(x,scl,thr),'LineWidth',1.5);
%     ylim([0 1])
%     xlabel('GBE ($J m^{-2}$)','Interpreter','latex')
%     ylabel('Mixing fraction (f)','Interpreter','latex')
%     axis square
%     papertext(2)
end
if dispQ
    get_errmetrics(ytrue,ypred,'dispQ',true);
end

end
%% CODE GRAVEYARD
%{
[ypred,ysd] = deal(zeros(size(ppts2,1),1));
% ysd(ids2) = ypredtmp2(ids2);
% ysd(~ids2) = ysdtmp(~ids2);

% Xsub = X2(ids2,:);
% ytruesub = ytrue(ids2);
%}