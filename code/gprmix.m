function [ypred,ysd] = gprmix(gprMdl,X,y,X2,ytrue,thr,nv)
arguments
    gprMdl
    X double
    y(:,1) double
    X2 double
    ytrue = []
    thr(1,1) double = 1.1
    nv.scl(1,1) double = 30
    nv.plotQ(1,1) logical = false
    nv.dispQ(1,1) logical = false
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

% account for case where no true values supplied
if isempty(ytrue)
    ytrue = nan(size(X2,1),1);
end

%% split input data with property values below threshold
ids = find(y < thr);
X = X(ids,:);
y = y(ids);

%% original model
%predict on all test points using original model
[ypredtmp,ysdtmp] = predict(gprMdl,X2);

%find predictions below threshold and populate ypred
ids2 = ypredtmp < thr;
ypred = zeros(size(X2,1),1);
ypred(~ids2) = ypredtmp(~ids2);

%% lower subset model
%train on input points with properties less than threshold
if ~isempty(X)
gprMdl2 = fitrgp(X,y,'PredictMethod','exact','KernelFunction','exponential');
[ypredtmp2,ysdtmp2] = predict(gprMdl2,X2);
else
    [ypredtmp2,ysdtmp2] = deal([]);
end

% populate ypred with predictions above threshold
ypred(ids2) = ypredtmp2(ids2);
ypredtmp3 = ypred;

%% GPR mixture model
% sigmoid function values
sigfn = @(x,scl,xshift) 1./(1+exp(-scl*(x-xshift)));
A = sigfn(ypred,scl,thr);
B = 1-A;

% sigmoid mixing
ypred = A.*ypredtmp+B.*ypredtmp2;
ysd = A.*ysdtmp+B.*ysdtmp2;

%% plotting & error metrics
if plotQ
%     paperfigure(2,2);
    % new parity plot
    multiparity({...
        struct('ytrue',ytrue,'ypred',ypred),...
        struct('ytrue',ytrue,'ypred',ypredtmp3),...
        struct('ytrue',ytrue,'ypred',ypredtmp2),...
        struct('ytrue',ytrue,'ypred',ypredtmp)})
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