function [ypost,ypred,ysd,ytrue,ci,covmat,kfntmp,kfntmp2,mdl,gprMdl2list,X2] = egprm_fast(mdls,o2,ytrue,thr)
arguments
    mdls
    o2 = []
    ytrue = []
    thr(1,1) double = 1.1
end
% EGPRM  ensemble Guassian process regression mixture model with mixture threshold (thr) given GPR models (mdls)
K = length(mdls);

if iscell(mdls)
    mdls = [mdls{:}];
end

% unpack model components (X, y should be the same for all mdls)
X = mdls(1).mesh.ppts;
y = mdls(1).mesh.props;

X2 = cell(1,K);
for i = 1:K
    mdl = mdls(i);
    o2tmp = get_octpairs(o2,'oref',mdl.oref);
    o2tmp = normr(o2tmp);
    X2{i} = proj_down(o2tmp,mdl.projtol,mdl.usv,'zeroQ',mdl.zeroQ);
    % other gprmix inputs
    plotQ = false;
    dispQ = false;
    
    % GPR mixture model
    [ypredlist{i},ysdlist{i},cilist{i},covmatlist{i},kfntmplist{i},...
        kfntmp2list{i},gprMdl2list{i},ypredsigmoidlist{i},gprmMdllist{i}] = ...
        gprmix(mdl,X,y,X2{i},ytrue,thr,'plotQ',plotQ,'dispQ',dispQ);
end

%compile variables
ypredtmp = [ypredlist{:}];
ysdtmp = [ysdlist{:}];
citmp = [cilist{:}];
covmattmp = cat(3,covmatlist{:});

% covariance function mixture
for i = 1:K
    kfntmp = kfntmplist{i};
    kfntmp2 = kfntmp2list{i};
    ypredsigmoid = ypredsigmoidlist{i};
    gprmMdl = gprmMdllist{i};
    kfn{i} = @(ypredsigmoid) kfnmix(kfntmp,kfntmp2,ypredsigmoid,gprmMdl.scl,gprmMdl.thr);
end

% kfntmp = [kfntmplist{:}];
% kfntmp2 = [kfntmp2list{:}];

ypred = mean(ypredtmp,2);
ysd = mean(ysdtmp,2);
ci = mean(citmp,2);
covmat = mean(covmattmp,3);

covmat = nearestSPD(covmat);
covmat = nearestSPD(covmat); %twice to deal with numerical precision

l = max([zeros(size(ypred)),0.895255*ypred],[],2)-ypred;
u = 1.2444*ypred;
zerofloorQ = true;
n = 100;
ypost = tmvn(ypred,covmat,l,u,n,zerofloorQ);

end

%% CODE GRAVEYARD
%{
% X2 = mdls(1).data.ppts;
% ytrue = mdls(1).data.props;
%}