function [ypost,ypred,ysd,ytrue,ci,covmat,kfn,egprmMdl,gprMdl2list,X2,egprmMdlpars] = ...
    egprm(qm,nA,y,qm2,nA2,K,method,epsijk,NV)
arguments
    qm = [] %input misorientation quaternions
    nA = [] %input BP normals
    y(:,1) = [] %property values
    qm2 = [] %query misorientations
    nA2 = [] %query BP normals
    K(1,1) double = 10 % number of ensemble components, make == 1 for 1-component ensemble (i.e. no ensemble)
    method char {mustBeMember(method,{'gpr','egpr','gprm','egprm'})} = 'egprm'
    epsijk(1,1) double = 1
    NV.pgnum(1,1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
    NV.databary = [] %for use with bary methods
    NV.facetIDs = [] %for use with bary methods
    NV.ytrue = [] %user-specified "true" values for error calculations
    NV.modelparsspec = struct()
    NV.brkQ(1,1) logical = true %whether to compute BRK values as ytrue
    NV.mygpropts = {'PredictMethod','exact'} %for use with gpr methods 'gpr' or 'sphgpr'
    NV.r double = [] %for use with 'idw' method, alternatively set to [] for automatic estimation
    NV.uuid(1,8) char = get_uuid() %unique ID associated with this interpolation run
    NV.o = [] %input octonions, specify these or qm/nA pairs
    NV.o2 = [] %query octonions, specify these or qm2/nA2 pairs
    NV.thr(1,1) double = 1.1 %threshold for GPR mixture
    NV.scl(1,1) double = 30 %controls width/steepness of mixing sigmoid function
    NV.mdls = [] %models from ensembleVFZO.m can be supplied to make it faster
    NV.egprmMdl = [] %an egprmMdl can be supplied to make the interpolation faster (skip symmetrization)
    NV.mixQ(1,1) logical = true %whether or not to perform the "mixture" part of EGPRM
    NV.postQ(1,1) logical = false %whether to sample from posterior distribution (expensive for e.g. 10000+ GBs)
    NV.sig(1,1) double = 0
    NV.dispQ(1,1) logical = false
    NV.egprmDispQ(1,1) logical = true
    NV.KdispQ(1,1) logical = true
    NV.covK(1,1) double {mustBeInteger} = 1
end
% EGPRM  ensemble Guassian process regression mixture model with mixture threshold (thr) given GPR models (mdls)
egprm_starttime = tic;
thr = NV.thr;
scl = NV.scl;
ytrue = NV.ytrue;
mdls = NV.mdls;
o2 = NV.o2;
pgnum = NV.pgnum;
egprmMdl = NV.egprmMdl;
if ~isempty(egprmMdl)
    mixQ = egprmMdl.mixQ;
    K = egprmMdl.K;
    thr = egprmMdl.thr;
    scl = egprmMdl.scl;
    cgprMdls2 = egprmMdl.cgprMdls2;
    covK = egprmMdl.covK;
else
    mixQ = NV.mixQ;
    cgprMdls2 = [];
    covK = NV.covK;
end
postQ = NV.postQ;
brkQ = NV.brkQ;
sig = NV.sig;
egprmDispQ = NV.egprmDispQ;
KdispQ = NV.KdispQ;

if contains(method,'gpr')
    ensembleMethod = 'gpr';
end

if egprmDispQ
    disp(['K: ' int2str(K) ', mixQ: ' int2str(mixQ)])
end

freshQ = false;
% perform the ensemble
if isempty(mdls) && isempty(egprmMdl)
    NVtmp = rmfield(NV,{'mdls','thr','egprmMdl','scl','mixQ','postQ','egprmDispQ'});
    NVpairs = namedargs2cell(NVtmp);
    [~,~,~,mdls,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,ensembleMethod,NVpairs{:});
    freshQ = true;
elseif ~isempty(egprmMdl)
    mdls = egprmMdl.mdls;
end
if iscell(mdls)
    mdls = [mdls{:}];
end

% K = length(mdls);

% unpack model components (y should be the same for all mdls)
y = mdls(1).mesh.props;

%initialize
[X2,ypredlist,ysdlist,cilist,covmatlist,kfntmplist,kfntmp2list,...
    gprMdl2list,ypredsigmoidlist,gprmMdllist] = deal(cell(1,K));

if isempty(o2)
    o2 = five2oct(qm2,nA2);
end

if isempty(egprmMdl) && brkQ
    ytrue = GB5DOF_setup(o2(:,1:4),o2(:,5:8),[0 0 1],'Ni',epsijk);
end

%Ensemble

for i = 1:K
    if KdispQ
        disp(' ')
        disp([int2str(i) '-th ensemble component symmetrization'])
    end
    mdl = mdls(i);
    if freshQ
        X2{i} = mdl.data.ppts;
    else
        o2tmp = get_octpairs(o2,'oref',mdl.oref,'pgnum',pgnum);
        o2tmp = normr(o2tmp);
        X2{i} = proj_down(o2tmp,mdl.projtol,mdl.usv,'zeroQ',mdl.zeroQ);
    end
    if isfield(mdl,'cgprMdl')
        gprMdl = mdl.cgprMdl;
    else
        gprMdl = mdl.gprMdl;
    end
    [ypredlist{i},ysdlist{i},cilist{i}] = predict(gprMdl,X2{i});
    kfntmplist{i} = gprMdl.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl.Impl.ThetaHat);
    covmatlist{i} = kfntmplist{i}(X2{i},X2{i});
end

if isempty(egprmMdl)
    gprMdl2list = [];
else
    if isfield(egprmMdl,'gprMdl2list')
        gprMdl2list = egprmMdl.gprMdl2list;
    else
        gprMdl2list = egprmMdl.cgprMdls2;
    end
end
%compile variables
ypredtmp = [ypredlist{:}];
ysdtmp = [ysdlist{:}];
citmp = cat(3,cilist{:});
covmattmp = cat(3,covmatlist{:});
clear covmatlist

% Homogenize (i.e. combine) the ensemble components
homtype = 'median'; %homogenization type
switch homtype
    case 'normal'
        num = (1:K).';
        pd = fitdist(num,'normal');
        ypdf = pdf(pd,num);
        ypdf = ypdf./sum(ypdf); %sum to 1
        fn = @(ypred,dim) sum(ypred.*(ypdf.'),dim); %mean(ypred.*ypdf,dim);
        fn2 = @(mat,dim) mean(pagemtimes(mat,reshape(ypdf,1,1,K)),dim);
    case 'mean'
        fn = @(ypred,dim) mean(ypred,dim);
        fn2 = fn;
    case 'median'
        fn = @(ypred,dim) median(ypred,dim);
        fn2 = fn;
    case 'min'
        fn = @(x,y) min(x,[],y);
        fn2 = fn;
    case 'max'
        fn = @(x,y) max(x,[],y);
        fn2 = fn;
end
yprede = fn(ypredtmp,2);
ysde = fn(ysdtmp,2);
cie = fn2(citmp,3);
covmate = fn2(covmattmp,3);
clear covmattmp %10000x10000x10 array = ~8 GB
kfne = @(X1,X2) kfnavg(kfntmplist,X1,X2);

if ~mixQ
    ypred = yprede;
    ysd = ysde;
    ci = cie;
    covmat = covmate;
    kfn = kfne;
else
    %% GPR Mix Setup
    % lower subset models
    [ypredlist2,ysdlist2,cilist2,covmatlist2,gprMdl2listtmp,X] = deal(cell(1,K));
    for i = 1:K
        %train on input points with properties less than threshold
        mdl = mdls(i);
        if isempty(gprMdl2list) && isempty(cgprMdls2)
            thr2 = thr+0.1;
            ids = find(yprede <= thr2);

            ysub = y(ids);
            X{i} = mdl.mesh.ppts(ids,:);
            gprMdl2 = fitrgp(X{i},ysub,'PredictMethod','exact','KernelFunction','exponential');
            gprMdl2listtmp{i} = gprMdl2;
        else
            %load from either cgprMdl or gprMdl
            if isempty(gprMdl2list)
                unpacklist = cgprMdls2;
            else
                unpacklist = gprMdl2list;
            end
            if iscell(unpacklist)
                gprMdl2 = unpacklist{i};
            else
                gprMdl2 = unpacklist;
            end
        end
        [ypredlist2{i},ysdlist2{i},cilist2{i}] = predict(gprMdl2,X2{i});
        kfntmp2 = gprMdl2.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl2.Impl.ThetaHat);
        covmatlist2{i} = kfntmp2(X2{i},X2{i});
    end
    gprMdl2list = gprMdl2listtmp;
    %compile variables
    ypredtmp2 = [ypredlist2{:}];
    ysdtmp2 = [ysdlist2{:}];
    citmp2 = cat(3,cilist2{:});
    covmattmp2 = cat(3,covmatlist2{:});
    clear covmatlist2
    
    % homogenize the lower subset models (e2 stands for ensemble of lower subsets)
    yprede2 = fn(ypredtmp2,2);
    ysde2 = fn(ysdtmp2,2);
    cie2 = fn2(citmp2,3);
    covmate2 = fn2(covmattmp2,3);
    clear covmattmp2
    kfne2 = @(X1,X2) kfnavg(kfntmplist2,X1,X2);
    
    %find predictions below threshold and populate ypred
    ids2 = yprede <= thr;
    ypredsigmoid = zeros(size(X2{1},1),1);
    ypredsigmoid(~ids2) = yprede(~ids2);
    
    % populate ypred with predictions above threshold
    ypredsigmoid(ids2) = yprede2(ids2);
    
    %% GPR mixture
    % sigmoid function values
    A = sigfn(ypredsigmoid,scl,thr);
    B = 1-A;
    
    % sigmoid mixing
    ypred = A.*yprede+B.*yprede2;
    ysd = A.*ysde+B.*ysde2;
    ci = A.*cie+B.*cie2;
    
    % sigmoid mix covariance matrices
    Amat = 1/2*(A+A.');
    Bmat = 1/2*(B+B.');
    covmat = Amat.*covmate+Bmat.*covmate2;
    
    % sigmoid mix covariance function
    kfn = @(X1,X2,ypredsigmoid) kfnmix(kfne(X1,X2),kfne2(X1,X2),ypredsigmoid,scl,thr);
end
egprm_runtime = toc(egprm_starttime);

posterior_starttime = tic;
if postQ
    % get nearest symmetric positive definite matrix
    covmat = nearestSPD(covmat); %also takes a while
    covmat((covmat > -1e-12) & (covmat < 0)) = 1e-12;
%     covmat = nearestSPD(covmat); %twice to deal with numerical precision, small negative numbers
    %alternative to second nearestSPD command, could do covmat(covmat < 0) = 1e-12; or similar
    
    %% Posterior Sampling
    l = max([zeros(size(ypred)),0.895255*ypred],[],2)-ypred; %lower bound
    u = 1.2444*ypred; %upper bound
    zerofloorQ = true;
    n = 100; %number of samples from posterior distribution
    
    ypost = tmvn(ypred,covmat,l,u,n,zerofloorQ); %takes a long time for 10000^2 matrix
else
    ypost = [];
    l = [];
    u = [];
    zerofloorQ = [];
    n = [];
end
posterior_runtime = toc(posterior_starttime);

%% Package EGPRM Model
if isempty(egprmMdl)
    p = gcp;
    cores = p.NumWorkers;
    
    mesh = mdls(1).mesh;
    projQ = mdls(1).projQ;
    projtol = mdls(1).projtol;
    usv = mdls(1).usv;
    zeroQ = mdls(1).zeroQ;
    oref = mdls(1).oref;
    oreflist = vertcat(mdls.oref);
    
    egprmMdl = var_names(ypred,ysd,ytrue,ci,covmat,l,u,zerofloorQ,n,ypost,...
        thr,scl,mdls,o2,mesh,oref,oreflist,projQ,projtol,...
        usv,zeroQ,gprMdl2list,mixQ,K,covK,cores,pgnum,sig,brkQ,...
        posterior_runtime,egprm_runtime);
    method = 'gpr';
    if mixQ
        method = [method 'm'];
    end
    if K >= 2
        method = ['e' method];
    end
    egprmMdl.method = method;
    egprmMdl.interpfn = @() egprm();
    egprmMdl.mdlcmd = @() egprm();
    
    egprmMdlpars = var_names(ypred,ysd,ytrue,ci,l,u,zerofloorQ,n,ypost,...
       thr,scl,o2,mesh,oref,oreflist,projQ,projtol,usv,...
        zeroQ,mixQ,K,covK,cores,pgnum,sig,brkQ,...
        posterior_runtime,egprm_runtime);
    
    % deal with pure 'gpr' case (for e.g. tunnelplot.m compatibility)
    if ~mixQ && (K == 1)
        cgprMdl = mdls(1).cgprMdl;
        if isfield(egprmMdl,'gprMdl')
            egprmMdl.gprMdl = mdls(1).gprMdl;
            egprmMdl.cgprMdl = cgprMdl;
            egprmMdlpars.cgprMdl = cgprMdl;
        else
            egprmMdl.cgprMdl = cgprMdl;
            egprmMdlpars.cgprMdl = cgprMdl;
        end
    end
    
    %remove memory-intensive variables from "parameter" structure
    if isfield(mdls,'gprMdl')
        mdls = rmfield(mdls,'gprMdl');
    end
    if isfield(mdls,'interpfn')
        %saving interpfn is significant contrary to what `whos interpfn` shows
        mdls = rmfield(mdls,'interpfn');
    end
    egprmMdlpars.mdls = mdls;
    
    egprmMdlpars.method = method;
    
    if mixQ
        cgprMdls2 = cellfun(@compact,gprMdl2list,'UniformOutput',false);
        egprmMdl.cgprMdls2 = cgprMdls2;
        egprmMdlpars.cgprMdls2 = cgprMdls2;
    end
    
end
end

%% CODE GRAVEYARD
%{
% X2 = mdls(1).data.ppts;
% ytrue = mdls(1).data.props;
% kfntmp = [kfntmplist{:}];
% kfntmp2 = [kfntmp2list{:}];

%         t = n2c(ypredtmp.');
%         pd = cellfun(@(t) fitdist(t,'normal'),t);
%         mu = [pd.mu].';
%         sigma = [pd.sigma].';
%         ypdf = normpdf(ypredtmp,mu,sigma);
%         ypdf = ypdf./sum(ypdf,2); %each row sums to 1
%         fn = @(ypred,dim) sum(ypdf.*ypred,dim);


% % GPR mixture model
% [ypred,ysd,ci,covmat,kfntmp,kfntmp2,gprMdl2,ypredsigmoid,gprmMdl] = ...
%     gprmix(mdl,X,y,X2{i},ytrue,thr,'plotQ',plotQ,'dispQ',dispQ,'gprMdl2',gprMdl2);


% GPR mixture model
% sigmoid function values
A = sigfn(ypred);
B = 1-A;
% sigmoid mixing
ypred = A.*ypredtmp+B.*ypredtmp2;
ysd = A.*ysdtmp+B.*ysdtmp2;
ci = A.*citmp+B.*citmp2;
covmat = A.*covmattmp+B.*covmattmp2;

ypredsigmoid,kfntmplist,kfntmp2list

% ci = zeros(size(citmp));
% ysd = zeros(size(ysdtmp));
% covmat = zeros(size(covmattmp));

% cie = zeros(size(ypredtmp,1),2);
% cie(:,1) = fn(citmp(:,1),2);
% cie(:,2) = fn(citmp(:,2),2);

% cie2 = zeros(size(ypredtmp2,1),2);
% cie2(:,1) = fn(citmp2(:,1),2);
% cie2(:,2) = fn(citmp2(:,2),2);


% % covariance function mixture
% for i = 1:K
%     kfntmp = kfntmplist{i};
%     kfntmp2 = kfntmp2list{i};
%     ypredsigmoid = ypredsigmoidlist{i};
%     gprmMdl = gprmMdllist{i};
%     kfn{i} = @(ypredsigmoid) kfnmix(kfntmp,kfntmp2,ypredsigmoid,gprmMdl.scl,gprmMdl.thr);
% end

% o = mdls(1).mesh.pts;
% X = mdls(1).mesh.ppts;



    %     if ~isempty(egprmMdl)
    % %         gprMdl2 = egprmMdl.gprMdl2list{i};
    %     else
    %         gprMdl2 = [];
    %     end



            % osub = o(ids,:);
                        % ysub = mdls(1).mesh.props(ids,:); %should be same across all mdls
            %         otmp = get_octpairs(sqrt2norm(osub),'oref',mdl.oref);
            %         otmp = normr(otmp);
            %         X{i} = proj_down(otmp,mdl.projtol,mdl.usv,'zeroQ',mdl.zeroQ);


    %     % GPR mixture model
    %     [ypredlist{i},ysdlist{i},cilist{i},covmatlist{i},kfntmplist{i},...
    %         kfntmp2list{i},gprMdl2list{i},ypredsigmoidlist{i},gprmMdllist{i}] = ...
    %         gprmix(mdl,X,y,X2{i},ytrue,thr,'plotQ',plotQ,'dispQ',dispQ,'gprMdl2',gprMdl2);


    % other gprmix inputs
    plotQ = false;
    dispQ = false;


% o = NV.o;
%}