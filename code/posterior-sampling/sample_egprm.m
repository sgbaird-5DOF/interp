function [ypost,postMdls,o2] = sample_egprm(mdl,npts2,npostpts,l,u,K,nv)
arguments
    mdl struct
    npts2(1,1) double = 10000
    npostpts(1,1) double = 100
    l(:,1) double = [] %lower bound
    u(:,1) double = [] %upper bound
    K(1,1) double = 10
    nv.method char {mustBeMember(nv.method,{'mvrandn','slicesample'})} = 'mvrandn'
    nv.zerofloorQ(1,1) logical = true
    nv.nearestSPD_Q(1,1) logical = true
end

disp('posterior sampling')

%get new points that are uniformly distributed in 5DOF space
o2 = get_ocubo(npts2); %posterior points

[~,mu,~,~,~,covmat] = egprm('o2',o2,'egprmMdl',mdl);

%parse lower and upper bounds
if isempty(l)
    l = 0.895255*mu;
elseif isscalar(l)
    l = mu*l;
end
if isempty(u)
    u = 1.2444*mu;
elseif isscalar(u)
    u = mu*u;
end
l = l-mu;
u = u-mu;

%% Posterior Sampling
nvpairs = namedargs2cell(nv);
ypost = tmvn(mu,covmat,l,u,npostpts,nvpairs{:}); %takes a long time for 10000^2 matrix

postMdls = cell(1,npostpts);
for i = 1:npostpts
    yp = ypost(:,i);
    postMdls{i} = egprm([],[],yp,[],[],K,'o',o2,'mixQ',false);
end

end

%% CODE GRAVEYARD
%{
% l = zeros(size(mu));
% u = inf*ones(size(mu));


    l(:,1) double = max([zeros(size(mu)),0.895255*mu],[],2)-mu; %lower bound
    u(:,1) double = 1.2444*mu; %upper bound

egprmMdltmp = egprmMdl;

% y3 = egprm();
% X3 = proj_down(o2,projtol,usv);
% covmat = mdl.kfn(o2,o2); %assumes use of reference octonions from oreflist

%max([zeros(size(mu)),0.895255*mu],[],2)-mu;

%     ypost = mvnrnd_trn(l.',u.',mu.',covmat,n);

% starttime = tic;
% runtime = toc(starttime);

%}