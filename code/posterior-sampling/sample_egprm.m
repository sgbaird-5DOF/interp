function [ypost, runtime] = sample_egprm(mdl,l,u,n,nv)
arguments
    mdl struct
    l(:,1) double = max([zeros(size(mu)),0.895255*mu],[],2)-mu; %lower bound
    u(:,1) double = 1.2444*mu; %upper bound
    n(1,1) double = 100
    nv.nearestQ(1,1) logical = true
    nv.method char = 'mvrandn'
    nv.zerofloorQ
end

starttime = tic;
disp('posterior sampling')

%get new points that are uniformly distributed in 5DOF space
pts3 = get_ocubo(npostpts); %posterior points



y3 = egprm();
X3 = proj_down(pts3,projtol,usv);
covmat = mdl.kfn(pts3,pts3); %assumes use of reference octonions from oreflist

%% Posterior Sampling
zerofloorQ = true;

ypost = tmvn(mu,covmat,l,u,n,method,zerofloorQ); %takes a long time for 10000^2 matrix
%     ypost = mvnrnd_trn(l.',u.',mu.',covmat,n);

runtime = toc(starttime);

end

%% CODE GRAVEYARD
%{
% l = zeros(size(mu));
% u = inf*ones(size(mu));


    l(:,1) double = max([zeros(size(mu)),0.895255*mu],[],2)-mu; %lower bound
    u(:,1) double = 1.2444*mu; %upper bound

egprmMdltmp = egprmMdl;
%}