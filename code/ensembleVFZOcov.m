function cov = ensembleVFZOcov(XN,XM,theta,usv,nv)
arguments
    XN(:,:) double
    XM(:,:) double
    theta(1,:) double
    usv = struct.empty
    nv.K(1,1) double = []
    nv.orefs = []
    nv.KernelFunction = 'squaredexponential'
end
% ENSEMBLEVFZODIST  compute a minimized GB covariance matrix across K VFZs based on reference octonions (orefs)

%unpack
K = nv.K;
orefs = nv.orefs;
KernelFunction = nv.KernelFunction;

if ~isempty(orefs) && ~isempty(K)
    error('only orefs or K should be specified')
end
if isempty(orefs) && isempty(K)
    K = 10;
end
if isempty(orefs)
    orefs = get_ocubo(K,'random',[],10);
end

%loop through VFZs
pd = [];
for i = 1:K
    oref = orefs(i,:);
    
    % symmetrize both sets of points
    o = proj_up(XN,usv);
    o = sqrt2norm(o,'quat');
    o = get_octpairs(o,'oref',oref);
    
    o2 = proj_up(XM,usv);
    o2 = sqrt2norm(o2,'quat');
    o2 = get_octpairs(o2,'oref',oref);
    
    %keep taking minimum distances
    pd = cat(3,pd,pdist2(o,o2));
    pd = min(pd,3);
end

% convert from distances to covariances
cov = d2cov(pd,theta,KernelFunction);

end