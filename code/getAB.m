function [A,B] = getAB(pts)
arguments
    pts
end
% GETAB  get octonion points that are far from each other within ppts and return the corresponding pts
%% pairwise distance
npts = size(pts,1);
if npts > 20000
    ids = 1:20000;
else
    ids = 1:npts;
end

%     pts = mdls{1}.mesh.pts;
pd = squareform(pdist(pts(ids,:)));
[mx,id] = max(pd,[],'all','linear');

% max dimension
indlen = length(ids);
[row,col] = ind2sub([indlen,indlen],id);
A = pts(row,:);
B = pts(col,:);
end