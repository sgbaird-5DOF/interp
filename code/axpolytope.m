function pts = axpolytope(d,varargin)
% AXPOLYTOPE Indices of a polytope with vertices on every axis and every intersection of every plane formed by axes.
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-17
% 
% Inputs:
%		d			===	dimensionality
%
%		intlist	===	list of numbers (max range from 0 to d) for which
%							points will be formed that have a number of axis plane
%							intersections that is a member of the set given in
%							intlist. Default is 1:d (2020-07-17)
%
% Outputs:
%		pts	=== vertices of axis based polytope
%
% Dependencies:
%		allcomb.m
%
%		normr.m
%
% Notes:
%		There's probably a precise geometric term for this. No convex hull
%		for this one since the number of points gets very large in high
%		dimensions.
%
%--------------------------------------------------------------------------
usual = 1;
if nargin - 1 == 1
	intlist = varargin{1};
else
	intlist = 1:d;
end
% all vertices will contain some combination of at least one element of
% base
base = {[-1 0 1]};

% make a cell array containing d-copies of base, one per cell
baselist = repelem(base,d);

% create initial point list
pts=allcomb(baselist{:});

%remove certain points based on user-specification
[row,~] = find(~ismember(sum(abs(pts),2),intlist));

pts(row,:) = [];

%normalize points to fall on unit sphere
pts = normr(pts);

end

%----------------------------CODE GRAVEYARD--------------------------------
%{
% pts = setdiff(pts,[zeros(1,d); ones(1,d)],'rows')

% % compute convex hull
% K = convhulln(pts);

%}