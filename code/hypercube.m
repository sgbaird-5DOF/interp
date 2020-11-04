function [pts,K] = hypercube(d)
% HYPERCUBE  Calculate vertices and convex hull triangulation of a hypercube
%  in d-dimensions centered at origin and edge length of 1.
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-06
% 
% Inputs:
%		d		=== dimensionality
%
% Outputs:
%		pts	=== vertices of the hypercube
%
%		K		=== convex hull triangulation of the hypercube
%
% Dependencies:
%		allcomb.m (available on FEX)
%
%		normr.m	(built-in via Computer Vision toolbox, or use an alternative
%					I wrote if you don't have the Computer Vision toolbox)
%
%--------------------------------------------------------------------------
%vertices and triangulation
sidelength = 2; %of hypercube centered at origin
t = {sidelength/2*[-1 1]};
tlist = repelem(t,d);

normQ = true;

pts = allcomb(tlist{:});
if normQ
	pts = normr(pts);
end
K = convhulln(pts);

end