function [pts,K] = orthoplex(d)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-06
%
% Description: Find indices and convex hull triangulation of an orthoplex
% (i.e. n-dimensional analogue of octahedron).
% 
% Inputs:
%		d		=== dimensionality
%
% Outputs:
%		pts	=== vertices of orthoplex
%
%		K		=== triangulation of orthoplex
%
% Dependencies:
%
%--------------------------------------------------------------------------
%vertices and triangulation
pts = [eye(d);-eye(d)];
K = convhulln(pts);

end