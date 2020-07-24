function [c,Pnew] = sphbary(data,vertices,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: compute spherical barycentric coordinates of facet or simplex [1]
% 
% Inputs:
%		data		===	datapoint of interest on hypersphere which defines the
%							tangent hyperplane to the hypersphere
%
%		vertices ===	rows of vertices of facet that have been projected to
%							tangent hyperplane defined by 'a'
%
%		Pnew		===	varargin{1}, vertices projected onto tangent
%							hyperplane
%
% Outputs:
%		c			=== spherical barycentric coordinates of 'a'
%
% Dependencies:
%		projfacet2hyperplane.m
%			-projray2hyperplane.m
%
% References:
%		[1] T. Langer, A. Belyaev, H.-P. Seidel, Spherical barycentric
%		coordinates, Proc. Fourth Eurographics Symp. Geom. Process. (2006)
%		81–88. http://portal.acm.org/citation.cfm?id=1281957.1281968.
%
%		https://math.stackexchange.com/q/1256236/798661
%
%--------------------------------------------------------------------------
if nargin == 3
	Pnew = varargin{1};
end

if exist('Pnew','var') == 0
	Pnew = projfacet2hyperplane(data,vertices);
end

c = (Pnew.'\data.').';

end


%-----------------------------NOTES----------------------------------------
%{
Maybe sphbary should have the computation of a and P (from
projfacet2hyperplane.m) inside of it. 2020-07-09, now it does
%}