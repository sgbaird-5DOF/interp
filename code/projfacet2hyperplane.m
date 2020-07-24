function newvertices = projfacet2hyperplane(nvec,vertices)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: Take a list of vertices defining a facet with vertices on
% the unit hypersphere and project them onto a tangent hyperplane at a
% user-defined point, nvec. Useful for computing spherical barycentric
% coordinates (sphbary.m).
%
% Inputs:
%	nvec			===	normalized normal to desired tangent hyperplane (i.e.
%							point on unit hypersphere)
%
%	vertices		===	rows of vertices of facet to project from
%							hypersphere to tangent hyperplane @ nvec
%
% Outputs:
%	newvertices ===	new vertices of facet (projected from hypersphere to
%							tangent hyperplane @ nvec)
%
% Dependencies:
%	projray2hyperplane.m
%
% References:
%		https://math.stackexchange.com/q/1256236/798661
%--------------------------------------------------------------------------

%% project vertices onto hyperplane
%initialize
newvertices = zeros(size(vertices));

nvertices = size(vertices,1);

for i = 1:nvertices
	vtx = vertices(i,:); %vertex
	newvertices(i,:) = projray2hyperplane(nvec,vtx);
end

end %projfacet2hyperplane