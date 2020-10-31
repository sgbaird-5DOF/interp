function	meshpts = hypersphereSetup(n,d,subdivType)
arguments
    n
    d
    subdivType {mustBeMember(subdivType,{'orthant','hypersphere'})} = 'hypersphere'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Description:
%
% Inputs:
%		n				===	length of one side of d-dimensional lattice to
%								generate hypersphere points on. In current
%								implementation (2020-07-01), n needs to be odd in
%								order to have a symmetric sampling.
%
%		d				===	dimensionality
%
%		subdivType	===	type of subdivision to use (e.g. orthant or
%								hypersphere)
%
% Outputs:
%		meshpts		===	mesh points
%
% Dependencies:
%		hypersphere.m
%
% Note: n == 3, subdivType == 'hypersphere' gives an n-orthoplex or
% hyperoctrahedron, i.e. a simplex with 2^n facets where each orthant
% contains one facet. The default for subdivType == 'orthant' gives the
% non-negative orthant.
%--------------------------------------------------------------------------

if mod(n,2) == 0
	warning('n is even. Centering will be off. Sphere points will be asymmetric.')
end

disp(['creating hypersphere in ' int2str(d) 'D using sidelength of ' int2str(n) ' pixels'])

%get hypersphere lattice (hypersphere points == 1, 0 otherwise)
sz = repelem(n,d);
s = hypersphere(sz);

%initalize
mysub = cell(d,1);

%extract indices of hypersphere points
npts = sum(any(s),'all');
ind = find(s);

%extract indices in order to get un-normalized vectors
[ mysub{:} ] = ind2sub(size(s),ind); %named mysub (as in subscripts) because sub is a built-in function

%center the un-normalized vectors
subcenter = cellfun(@(sub) (sub-ceil(n/2))/(n/2),mysub,'UniformOutput',false);

smat = horzcat(subcenter{:}); %matrix form instead of cell form

switch subdivType
	case 'orthant'
		smat = smat(all(smat >= 0, 2),:); %subdivide the hypersphere into a hyperquadrant
	case 'hypersphere'
		%smat = smat;
end

%normalized vectors such that every point is on unit sphere (assumes centering at 0,0,0)
meshpts = smat./repelem(vecnorm(smat,2,2),1,d);
end


%------------------------CODE GRAVEYARD------------------------------------
%{
switch d
	case embedcases
		sz = [1,repelem(n,d)];
		matsize = repelem(n,d+1);
		s = hypersphere(sz,matsize);
		%initalize
		mysub = cell(d+1,1);
	otherwise


end


switch d
	case embedcases
		%normalized vectors such that every point is on unit sphere (assumes centering at 0,0,0)
		meshpts = smat./repelem(vecnorm(smat,2,2),1,d+1);

	otherwise

end
%}
