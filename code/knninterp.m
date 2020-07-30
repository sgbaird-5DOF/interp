function interpvals = knninterp(meshpts,meshprops,datapts,K)
arguments
	meshpts double {mustBeReal,mustBeFinite}
	meshprops(:,1) double {mustBeReal,mustBeFinite}
	datapts double {mustBeReal,mustBeFinite}
	K(1,1) double {mustBeInteger} = size(meshpts,2)+50;
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: compute a linear interpolation of data using a hyperplane
% fitted to k-nearest neighbors of mesh (points and property values).
%
% Inputs:
%		a -	a
%
% Outputs:
%		b -	b
%
% Usage:
%		c = b(a);
%
% Dependencies:
%		*
%
% Notes:
%	griddatan could compute interpolated values based on a delaunay
%	triangulation; however, this is infeasible if the dimensionality and/or
%	number of points is too high. Thus, it is broken up into a nearest
%	neighbor search and then defining a triangulation & interpolation
%	locally.
%
% See also:
%  griddatan
%  knnsearch
%--------------------------------------------------------------------------

ndatapts = size(datapts,1);

%find k-nearest neighbors
IDs = knnsearch(meshpts,datapts,'K',K,'IncludeTies',true);

%initialize
interpvals = NaN(ndatapts,1);

%get interpolated values
for i = 1:ndatapts %parfor compatible
	%unpack id list
	ID = IDs{i};
	
	%griddatan setup
	x = meshpts(ID,:);
	v = meshprops(ID,:);
	xq = datapts(i,:);
	
	%project meshpts onto tangent hyperplane
	x = projfacet2hyperplane(xq,x);
	
	a = proj_down([x;xq],1e-3);
	x = a(1:end-1,:);
	xq = a(end,:);
	
	%interpolation
	interpvals(i) = griddatan(x,v,xq);
	
	if isnan(interpvals(i))
		warning(['no intersection found for i == ' int2str(i) '. Taking nearest neighbor interpolation value.'])
		interpvals(i) = griddatan(x,v,xq,'nearest');
	end
	
	if i == 7
		figure
		triplot(delaunayn(x),x(:,1),x(:,2))
		hold on
		plot(xq(1),xq(2),'*')
		1+1;
	end
	
end

end %knninterp