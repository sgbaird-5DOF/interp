function [newpts,TRI,varargout] = facet_subdiv(pts,nint,delaunayQ,convhullQ)
arguments
	pts {mustBeReal,mustBeFinite}
	nint(1,1) double = 1
	delaunayQ(1,1) logical = true
	convhullQ(1,1) logical = false
end
% FACET_SUBDIV  Project a facet from n-dimensional space to a simplex in n-1
% dimensions, subdivide the simplex, compute the triangulation, and
% reproject back to a facet in n-dimensional space
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-27
%
% Inputs:
%		pts			===	n-dimensional vertices of facet
%
%		nint			===	# of intervals/levels with which to subdivide
%								facet, nint == 1 does nothing, nint == 2 subdivides
%								once, nint == 3 subdivides twice, etc.
%
%		delaunayQ	===	whether or not to compute and output a delaunay
%								triangulation along with the subdivided points.
%
% Outputs:
%		newpts		===	subdivided points (including pts) projected back
%								into n-dimensional space (rows of points)
%
%		TRI			===	delaunay triangulation of facet in d-1 dimensional
%								space. Used to index newpts in n-dimensional space.
%
% Dependencies:
%		simplex_subdiv.m
%
% See also: 
%	SphereUniformSamplingToolbox (FEX)
%		TriQuad.m, SubdivideSphericalMesh.m, TriangleRayIntersection.m,
%		TriangleRayIntersection_tutorial.m
%
%	ManifoldToolbox in Gerardus Project
%		delaunay_sphere.m, proj_on_sphere.m, remove_vertex_from_tri.m,
%		closest_trifacet.m, surface_tridomain.m trihemisphere.m,
%		trifacet_area3D.m, plot_dmatrix.m
%--------------------------------------------------------------------------
d = size(pts,2);

[projpts,usv] = proj_down(pts,1e-6,struct.empty,'zeroQ',false);

%subdivide the facet turned into simplex
if nint > 1
	subdivpts = simplex_subdiv(projpts,nint);
else
	subdivpts = projpts;
end

nnew = size(subdivpts,1); %number of new points

%compute the delaunay triangulation
if delaunayQ
	TRI = delaunayn(subdivpts);
else
	TRI = 1:d;
end

%compute the convex hull
if convhullQ && (nnew > d-1)
	varargout{1} = convhulln(subdivpts);
else
	varargout{1} = 1:d-1;
end

newpts = proj_up(subdivpts,usv);

end


%------------------------------CODE GRAVEYARD------------------------------
%{
% if abs(S(end,end)) < tol
% 	ndegdim = 1;
% else
% 	ndegdim = 0;
% 	warning('ndegdim == 0')
% end




% if size(projpts,2) < size(pts,2) - 1
% 	disp('more than one degenerate dimension')
%  	newpts = [];
%  	TRI = [];
%  	return
% end

%project to d-1 dimensional space
% [U,S,V]=svd(pts-mean(pts),0);

% ndegdim = sum(abs(diag(S)) < tol); %number of degenerate dimensions
% 2020-07-22 changed to only go down a single dimension or none

% if ndegdim > 1
% % 	disp('more than one degenerate dimension')
%  	newpts = [];
%  	TRI = [];
%  	return
% end

%project to lower dimension
% projpts = U*S(:,1:d-ndegdim);

%find columns with constant values (i.e. degenerate)
% avgpts = mean(projpts);
% ids = find(all(projpts - avgpts < tol,1));

% if ~isempty(ids)
% 	splice_vals = avgpts(ids);
% 	projpts(:,ids) = []; %remove constant columns
% 	idstemp = ids + 1 - (1:length(ids));
% 	zeroQ = true;
% 	d = size(pts,2);
% 
% else
% 	zeroQ = false;
% end


% [a,usv] = proj_down([0 0 0; pts],1e-6);
% 
% zeropt = a(1,:);
% projpts = a(2:end,:)-zeropt;

% [projpts,usv] = proj_down(pts,tol,struct.empty,'zeroQ',true);


	% Create rotation matrix
% 	theta = 22.5; % to rotate 90 counterclockwise
% 	R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% 	TRI = delaunayn((R*subdivpts.').');


% 	t=n2c(subdivpts);
% 	triplot(TRI,t{:})



% if zeroQ
% 	%add columns of zeros back
% 	for i = 1:length(idstemp)
% 		id = idstemp(i);
% 		splice_cols = splice_vals(i)*ones(nnew,1);
% 		subdivpts = [subdivpts(:,1:id-1) splice_cols subdivpts(:,id:end)];
% 	end
% end
%project back into d-dimensional space
% newpts = padarray(subdivpts,[0 ndegdim],'post')*V'+mean(pts);

npts = size(pts,1);

if nargin == 4
	convhullQ = varargin{1};
else
	convhullQ = false;
end

%}

