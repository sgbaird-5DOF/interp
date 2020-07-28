function [newpts,A,R,TRI,af,sphpts] = meshBP(q,nint,ctrcuspQ,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-30
%
% Description: generate boundary plane (BP) fundamental zone mesh for a
% quaternion and rotate quaternion so that the BP normal matches the
% quaternion.
%
% Inputs:
%		q			===	quaternion (according to convention in [1])
%
%		nint		===	number of intervals (i.e. subdivisions),
%							where n == 1 does zero subdivisions, n == 2 does one
%							subdivision, etc.
%		ctrcuspQ ===	whether or not to center cusps (true) or have them at
%							the boundaries (false)
%
% Outputs:
%		newpts	===	newsphpts, but in Cartesian coordinates instead of
%							spherical. Rows of points
%
%		A			===	% symmetry axes, [ax,ay,az], column vectors, to
%							unalign or align the BP normal with the symmetry axes.
%
%		TRI		===	triangulation of newsphpts. See delaunay() or
%							delaunayn() for more details.
%
%		af			===	final arc position (of first quadrant)
%
%	Note:
%		newpts and TRI currently (2020-07-01) implemented in struct() format
%		with main and sub corresponding to main spherical triangle and
%		corresponding subdivision.
%
% Dependencies:
%
%		findgeometry.m
%
%		symaxis.m
%
%		sphtri_subdiv.m
%			-facet_subdiv.m
%				--simplex_subdiv.m
%
%		orthoplex.m
%
% References:
%		[1] S. Patala, C.A. Schuh, Symmetries in the representation of grain
%		boundary-plane distributions, Philos. Mag. 93 (2013) 524–573.
%		https://doi.org/10.1080/14786435.2012.722700.
%--------------------------------------------------------------------------

if nargin == 4
	geometry = varargin{1};
else
	geometry = findgeometry(q); %geometry in misorientation FZ, such as 'line OB'
	geometry = geometry{1};
end
% geometry = 'twosphere'; %just for testing
[A,R] = symaxis(q,geometry);

%% assign arclength & make coarse grid
switch geometry
	%consider adding these values to 'misFZfeatures.mat' via plotFZrodriguez_test.m
	case {'OAB','OBCE','OADE','CDE'}
		%surfaces
		ptgrp = 'C2h';
		alen = pi;
	case {'OB','CE','ED'}
		%lines
		ptgrp = 'D2h';
		alen = pi/2;
	case 'OE'
		ptgrp = 'D3d';
		alen = pi/6;
	case 'OA'
		ptgrp = 'D4h';
		alen = pi/4;
	case 'AC'
		ptgrp = 'C2h';
		alen = pi;
		
	case {'B','D'}
		%points
		ptgrp = 'D2h';
		alen = pi/2;
	case 'E'
		ptgrp = 'D6h';
		alen = pi/6;
	case 'A'
		ptgrp = 'D8h';
		alen = pi/8; %maybe not be correct: I thought D8h wasn't possible for a lattice, but maybe OK since it's a GB
	case 'C'
		ptgrp = 'D4h';
		alen = pi/4;
	case 'O'
		ptgrp = 'Oh';
		alen = pi/4; %not accurate (more area than a true BP FZ), need to adjust elevation angle
	case 'interior'
		ptgrp = 'C1';
		alen = 2*pi;
	case 'twosphere'
		alen = 2*pi+1e-6;
end

if ~ctrcuspQ
	a0 = 0; %beginning of arc
else
	a0 = alen/2; %center starting point halfway through region (so that mirror planes, etc. are encompassed)
end

% construct top-level spherical triangles
if strcmp(geometry,'O')
		af = a0+alen; %end of arc
	sphpts{1} = [...
		a0		  	pi/2				%pt1
		a0			pi/4				%pt2
		af			atan(sqrt(2))	%pt3
		]; %[az; el]
	
elseif (alen <= pi/2) && (alen > 0)
	%single octant
	af = a0+alen; %end of arc
	sphpts{1} = [...
		a0		  	pi/2	%pt1
		a0			0  	%pt2
		af			0		%pt3
		]; %[az; el]
	
elseif (alen > pi/2) && (alen <= pi)
	%two octants
	af = a0+alen/2;
	sphpts{1} = [...
		a0			 	pi/2	%pt1
		a0				0  	%pt2
		af				0		%pt3
		]; %[az; el]
	sphpts{2} = sphpts{1};
	sphpts{2}(:,1) = sphpts{1}(:,1) + alen/2*repelem(1,3,1);
	
elseif (alen > pi) && (alen <= 2*pi)
	%four octants (hemisphere)
	af = a0+alen/4;
	sphpts{1} = [...
		a0				pi/2	%pt1
		a0				0  	%pt2
		af				0		%pt3
		]; %[az; el]
	for i = 2:4
		sphpts{i} = sphpts{i-1};
		sphpts{i}(:,1) = sphpts{i-1}(:,1) + alen/4*repelem(1,3,1);
	end
	
elseif alen > 2*pi 
	%create full sphere
	af = [];
	[orthopts,orthoK] = orthoplex(3);
	tmp = num2cell(orthopts,1);
	[az,el,~] = cart2sph(tmp{:});
	sphtemp = [az el];
	for i = 1:size(orthoK,1)
		vertices = orthoK(i,:);
		sphpts{i} = sphtemp(vertices,:);
	end
end

%% subdivide spherical triangles (for each quadrant if necessary)
n = length(sphpts);

%initialize struct (consider switching from .main and .sub struct notation
%to cell(1,2), where 1 === main, 2 === sub for multi-level subdivisions)
pts(n) = struct();
TRI(n) = struct();
newpts(n) = struct();

pts(1).main = [];
TRI(1).main = [];
newpts(1).main = [];

if n == 1
	pts(1).sub = [];
	TRI(1).sub = [];
	newpts(1).sub = [];
else
	pts(1).sub = cell(1,n);
	TRI(1).sub = cell(1,n);
	newpts(1).sub = cell(1,n);
end

%subdivisions and rotations

for i = 1:n
	%get main spherical triangle and subdivisions (if applicable)
	[~,pts(i).main,TRI(i).main] = sphtri_subdiv(sphpts{i},1);
	[~,pts(i).sub,TRI(i).sub] = sphtri_subdiv(sphpts{i},nint);
	
	%rotate GB normal points to match quaternion
	newpts(i).main = (R*pts(i).main.').';
	newpts(i).sub = (R*pts(i).sub.').';
end

end


%{
%---------------------------CODE GRAVEYARD---------------------------------
%recompute triangulation (note, delaunay is already called deeper inside
%sphtri_subdiv, so calling it again isn't ideal)
TRI = convhull(pts); % since it's not a full sphere, you get spurious facets


if n == 1
	[~,pts(1).main,TRI(1).main] = sphtri_subdiv(sphpts{1},1);
	[~,pts(1).sub,TRI(1).sub] = sphtri_subdiv(sphpts{1},nint);
	newpts(1).main = (R*pts(1).main.').';
	
else




% construct top-level spherical triangles
if (alen <= pi/2) && (alen > 0)
	%single octant
	af = a0+alen; %end of arc
	sphpts{1} = [...
	a0		  	pi/2	%pt1
	a0			0  	%pt2
	af			0		%pt3
	]; %[az; el]
	
elseif (alen > pi/2) && (alen <= pi)
	%two octants
	af = a0+alen/2;
	sphpts{1} = [...
		a0			 	pi/2	%pt1
		a0				0  	%pt2
		af				0		%pt3
		]; %[az; el]
	sphpts{2} = sphpts{1}(3,:) + [alen/2,0];
	
elseif (alen > pi) && (alen <= 2*pi)
	%four octants (hemisphere)
	af = a0+alen/4;
	sphpts{1} = [...
		a0				pi/2	%pt1
		a0				0  	%pt2
		af				0		%pt3
		]; %[az; el]
	for i = 2:4
		sphpts{i} = sphpts{1}(3,:) + i*[alen/4,0];
	end
end

%% subdivide spherical triangles (for each quadrant if necessary)
n = length(sphpts);

%initialize struct (consider switching from .main and .sub struct notation
%to cell(1,2), where 1 === main, 2 === sub for multi-level subdivisions)
pts(n) = struct();
TRI(n) = struct();
newpts(n) = struct();

pts(1).main = [];
TRI(1).main = [];
newpts(1).main = [];

if n == 1
	pts(1).sub = [];
	TRI(1).sub = [];
	newpts(1).sub = [];
else
	pts(1).sub = cell(1,n);
	TRI(1).sub = cell(1,n);
	newpts(1).sub = cell(1,n);
end

%subdivisions and rotations

catsphpts = vertcat(sphpts{:});
npts = size(catsphpts,1);

catsphptsCell = num2cell([catsphpts repelem(1,npts,1)],1);
[x,y,z] = sph2cart(catsphptsCell{:});

K = sphconvhulln([x y z]);

nfacets = size(K,1);

for i = 1:nfacets
	vtxIDs = K(i,:);
	ptsTemp = catsphpts(vtxIDs,:);
	%get main spherical triangle and subdivisions (if applicable)
	[~,pts(i).main,TRI(i).main] = sphtri_subdiv(ptsTemp,1);
	[~,pts(i).sub,TRI(i).sub] = sphtri_subdiv(ptsTemp,nint);
	
	%rotate GB normal points to match quaternion
	newpts(i).main = (R*pts(i).main.').';
	newpts(i).sub = (R*pts(i).sub.').';
end



%subdivisions and rotations
nptstot1 = 0;
nptstot2 = 0;
for i = 1:n
	%get main spherical triangle and subdivisions (if applicable)
	[~,pts(i).main,TRI(i).main] = sphtri_subdiv(sphpts{i},1);
	TRI(i).main = TRI(i).main + nptstot1;
	nptstot1 = nptstot1 + size(pts(i).main,1);
	
	[~,pts(i).sub,TRI(i).sub] = sphtri_subdiv(sphpts{i},nint);
	TRI(i).sub = TRI(i).sub + nptstot2;
	nptstot2 = nptstot2 + size(pts(i).sub,1);
	
	%rotate GB normal points to match quaternion
	newpts(i).main = (R*pts(i).main.').';
	newpts(i).sub = (R*pts(i).sub.').';
end


	% full sphere
	% 	af = a0+alen/4;
	% 	sphpts{1} = [...
	% 		a0				pi/2	%pt1
	% 		a0				0  	%pt2
	% 		af				0		%pt3
	% 		]; %[az; el]
	% 	for i = 2:4
	% 		sphpts{i} = sphpts{i-1};
	% 		sphpts{i}(:,1) = sphpts{i-1}(:,1) + alen/4*repelem(1,3,1);
	% 	end
	%
	% 	for i = 5:8
	% 		sphpts{i} = -sphpts{i-4};
	% 	end

[alen,af,sphpts] = getBPFZ(geometry,ctrcuspQ)
%}
