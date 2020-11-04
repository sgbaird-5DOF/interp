function [newsphpts,newpts,TRI] = sphtri_subdiv(sphpts,nint)
% Subdivision scheme that segments the triangle formed by vertices of a spherical triangle
%  by partitioning the edges of the triangle, then projects those points
%  back onto the sphere. This is not an equal area approach (2020-06-27).
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-27
%
% Input:
%		sphpts		===	points on the unit sphere in spherical coordinates,
%		e.g. sphpts = [0 pi/2; 0 0; pi/2 0]; % [azimuthal elevation]. Maximum
%		arclength between any two points must be < pi. Anisotropies in the
%		arclengths will be propagated by the subdivision technique.
%
%		nint			===	number of intervals with which to subdivide the
%		spherical triangle. nint == 1 does nothing, nint == 2 creates 4 new
%		spherical triangles, nint == 2 creates 9 new spherical triangles.
%		Note: no matter what nint is, it will only perform one level of
%		subdivision (i.e. a single subdivision step, and a single projection
%		of points back onto the sphere).
%
% Output:
%
%		newsphpts	===	new points on the unit sphere after subdivision, in
%		same format as sphpts: [azimuthal elevation]
%
%		newpts		===	Cartesian coordinates of newsphpts (r == 1), in
%		[x,y,z] format, with rows of points.
%
%		TRI			===	triangulation of newsphpts. See delaunay() or
%		delaunayn() for more details.
%
%	Dependencies:
%		facet_subdiv.m
%			simplex_subdiv.m
%--------------------------------------------------------------------------

npts = size(sphpts,1);

%convert to cartesian
azimuth = sphpts(:,1);
elevation = sphpts(:,2);
r = repelem(1,npts,1);

[x,y,z] = sph2cart(azimuth,elevation,r);
pts = [x,y,z]; %rows of pts

%subdivide
delaunayQ = true;
[newpts,TRI] = facet_subdiv(pts,nint,delaunayQ);

%project back to unit sphere
newpts = newpts./vecnorm(newpts,2,2);

%convert back to spherical
xnew = newpts(:,1);
ynew = newpts(:,2);
znew = newpts(:,3);

[aznew,elnew,rnew] = cart2sph(xnew,ynew,znew);

newsphpts = [aznew,elnew];
end

%{
-------------------------------CODE GRAVEYARD------------------------------
%P = eye(3); %rows of points defining octant which is a spherical triangle
%[azimuth, elevation, r] = cart2sph(P(:,1),P(:,2),P(:,3));

azimuth   =	[0			0		pi/2];
elevation = [0			pi/2	0	 ];
r			 = [1			1		1	 ];

%}


