%spherical convhulln test
clear; close all
d = 3;

seed = 10;
rng(seed);

load_type = 'uppHem'; %'misFZ', 'uppHem','nonnegOrth','octant', '2xoctant','sphCone', 'axPlane', 'axPlane2'

switch load_type
	case 'misFZ'
		%load BP FZ
		load('misFZfeatures.mat');
		%'interior' - hemisphere,'AC' - two adjacent octants,'B' - octant,
		%'OA' - half an octant, 'C' - 1/3 octant, 'A' - 1/4 octant
		feature = 'A'; %'interior','AC','B','OA','C','A', etc.
		q = qlist.(feature);
		nint = 2;
		ctrcuspQ = false;
		[pts,A,R,TRI,af] = meshBP(q,nint,ctrcuspQ);
		pts = vertcat(pts.sub);
		
	case 'uppHem'
		%load random points in upper hemisphere
		pts = rand(100,d)-0.5;
		pts = pts(pts(:,1) > 0,:);
		
	case 'uppHem+'
		pts = rand(100,d)-0.5;
		tol = -1e-7; % once it falls below -tol2 from projray2hypersphere.m, the hemisphere+ closes
		pts = pts(pts(:,1) >= tol,:);
		pts(35:38,1) = tol;
		
	case 'sphCone'
		%load random points in upper hemisphere and add a point for the cone
		pts = rand(100,d)-0.5;
		pts = pts(pts(:,1) >= 0,:);
		pts = [pts;-1 repelem(0,d-1)];
		
	case 'nonnegOrth'
		%load random points in non-negative orthant
		pts = rand(100,d);
		
	case 'octant'
		%load octant
		pts = eye(d);
		
	case '2xoctant'
		pts = [eye(d);-eye(d)];
		keepPtIDs = (pts(:,1) >= 0) & (pts(:,2) >= 0);
		pts = pts(keepPtIDs,:);
		
	case 'axPlane'
		pts = [eye(d)-1; -eye(d)+1];
		
	case 'axPlane2'
		pts = [...
			1 1 0;
			1 0 1];
		for i = 1:3
			pts = [pts; pts(1,:)+0.5*rand(1,size(pts,2))];
		end
		
	case 'axPoly'
		pts = axpolytope(3);
		
end

pts = normr(pts);

rmUndercutsQ = true;

maxnormQ = true;

if rmUndercutsQ
	K = sphconvhulln(pts,maxnormQ);
else
	K = convhulln(pts);
end

%plotting
if d == 3
	ptscell = num2cell(pts,1);
	trisurf(K,ptscell{:})
	axis equal
	ax = gca;
	ax.View = [-45, -45];
	%ax.View = [25,10];
	[x,y,z] = sphere(40);
	scl = 0.5;
	x = scl*x;
	y = scl*y;
	z = scl*z;
	hold on
	ax2 = surf(x,y,z,'EdgeColor','none');
	%alpha(ax2,0.75);
end

axis equal tight vis3d off

% camlight


%Yes! 2020-07-04, 7:45 PM