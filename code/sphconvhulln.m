function K = sphconvhulln(pts,maxnormQ)
arguments
	pts double {mustBeFinite,mustBeReal}
	maxnormQ(1,1) logical = true
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-06
%
% Description: Compute the convex hull of a section of a (or a full)
% hypersphere (note [1]) without producing "undercut" facets (note
% [3]). Hypercube points are placed in the orthants (defined by an
% orthoplex) that don't contain datapoints and are added to the input to
% convhulln(). After computing the triangulation, facets connected to the
% hypercube points are removed (note [2]) and the numbering is made
% consistent with the original input.
%
% Example usage: K = sphconvhulln(pts); K = sphconvhulln(pts,maxnormQ);
% 
% Inputs:
%		pts	===	rows of cospherical points that define a hypersphere
%						surface.
%
%		maxnormQ ===	whether to use intersect_facet with maxnormQ == true.
%							Could take a while if too many points/too high of
%							dimension.
%
% Outputs:
%		K		===	triangulation of pts without "undercut" facets.
%
% Dependencies:
%		hypercube.m
%			-allcomb.m
%			-normr.m
%
%		orthoplex.m
%
%		intersect_facet.m
%			-projray2hypersphere.m
%				--numStabBary.m (optional)
%
%		normr.m
%
% Notes:
%		[1] Assumes that all input points (pts) already define a convex hull. Not
%		tested without this assumption, but an easy solution is to parse the
%		input with a call to convhulln() before passing it into
%		sphconvhulln().
%
%		[2] If hypercube or orthoplex points are already in the set of pts, then
%		facets connected to these are kept.
%
%		[3] An example of an undercut facet would be the largest face on a
%		hemisphere (i.e. the planar face that would "close" the hemisphere).
%		See sphconvhulln_test.m with rmUndercutQ toggled to false and true in
%		separate runs.
%
%		I could deal with points lying on axis plane intersections by
%		rotating the orthoplex slightly if such points exist.
%
%		Could check to see if points fall on a hemisphere by computing the
%		regular convex hull, taking a few facet averages, and checking to see
%		if any of them have two averages. Might not be infallible, but would
%		probably do the trick. There's also an orthant check (make sure that
%		no datapoint has a corresponding datapoint in the opposite orthant,
%		and that less than 2^n orthants have intersections. Might not work if
%		it's exactly a hemisphere.
%
%		maxnormQ == false might not be working correctly (2020-07-29)
%
%--------------------------------------------------------------------------
%% setup
d = size(pts,2); %dimension

% helper functions
prec = 12;
r = @(a) round(a,prec);
tol = 1e-6;
mymembercheck = @(a,b) ismembertol(r(a),r(b),tol,'ByRows',true);

% check complexity
if d >= 7 && size(pts,1) > 400
	warning(['d = ' int2str(d) ' and npts = ' ...
		int2str(size(pts,1)) '. Initial convex hull input may be too big.'])
	slurmQ = true;
	if ~slurmQ
		m = input('Continue? y/n:','s');
		if ~strcmp(m,'y') && ~strcmp(m,'Y')
			return
		end
	end
end

if maxnormQ
	%check to see if data falls on less than hemisphere
% 	if size(pts,2) < 5
% 		opts = {'QJ'};
% 	else
% 		opts = {'Qx','QJ'};
% 	end
% 	K = convhulln(pts,opts);

	K = convhulln(pts);
	
	subhemiQ = false;
	k = 0;
	while ~subhemiQ && (k < 10)
		rid = randi(size(K,1));
		avg1 = mean(pts(K(rid,:),:));
% 		intTemp = intersect_facet(pts,K,avg1,tol);
		[~,~,~,intTemp,~] = projray2hypersphere(pts,K,avg1,1e-6,false);
		subhemiQ = length(intTemp) > 1;
		k = k+1;
	end
	
	if subhemiQ
		
		projpts = projfacet2hyperplane(normr(mean(pts)),pts);
		
		a = proj_down(projpts,1e-6);
		
		K = delaunayn(a);
		
		%compute convex hull with extra point
		%extrapt = 0.1*normr(mean(pts)); %assumes points fall on less than a hemisphere
		%K = convhulln([pts;extrapt]);
		
		%delete everything connected to extra point
		%npts = size(pts,1);
		%[row,~] = find(K == npts+1);
		%K(row,:) = [];
		
	else
		K = convhulln(pts); %compute regular convex hull
	end
	
	return
end

if maxnormQ == false
	disp('maxnormQ == false might not be working correctly (2020-07-29)')
end

%% add orthoplex points (if not in pts) to prevent "undercut facets"

% pts = uniquetol(pts,1e-6,'ByRows',true);

%normalize the points to be on a unit sphere
if norm(pts(1,:)) ~= 0
	pts = normr(pts);
end

% pts(abs(pts) < 1e-12) = 0;

%define hypercube vertices & triangulation, 1 point per orthant, no vertices on major axes
[hcubePts,~] = hypercube(d);

%define orthoplex vertices & triangulation, 1 point on each major axis
[orthoPts,orthoK] = orthoplex(d);

% don't remove hypercube points that are already in pts (i.e. a check for
% duplicates, and don't remove duplicate pts & facets connected to them)
nonDupPtIDs = find(~mymembercheck(hcubePts,pts));
nonDupPts = hcubePts(nonDupPtIDs,:);

% don't remove orthoPts (& facets) that are a part of pts
nonDupPtIDsOrtho = find(~mymembercheck(orthoPts,pts));
nonDupPtsOrtho = orthoPts(nonDupPtIDsOrtho,:);
	
intcheckPts = nonDupPtsOrtho; % intersection check points

% find orthoplex points that don't intersect with traditional convex hull
if size(pts,1) <= d
	K = 1:d;
else
	K = convhulln(pts);
end

temp = intersect_facet(pts,K,intcheckPts,tol); %facets formed by pts that nonDupPtsOrtho intersect 
rmintcheckIDs = find(cellfun(@isempty,temp)); % nonDupPtsOrtho points that do not intersect K
% rmOrthoIDs = setdiff(rmOrthoIDs,nonDupPtIDsOrtho);
rmintcheckPts = intcheckPts(rmintcheckIDs,:);


%% find orthants that hcubePts and pts lie in, respectively
hcubeintIDs = intersect_facet(r(orthoPts),orthoK,r(hcubePts),tol);
hcubeintIDs = vertcat(hcubeintIDs{:});

intfacetIDs = intersect_facet(r(orthoPts),orthoK,pts,tol);
intfacetIDs = vertcat(intfacetIDs{:});

orthants = unique(intfacetIDs);

%% IDs of hcubePts to be removed from convex hull
%	along with facets that contain those points i.e. find hcubePts that have
%	a datapoint (from pts) in the same orthant. Don't include these in the
%	convexhull.
rmIDs = find(ismember(hcubeintIDs,orthants));

% make sure that you don't delete one of the "keep" points
rmIDs = setdiff(nonDupPtIDs,rmIDs);

%points to remove (along with facets connected to these) after triangulation
rmpts = [hcubePts(rmIDs,:); rmintcheckPts];

%% remove vertices in orthants
nptsrm = size(rmpts,1); %number of remove points

%catenated points (to prevent undercuts during convhulln() )
catpts = [rmpts;pts];

%% Compute convex hull
K = convhulln(catpts); %(without undercuts)

%% remove rmpts and facets connected to rmpts
[row,col] = find(ismember(K,1:nptsrm)); %facet rows to remove, with duplicates
row = unique(row);
K(row,:) = []; %delete facets
K = K - nptsrm; %change back to original numbering

if isempty(K)
	warning('Convex hull is empty. Too many facets removed.')
end

end

%%
%-----------------------------CODE GRAVEYARD-------------------------------
%{
%create orthoplex
K = convhull(orthpts);

% npts = length(pts);
% orthpts = [eye(d);-eye(d)]; % define orthoplex points

% hcubeorth = unique(hcubeintIDs);

%rmIDs = npts+1:npts+nptsrm-1; %IDs of points to remove referenced to catpts

catpts(rmIDs,:) = []; %delete vertices

[row,col] = find(ismember(K,rmIDs)); %facet rows to remove, with duplicates



ptscheck = pts;
ptscheck(nonDupPtIDsOrtho,:) = [];


% remove ortho pts that 
orthoPts(nonDupPtIDsOrtho,:) = [];
[row,~] = find(ismember(orthoK,nonDupPtIDsOrtho));
orthoK(row,:) = [];

orthants1 = dupPtIDsOrtho;

%% remove facets that have a projection magnitude of zero
%i.e. undercut hemisphere facets
[intfacetIDs,dataProj] = intersect_facet(pts,K,orthoPts);
intfacetIDs = vertcat(intfacetIDs{:});


hcubeorthoPts = vertcat(hcubePts,orthoPts);
hcubeorthoK = convhulln(hcubeorthoPts);
ptstempIDs = find(ismembertol(pts,orthoPts,1e-6,'ByRows',true));
intfacetIDs2 = intersect_facet(hcubePts,hcubeK,pts);

% keep track of orthoPts that are a part of pts
dupPtIDsOrtho = find(ismembertol(pts,orthoPts,1e-6,'ByRows',true));
% dupPtsOrtho = orthoPts(dupPtIDsOrtho,:);

for i = 1:length(intfacetIDs)
	tempIDs = intfacetIDs{i};
	dupPtIDsOrtho
end


%if all of the intersecting facets happen to lie on major axes, and 6 is
%the only one that doesn't have an intersection with 4 repeats..

% %ignore pts that lie on an axis or a plane defined by multiple axes (i.e.
% %vertices with at least one zero component, but not all zeros)
% intcheckIDs = find(ismembertol(round(pts,12),round(axpPts,12),1e-6,'ByRows',true));

intfacetIDs = intersect_facet(r(orthoPts),orthoK,pts,tol);

% %only consider pts that are axpPts as being in single orthant
% for i = 1:length(intcheckIDs)
% 	intcheckID = intcheckIDs(i);
% 	intfacetIDs{intcheckID} = intfacetIDs{intcheckID}(1);
% end


if ~isempty(find(mymembercheck(axpPts,pts),1));
	% don't remove axpPts (& facets) that are a part of pts
	nonDupPtIDsAx = find(~mymembercheck(axpPts,pts));
	nonDupPtsAx = axpPts(nonDupPtIDsAx,:);
else
	nonDupPtsAx = [];
end

% 	extrapt = zeros(1,d);



% 		%calculate average of facet vertices for each facet
% 		nfacets = size(K,1);
% 		avgpts = zeros(nfacets,d);
% 		for i = 1:nfacets
% 			facetptIDs = K(i,:);
% 			avgpts(i,:) = mean(pts(facetptIDs,:));
% 		end
% 		
% 		%project facet averages onto each facet
% 		%take intersection with largest norm
% 		intfacetIDs = intersect_facet(pts,K,avgpts,tol,maxnormQ);
% 		intfacetIDs = vertcat(intfacetIDs{:});
% 		K = K(intfacetIDs,:);


%define axis-based polytope
% axpPts = axpolytope(d);
	
% check for duplicates between axis-based polytope points and pts
% check = ~isempty(find(mymembercheck(axpPts,pts),1));
% k = 0;
% while check && k < 10
% 	k = k+1;
% 	% don't remove axpPts (& facets) that are a part of pts
% 	R = slightRot(d);
% 	hcubePts = (R*hcubePts.').';
% 	axpPts = (R*axpPts.').';
% 	check = ~isempty(find(mymembercheck(axpPts,pts),1));
% end

%}