clear; close all force
%---------------------------------run.m------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Description:
%
% References:
%		Check if ray intersects internals of D-facet:
%		https://math.stackexchange.com/q/1256236
%--------------------------------------------------------------------------

%{
Generate an approximated sphere using gridded matrix points and
hypersphere(), normalize the vectors, and compute a convex hull for
triangulation and interpolation. Barycentric coordinates can be output as
well, and these can be extended into spherical barycentric coordinates by
projecting all vertices onto a hyperplane tangent to the point of interest
and using the barycentric coordinates of that hyperplane.

n: odd values preferred to reduce chance of single points in separate plane
which could mess up centering. I think centering formula also needs to
change if using even numbers. n = 15 gives "out of memory" on 4 GB RAM
surface pro 4 i5 processor using 7D sphere embedded in 8D. n = 13 for same
gives ~0.9 GB matrix for s, and an extra few hundred MB for other variables
combined. 1.6E6 points for n = 13 in 7Dsphere@8D case.

convhulln() for a hypersphere with n = 5 and d = 7 (d === dimension), takes
~4 min (2982 vertices,1082756 facets) n = 7, d = 7 gives 30578 vertices,
and 2.8591e13 facets, or 2.64e7 x more facets Given 1e4 vertices in 7D, we
get 1e12 facets, and might take ~7 years to compute..

Once the dimension & # vertices gets high enough (e.g. 10,000 pts in 7D
cartesian, it would probably make more sense to use a k-nearest neighbor
approach and fit a hyperplane using griddatan. Or, I could do a subdivision
scheme once I've found the containing facet. This way, I can triangulate it
in a sort of r-tree indexing scheme. Start out with a coarse mesh and
refine it.
%}

tic
%% Setup

%set random number generator
seed = 10;
rng(seed);

%mesh and data types

addpathdir({'misFZfeatures.mat','PGnames.mat','nlt.m','q2rod.m',...
	'GBfive2oct.m','qu2ax.m'})

%'hypersphereGrid', 'Rohrer2009', 'Kim2011', '5DOF',
%'Olmsted2004','5DOF_vtx','5DOF_misFZfeatures',
%'5DOF_interior','5DOF_exterior', '5DOF_oct_vtx','5DOF_hsphext'
%'5DOF_exterior_hsphext'
meshMethod = '5DOF_exterior';
dataMethod = '5DOF_interior';
pseudoMethod = [meshMethod '_pseudo'];

%initialize
meshopts = struct();
dataopts = meshopts;

%mesh parameters
meshopts.res = 12.5;
meshopts.nint = 1; % 1 == zero subdivisions, 2 == one subdivision, etc.
meshopts.octsubdiv = 2;

%data parameters
dataopts.res = 10;
dataopts.nint = 2;
dataopts.octsubdiv = 1;

%psuedo mesh parameters
pseudoOpts.res = meshopts.res;
pseudoOpts.nint = meshopts.nint;
pseudoOpts.octsubdiv = 1;

T = true;
F = false;
meshloadQ = T; %just makes it easier to switch back and forth between true and false
dataloadQ = T;
pseudoloadQ = T;
meshdataloadQ = T; %whether to check for and load intersection & barycentric data from previous run

%% generate mesh
disp('generating mesh')
mesh = datagen_setup(meshMethod,meshopts,meshloadQ);

toc; disp(' ')

%% generate data
disp('generating data')
data = datagen_setup(dataMethod,dataopts,dataloadQ);
ndatapts = size(data.pts,1);

toc; disp(' ')

%% generate pseudomesh
disp('generating pseudo mesh')
pseudo = datagen_setup(pseudoMethod,pseudoOpts,pseudoloadQ);

toc; disp(' ')

%% find intersecting facet of datapoints
disp('find intersecting facet for each datapoint')

r = 0.1; %Euclidean distance
meshdata.fname = ['mesh_' mesh.fname(1:end-4) '_data_' data.fname(1:end-4) '_100r--' int2str(100*r) '.mat'];

%find symmetrically equivalent octonions inside convex hull
% 	xyz = mesh.Ktr.main.pts;
% 	tess = mesh.Ktr.main.K;

% 	xyz = mesh.pts;
% 	tess = mesh.sphK;
% 	pts = data.pts;
% 	usv = mesh.usv;
% 	five = data.five;
% 	savename = [meshdata.fname(1:end-4) '_oint.mat'];
% 	oint = inhull_setup(pts,usv,xyz,tess,five,savename);

%% find nearest vertex of mesh to datapoint
disp('nearest neighbor search for all datapoints')
% 	if ~isempty(mesh.usv) && (size(mesh.pts,2) == 8)
% 		ptstemp = proj_down(mesh.pts,1e-6,mesh.usv);
% 	else
% 		ptstemp = mesh.pts;
% 	end
[nnList,nndistList] = dsearchn(mesh.pts,data.pts);

%remember to have each quaternion normalized to 1 before using get_omega
meshtemp = sqrt2norm(mesh.pts(nnList,:),'oct');
datatemp = sqrt2norm(data.pts,'oct');
nndistList = get_omega(meshtemp,datatemp);
% 	disp('convert mesh to 5DOF FZ')
%	meshfiveFZ = tofiveFZ(mesh.five,mesh.pts);

% 	disp('convert data to 5DOF FZ')
%	datafiveFZ = tofiveFZ(data.five,datapts);

% 	d_all = vertcat(meshfiveFZ.d);
% 	Y = vertcat(datafiveFZ.d);
d_all = q2rod(disorientation(vertcat(pseudo.five.q),'cubic'));
Y = q2rod(disorientation(vertcat(data.five.q),'cubic'));

%check for points within specified distance
[idxlist,D] = rangesearch(d_all,Y,r);
%[idxlist,D] = knnsearch(d_all,Y,'K',1000);
%idxlist = num2cell(idxlist,2);

%find nearest vertex of pseudomesh to datapoint
psdata.pts = zeros(size(data.pts));
psdata.five(ndatapts) = struct;
psdata.wmin = zeros(size(data.pts,1),1);
%for i = 1:ndatapts
for i = []
	%unpack point
	pt = data.pts(i,:);
	
	idx = idxlist{i};
	disp([int2str(length(idx)) ' octonion inputs to GBdist'])
	%olist = pseudo.pts(idx,:);
	olist = mesh.pts(idx,:);
    %  		olist = pseudo.pts;
	%olist = mesh.pts;
	
	ptrep = repmat(pt,size(olist,1),1); % matrix of repeated rows (pt)
	olistnorm = sqrt2norm(olist,'oct');
	ptrepnorm = sqrt2norm(ptrep,'oct');
	
	omega = GBdist([olistnorm ptrepnorm],32,false,false);
	
	[mintemp,minID] = min(omega);
	ptstemp = olist(minID,:);
	disp(mintemp)
	
	psdata.wmin(i) = mintemp;
	psdata.pts(i,:) = ptstemp;
	
	if mod(i,10) == 0
		disp(['i = ' int2str(i)])
	end

end

%project mesh and data together
tol = 1e-3;
[a,usv] = proj_down([mesh.pts;data.pts],tol);

datapts = data.pts;
meshpts = mesh.pts;

projdownQ = T;
if projdownQ
	datapts = proj_down(datapts,tol,usv);
	meshpts = proj_down(meshpts,tol,usv);
end

normQ = T;
if normQ
	datapts = normr(datapts);
	meshpts = normr(meshpts);
end

%compute intersecting facet IDs (might be zero, might have more than one)
tol2 = 1e-6;
intfacetIDs = intersect_facet(meshpts,mesh.sphK,datapts,tol2);

toc; disp(' ')

%% projection, sph. barycentric coordinates & interpolation

%initialize
databary = NaN(size(datapts));
facetprops = databary;
datainterp = NaN(size(datapts,1),1);
nnID = [];
ilist = [];
nonintDists = [];
databaryTemp = cell(1,size(datapts,1));
for i = 1:ndatapts
	datapt = datapts(i,:); %use down-projected data (and mesh)
	baryOK = false; %initialize
	if ~isempty(intfacetIDs{i})
		%setup
		intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
		vtxIDs = mesh.sphK(intfacetID,:);
		facet = meshpts(vtxIDs,:); %vertices of facet
		facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
		prop = data.props(i,:);
		
		baryType = 'planar'; %'spherical', 'planar'
		%% barycentric coordinates
		switch baryType
			case 'spherical'
				databary(i,:) = sphbary(datapt,facet); %need to save for inference input
				nonNegQ = all(databary(i,:) >= 0);
				greaterThanOneQ = sum(databary(i,:)) >= 1-1e-12;
				baryOK = nonNegQ && greaterThanOneQ;
				
			case 'planar'
				[~,databaryTemp{i}] = intersect_facet(facet,1:7,datapt,1e-12,false);
				databary(i,:) = databaryTemp{i}(1,:);
				nonNegQ = all(databary(i,:) >= -1e-12);
				equalToOneQ = abs(sum(databary(i,:)) - 1) < 1e-6
				baryOK = nonNegQ && equalToOneQ;
		end
		
		if baryOK
			%% interpolate using bary coords
			datainterp(i) = dot(databary(i,:),facetprops(i,:));
		else
			disp([num2str(databary(i,:),2) ' ... sum == ' num2str(sum(databary(i,:)))]);
		end
	end
	if ~baryOK
		disp(['i == ' int2str(i) ...
			'; no valid intersection, taking NN with dist = ' num2str(nndistList(i))])
		nonintDists = [nonintDists;nndistList(i)];
		nndistList(i) = NaN; %to distinguish interp vs. NN distances in plotting
		nnID = [nnID;nnList(i)]; %nearest neighbor indices
		ilist = [ilist;i];
		% 			datainterp(i) = mesh.props(k(i));
	end
end
fpath = fullfile('data',meshdata.fname);
save(fpath)

disp(['# non-intersections: ' int2str(length(nnID)) '/' int2str(ndatapts)])

%% plotting

% parity plot
% xmin = min(data.props);
% xmax = max(data.props);


interpplot(meshdata.fname)


%%
toc
disp('end run')
disp(' ')
%%
%-----------------------------CODE GRAVEYARD-------------------------------
%{

%}
