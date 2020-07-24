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
meshopts.nint = 2; % 1 == zero subdivisions, 2 == one subdivision, etc.
meshopts.octsubdiv = 1;

%data parameters
dataopts.res = 10;
dataopts.nint = 2;
dataopts.octsubdiv = 1;

%psuedo mesh parameters
pseudoOpts.res = meshopts.res;
pseudoOpts.nint = meshopts.nint;
pseudoOpts.octsubdiv = 3;

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

r = 0.01; %Euclidean distance
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
	olist = pseudo.pts(idx,:);
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
datapts = proj_down(data.pts,1e-3,usv);
meshpts = proj_down(mesh.pts,1e-3,usv);

%compute intersecting facet IDs (might be zero, might have more than one)
tol2 = 1e-4;
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
for i = 1:ndatapts
	datapt = data.pts(i,:);
	if ~isempty(intfacetIDs{i})
		%setup
		intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
		vtxIDs = mesh.sphK(intfacetID,:);
		facet = mesh.pts(vtxIDs,:); %vertices of facet
		facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
		prop = data.props(i,:);
		
		%% spherical barycentric coordinates
		databary(i,:) = sphbary(datapt,facet); %need to save for inference input
		
		%% interpolate using sph. bary coords
		datainterp(i) = dot(databary(i,:),facetprops(i,:));
		
	else
		disp(['i == ' int2str(i) ...
			'; no intersection, taking NN with dist = ' num2str(nndistList(i))])
		nonintDists = [nonintDists;nndistList(i)];
		nndistList(i) = NaN;
		nnID = [nnID;nnList(i)]; %nearest neighbor indices
		ilist = [ilist;i];
		% 			datainterp(i) = mesh.props(k(i));
	end
end
save(meshdata.fname)

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

% mat = cell(d,1);
% for i = 1:d
% 	mat{i} = p
% end

%[x,y,z] = sphere;

%K = convhulln(smat);
%smat = smat(K);

	%plot3(0,0,0,'k*')
	%scale = 0.9;
	%surf(x*scale,y*scale,z*scale)
	%plot3(subcenter{:},'r.','MarkerSize',15)

%nmatTemp = cell([ndatapts,adjSize(1)]);
%p = cell([ndatapts,adjSize]);
%ddet = zeros([ndatapts,adjSize(1)]);

%snorm = smat./repelem(vecnorm(smat,2,2),1,d); %normalized vectors such that every point is on unit sphere (assumes centering at 0,0,0)
		%snorm = snorm + (0.01*rand(size(snorm))-0.005); %doesn't seem to make much of a difference
		%snorm = smat; %actually seems to take longer

%nn2 = tsearchn(snorm,TR,datanorm);

%dataList = repelem(1,d)+(0.05*rand(1,d)-0.025); %thought this might be a common thread for when no intersecting facets are found, but that doesn't seem to be the case.


	projMethod = 'numStabBary'; %'Bary', 'numStabBary'
	switch projMethod
		case 'Bary'

		case 'numStabBary'
			f = waitbar(0,['looping through ' int2str(adjSize(1)) ' neighboring facets']);
			%extract vertex points
			for j = 1:adjSize(1) %loop through facets
				for k = 1:adjSize(2) %loop through vertices of facet
					p{j,k} = snorm(facetsPtIDs(j,k),:);
				end
				
				%catenate vertices of each facet
				nmatTemp{j} = vertcat(p{j,:});
				
				if j == 13
					disp('')
				end
				%compute numerically stable barycentric coordinates
				lambda{j} = numStabBary(nmatTemp{j},datanorm);
				
				posQ(j) = all(lambda{j} >= 0);
				waitbar(j/adjSize(1),f);
			end
			close(f)

					if j == 244
						disp(' ')
					end

	%visual check that correct simplex was chosen
	% 	if d == 3
	% 		%adata = vertcat(a{i,j});
	% 		adata = a{id};
	% 		triPts = facetPts{i};
	% 		hold on
	% 		plot3(adata(:,1),adata(:,2),adata(:,3),'r*','MarkerSize',10)
	% 		plot3(triPts(:,1),triPts(:,2),triPts(:,3),'k.','MarkerSize',25)
	% 		hold off
	% 	end

		%trisurf(K(nn,:),snorm(:,1),snorm(:,2),snorm(:,3),0.5)
		%tetramesh(TR(nn,:),snorm)




%initialize
p = cell(ndatapts,1);
nmat = p;
nvec = p;
ddet = cell(ndatapts,1);
dmat = ddet;
	p = cell(adjSize);
	nmat = p;
	nvec = p;
	ddet = zeros(adjSize(1),1);
	dmat = cell(adjSize(1),1);

% nA = zeros(size(ropts,1),3);

% 		nA(:,1) = rand(size(ropts,1),1);
% 		nA(:,2) = 1-nA(:,1);

% nA = rand(size(ropts,1),3);

% %construct octonions
% for i = 1:size(ropts,1)
% 	qlist(i,:) = rod2q(ropts(i,:));
% 	nA(i,:) = nA(i,:)./norm(nA(i,:));
% 	o(i,:) = GBfive2oct(qlist(i,:),nA(i,:));
% end



%% load crystal symmetry
pgnum = 32; %30 === cubic symmetry 432, 32 === FCC (m-3m)

filenames = {'PGnames.mat','PGsymops.mat'};
for i = 1:length(filenames)
	file = dir(fullfile('**',filenames{i}));
	filepaths{i} = fullfile(file.folder,file.name);
	load(filepaths{i});
end

pgname = PG_names{pgnum};

qpt = Q{pgnum}; %choose point group symmetry
npt = length(qpt(:,1));

for i = 1:nmesh.pts
	%apply symmetry operators
	%qS{i} = qmult(qpt(i,:),q);
end

%% apply symmetry operators
%qSA = qmult(Si,qA);



checkAllFacetsQ = false; %check only facets attached to NN vertex (false) or check all facets (true)
if ~checkAllFacetsQ
	doublecheckFacetsQ = true; %check all facets for datapoints with no intersection found
end



delaunayQ = false; %not implemented correctly for true

if delaunayQ
	TR = delaunayn(mesh.pts); %takes longer than convhulln(), less simplices, one extra vertex
	disp(size(TR))
end

if delaunayQ
	nnList = dsearchn(mesh.pts,TR,dataList,[]);
else
	nnList = dsearchn(mesh.pts,dataList);
end

	if delaunayQ
		[row,col] = find(TR==nn);
		facetPtIDs = TR(row,:);
	else
		[row,col]=find(K==nn);
		facetPtIDs= K(row,:);
	end

%{
weighted average using barycentric coords (easy), octonion distance
(medium), or generalized barycentric coordinates of spherical triangle
(hard)
%}



if (d == 3 || d == 2) && strcmp(meshgenMethod,'hypersphereGrid')
	hold on
	
	plot3(mesh.pts(:,1),mesh.pts(:,2),mesh.pts(:,3),'b.','MarkerSize',15)
	
	if d == 2
		%account for d+1 embedding, assuming null dimension in first column
		trisurf(K,repelem(0,nmesh.pts),mesh.pts(:,2),mesh.pts(:,3))
	else
		trisurf(K,mesh.pts(:,1),mesh.pts(:,2),mesh.pts(:,3))
	end
	hold off
	
	axis equal
	ax = gca;
	ax.View = [8.3605e+01 1.3650e+01];
	%title(['[n x n x n], n = ',int2str(n),', number of pts in 7D sphere: ',...
	%	num2str(npts^(7/3),'%.1e')])
end

	if d == 3 && strcmp(meshgenMethod,'hypersphereGrid')
		hold on
		plot3([0 datanorm(:,1)*1.25],[0 datanorm(:,2)*1.25],[0 datanorm(:,3)*1.25],'k','LineWidth',2)
		plot3(mesh.pts(nn,1),mesh.pts(nn,2),mesh.pts(nn,3),'k*','MarkerSize',10)
		hold off
	end

%% Compute main convex hull
disp('computing convex hull')


disp(['# facets: ',int2str(size(K,1))])

toc; disp(' ')




%initalize
dataProj = cell(1,ndatapts); %projected data
facetPts = dataProj; %facet points
dataBary = dataProj; %barycentric data
facetID = dataProj; %intersecting facet ID



randlist = randi(ndatapts,1,iend+1-istart);

for i  = istart:iend
	j = randlist(i);
	disp(['---datapoint ',int2str(j)])
	datanorm = datameshList(j,:);
	nn = nnList(j);
	
	%find vertices of facets attached to NN vertex (or use all facets)
	[row,col]=find(K==nn);
	facetPtIDs= K(row,:);
	
	%compute projections
	[dataProj{i},facetPts{i},dataBary{i},facetID{i}] = projray2hypersphere(mesh.pts,facetPtIDs,datanorm);
	facetID{i} = row(facetID{i}); %convert from facetPtIDs index to K index
	
	mesh.ptsTemp = mesh.pts; %dummy variable to be able to sift through new NN's
	k = 0;
	while isempty(dataProj{i}) && k < nmesh.pts
		k = k+1;
		%remove previous NN
		mesh.ptsTemp(nn,:) = NaN(1,size(mesh.pts,2));
		
		%find next NN
		nn = dsearchn(mesh.ptsTemp,datanorm);
		
		%find facets attached to next NN
		[row,col]=find(K==nn);
		facetPtIDsNext= K(row,:);
		
		%compute projections
		[dataProj{i},facetPts{i},dataBary{i},facetID{i}] = projray2hypersphere(mesh.pts,facetPtIDsNext,datanorm);
	end
	
	if k > 1 && k < nmesh.pts-1
		disp(['intersection repeated up to ' int2str(k) '-th nearest neighhor'])
	elseif k == nmesh.pts
		disp('Looped through all facets.')
		if isempty(dataProj{i})
			disp('no intersecting facet found')
		end
	end
	disp(' ')
end

	%% project vertices onto hyperplane tangent to sphere at norm
	newvertices(i,:) = projfacet2hyperplane(data,vertices);


for i = istart:iend
	%get datapoint
	rand_id = randlist(i);
	data = datameshList(rand_id,:);
	vertices = facetPts{i};
	
	%% spherical barycentric coordinates
	sphbarycoords(i,:) = sphbary(data,vertices);
	
	%% interpolate using sph. bary coords
	datatrue = propList(rand_id);
	
% 	vtxIDs = K(row,:);
% 	meshdata = mesh.props(row
end

	rand_id = randlist(i);


% 	vtxIDs = K(row,:);
% 	meshdata = mesh.props(row

% meshinterp= zeros(size(meshList,1));

% databary = zeros(size(datameshList));


%[mesh.pts,~,K] = datagen(meshMethod,'mesh');




%note, I can get around the hyperorthant double-facet finding issue by
%taking the facet intersection point with the largest norm (i.e. furthest
%away from origin)


	dataprop = data.props(i);

% 	data = data.pts(i,:);


istart = 1; % set to other than 1 for troubleshooting specific datapoints
iend = 10;

res.data = 15;
nint.data = 1;
octsubdiv.data = 1;

res.mesh = 15;
nint.mesh = 1;
octsubdiv.mesh = 1;


ax = scatter(data.props,datainterp);
alpha(ax,0.1);


res = dataopts.res;
nint = dataopts.nint;
octsubdiv = dataopts.octsubdiv;
data.fname = get_fname(dataMethod,res,nint,octsubdiv);

if exist(data.fname,'file') ~= 0
	disp(data.fname)
	if ~strcmp(mesh.fname,data.fname)
		S = load(data.fname,'pts','props');
		data.pts = S.pts;
		data.props = S.props;
	else
		data.pts = S.pts;
		data.props = S.props;
	end
else
	%assumption is that mesh.props is a 1D array 2020-07-09
	[data.pts,data.props,data.sphK] = datagen(dataMethod,octsubdiv,'data',res,nint);
	
	% package for output
	pts = data.pts;
	props = data.props;
	temp = mesh.sphK;
	mesh.sphK = data.sphK;
	save(data.fname,'pts','props','mesh.sphK')
	mesh.sphK = temp; %reassign mesh.sphK to previous value
end

ndatapts = length(data.pts);
disp(['# datapoints = ',int2str(ndatapts)])


% set(gca,'xscale','log','yscale','log')


			% 			[k(i),datadist(i)] = dsearchn(mesh.pts,datapt);

% box on;
%
% set(gca,'boxstyle','full');

	for i = 1:size(five,1)
		%symmetrically equivalent octonions for each datapoint
		pgnum = 32;
		olist = osymsets(datapts,pgnum,usv);
		olist2 = vertcat(olist{:});
	end

	rmIDs = setdiff(1:length(pseudo),nnList2);



		%generate symmetrically equivalent octonions
		% 		osyms = osymsets(pt,32);
		% 		osyms = osyms{1};
		% % 		osyms = proj_down(osyms,1e-6,usv);
		
		% 		idx = idxlist{i};
				% 		olist = mesh.pts(idx,:);
		% 		olist = proj_up(mesh.pts,mesh.usv);

% if (exist(meshdata.fname,'file') ~= 0) && meshdataloadQ
% 	disp(meshdata.fname)
% 	load(meshdata.fname,'mesh','data','datainterp','nnID','psdata','ilist','meshopts','dataopts','nndistList','nonintDists')
% else


	%find minimum euclidean distance NN
	% 		[nnList2,nndistList2] = dsearchn(pseudo.pts,osyms); %assumption that Euclidean NN == spherical NN
	% 		[mymin,nnIDtemp] = min(nndistList2);
	% 		nnID = nnList2(nnIDtemp);
	
	%package output
	% 		psdata.pts(i,:) = pseudo.pts(nnID,:); % pseudo data
	% 		psdata.five(i) = pseudo.five(nnList2);
%}
