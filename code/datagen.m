function [meshList,propList,K,five,usv,Ktr] = ...
	datagen(sampleMethod,sampleType,opts)
arguments
	sampleMethod char = 'ocubo'
	sampleType char = 'mesh'
	opts struct = struct()
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Description:
%
% Inputs:
%		samplingMethod	===	type of data generation scheme to use.
%
%			'random'
%				randomly generates data from a uniform distribution in the
%				non-negative orthant.
%
%			'Kim2011'
%				reads in the iron grain boundary energy simulation data from
%				[1].
%
%			'5DOF'
%				produces a 5DOF mesh with tetrahedra in the misorientation FZ
%				and spherical triangles in the BP FZ.
%
%		sampleType		===	type of sampling scheme
%										'mesh' -- generate mesh only
%
%										'data' -- generate mesh and property values
%
% Outputs:
%		meshList		===	rows of vertices of mesh.
%
%		propList		===	column of property values (or empty array if
%								sampleType == 'mesh')
%
% Dependencies:
%		Kim2011_FeGBEnergy.txt
%		mesh5DOF.m
%		allcomb.m (if using sampleMethod == 'Rohrer2009')
%		addpathdir.m
%
% References
%		[1] H.K. Kim, W.S. Ko, H.J. Lee, S.G. Kim, B.J. Lee, An
%		identification scheme of grain boundaries and construction of a grain
%		boundary energy database, Scr. Mater. 64 (2011) 1152–1155.
%		https://doi.org/10.1016/j.scriptamat.2011.03.020.
%
% Notes:
%		Re-work argument inputs using "arguments" validation syntax
%--------------------------------------------------------------------------

%% setup
%get filenames, if any
switch sampleMethod
	case 'Kim2011'
		filelist = {'Kim2011_FeGBEnergy.txt'};
	case 'Olmsted2004'
		filelist = {'olm_octonion_list.txt','olm_properties.txt'};
end

%add filename directories to path, if any
if contains(sampleMethod,{'Kim2011','Olmsted2004'})
	%add folders to path
	addpathdir(filelist)
end

%% generate data
switch sampleMethod
	case 'random'
		% 		seed = 10;
		% 		rng(seed)
		rng('shuffle')
		
		%transpose to preserve rng sequence of pts with different number of datapts
		datatemp = rand(d+1,ndatapts).';
		meshList = zeros(d,ndatapts);
		
		for i = 1:ndatapts
			meshList(i,:) = datatemp(i,1:end-1)/vecnorm(datatemp(i,1:end-1));
		end
		propList = datatemp(:,end);
		
	case 'Kim2011'
		meshTable = readtable(filelist{1},'HeaderLines',9,'ReadVariableNames',true);
		varNames = meshTable.Properties.VariableNames; %Euler angles (misorientation), the polar & azimuth (inclination), GBE (mJ/m^2)
		datatemp = table2array(meshTable);
		
		meshList = datatemp(:,1:end-1);
		propList = datatemp(:,end);
		
	case 'Rohrer2009'
		%for now (2020-07-02) just gets the meshpoints of the grid they used,
		%rather than the triple junction, and then later assigns GB energy
		%based on BRK function
		
		step = 10; %degrees
		phi1_list = 0:step:90;
		cap_phi_list = cos(phi1_list);
		phi2_list = phi1_list;
		theta_list = cap_phi_list;
		phi_list = phi1_list;
		
		meshList = allcomb(phi1_list,cap_phi_list,phi2_list,theta_list,phi_list);
		
		%initialize
		npts = length(meshList);
		five(npts) = struct;
		five(1).q = [];
		five(1).nA = [];
		
		%convert to q,nA representation
		for i = 1:length(meshList)
			phi1 = meshList(i,1);
			cap_phi = meshList(i,2);
			phi2 = meshList(i,3);
			theta = meshList(i,4);
			phi = meshList(i,5);
			
			%the next few lines need to be checked to get them consistent with
			%convention used in Rohrer for angles
			five(i).q = rod2q([phi1,cap_phi,phi2]);
			[x,y,z] = sph2cart(theta,phi);
			five(i).nA = [x,y,z];
		end
		
	case 'Olmsted2004'
		
		meshList = readmatrix(filelist{1},'NumHeaderLines',1,'Delimiter',' ');
		meshList = meshList(:,1:8); %octonion representation, gets rid of a random NaN column..
		
		% 		propList = readmatrix(filelist{2},'NumHeaderLines',1,'Delimiter',' '); %properties
		% 		propList = propList(:,1); %1st column == GB energy
		
		five = GBoct2five(meshList,false);
		
	case '5DOF'
		%resolution in misorientation FZ
		% 		resDegrees = 10;
		%subdivisions of spherical triangles of BPs, n == 1 does nothing, n
		%== 2 does one subdivision, n == 3 does 2 subdivisions etc.
		% 		nint = 1;
		featureType = 'resolution';
		
	case {'5DOF_interior','5DOF_interior_pseudo'}
		featureType = 'interior';
		
	case {'5DOF_exterior','5DOF_exterior_pseudo'}
		featureType = 'exterior';
		
	case '5DOF_vtx'
		featureType = 'vertices';
		
	case '5DOF_vtx_deleteO'
		featureType = 'vtx_deleteO';
		
	case '5DOF_vtx_deleteOz'
		featureType = 'vtx_deleteOz';
		
	case '5DOF_ext_deleteO'
		featureType = 'ext_deleteO';
		
	case '5DOF_misFZfeatures'
		featureType = 'misFZfeatures';
		
	case '5DOF_oct_vtx'
		featureType = 'vtx_deleteOz';
		%the rest gets taken care of later
		
	case {'5DOF_hsphext','5DOF_hsphext_pseudo'}
		featureType = 'vtx_deleteOz';
		
	case {'5DOF_exterior_hsphext','5DOF_exterior_hsphext_pseudo'}
		featureType = 'exterior';
		
	case {'ocubo','ocubo_hsphext'}
		featureType = 'ocubo';
		%unpack options
		n = opts.ocuboOpts.n;
		method = opts.ocuboOpts.method;
		sidelength = opts.ocuboOpts.sidelength;
		seed = opts.ocuboOpts.seed;
		
		%get cubochorically sampled octonions
		meshList = get_ocubo(n,method,sidelength,seed);
		five = GBoct2five(meshList,false);
end

%5DOF cases
if contains(sampleMethod,'5DOF')
	ctrcuspQ = false;
	[five,sept] = mesh5DOF(featureType,ctrcuspQ,opts.res,opts.nint);
	meshList = sept;
end

sz = size(meshList); %store size of degenerate meshList
%reduce to unique set of points
[~,IA] = uniquetol(round(meshList,6),1e-3,'ByRows',true);
meshList = meshList(IA,:);
%pare down five to a unique set of points
five = five(IA);

if contains(sampleMethod,'5DOF_')
	savename = [sampleMethod(6:end) '_pairmin.mat'];
else
	savename = [sampleMethod '_pairmin.mat'];
end

if size(meshList,2) == 7
	meshList = [meshList zeros(size(meshList,1),1)];
end

if contains(sampleMethod,'oct_vtx')
	NVpairs = {'o2addQ',true,'method','pairwise','wtol',1e-3}; %method can be 'standard' or 'pairwise'
else
	NVpairs = {'o2addQ',false,'method','pairwise','wtol',1e-3};
end

method = 1;
switch method
	case 1
		[meshList,~,five,~] = get_octpairs(meshList,savename,NVpairs{:}); %find a way to not call this for 'data'
		
	case 2
		o1tmp = meshList(1,:);
		o2 = meshList(2:end,:);
		o1rep = repmat(o1tmp,size(o2,1),1);
		[~,o2] = GBdist4(o1rep,meshList(2:end,:),32,'norm');
		meshList = [o1tmp;vertcat(o2{:})];
		five = GBoct2five(meshList);
end

[meshList,usv] = proj_down(meshList,1e-4,struct.empty,'zeroQ',true);

% 	if strcmp(sampleMethod,'ocubo') %might need a way to correlate back to original dataset for e.g. Rohrer2009
%reduce to unique set of points
[~,IA] = uniquetol(round(meshList,6),1e-3,'ByRows',true);
meshList = meshList(IA,:);
%pare down five to a unique set of points
five = five(IA);
% 	end
if size(meshList,2) == 7
	projupQ = true;
else
	projupQ = false;
end
%% Subdivide octonions, convex hull
if contains(sampleMethod,'hsphext')
	[Ktr,K,meshList] = hsphext_subdiv(meshList,opts.octsubdiv);
	
elseif opts.octsubdiv > 1
	[Ktr,K,meshList] = hypersphere_subdiv(meshList,[],opts.octsubdiv); %originally had sphK
	
elseif (exist('sphK','var') ~= 0)
	%create K if it exists & is empty
	if isempty(opts.sphK)
		if contains(sampleMethod,'hsphext')
			[Ktr,K,meshList] = hsphext_subdiv(meshList,1);
		else
			K = sphconvhulln(meshList);
		end
	else
		K = opts.sphK;
	end
	
else
	K = sphconvhulln(meshList);
end

if projupQ
	%project up
	meshList = proj_up(meshList,usv);
end

if contains(sampleMethod,'hsphext') || opts.octsubdiv > 1
	meshList = get_octpairs(meshList);
	five = GBoct2five(meshList,false);
end

%package geometry into "five" (e.g. 'A', 'O', 'AC', etc.)
geomQ = true;
if geomQ
	geometry = findgeometry(disorientation(vertcat(five.q),'cubic'));
end
[five.geometry] = geometry{:};

%normalize octonions
assert(size(meshList,2) == 8,['meshList should have 8 columns, not ' int2str(size(meshList,2))])
meshList = normr(meshList);

disp(' ')
disp(['total of ' int2str(size(meshList,1)) ' after oct subdivision.'])

%% Compute GB Energies
if strcmp(sampleType,'data')
	if contains(sampleMethod,{'5DOF','ocubo','Rohrer2009','Olmsted2004'})
		disp('GB5DOF')
		propList = GB5DOF_setup(five);
	else
		propList = [];
		warning('propList is empty. Ok if sampleType == mesh. Otherwise, did you add to switch case, or is "five" empty?')
	end
end
if exist('usv','var') == 0
	usv = [];
end
if exist('Ktr','var') == 0
	Ktr = [];
end

end %datagen
%----------------------------CODE GRAVEYARD--------------------------------
%{




	case '5DOF'
		%resolution in misorientation FZ
% 		resDegrees = 10;
		%subdivisions of spherical triangles of BPs, n == 1 does nothing, n
		%== 2 does one subdivision, n == 3 does 2 subdivisions etc.
% 		nint = 1;
		featureType = 'resolution';
		
	case '5DOF_interior'
% 		resDegrees = 10;
% 		nint = 3;
		featureType = 'interior';
		
	case '5DOF_exterior'
% 		resDegrees = 10;
% 		nint = 2;
		featureType = 'exterior';
		
	case '5DOF_vtx'
% 		resDegrees = 5;
% 		nint = 1;
		featureType = 'vertices';
		
	case '5DOF_vtx_deleteO'
% 		resDegrees = 5;
% 		nint = 1;
		featureType = 'vtx_deleteO';
		
	case '5DOF_ext_deleteO'
% 		resDegrees = 5;
% 		nint = 1;
		featureType = 'ext_deleteO';
		
	case '5DOF_misFZfeatures'
% 		resDegrees = 5;
% 		nint = 1;
		featureType = 'misFZfeatures';


	if isempty(sphK)
		[Ktr,K,meshList] = hypersphere_subdiv(meshList,[],octsubdiv);
	else
		[Ktr,K,meshList] = hypersphere_subdiv(meshList,sphK,octsubdiv);
	end


sz = size(meshList); %store size of degenerate meshList
%reduce to unique set of points
[meshList,IA] = uniquetol(meshList,1e-6,'ByRows',true);

%pare down five to a unique set of points
[row,~] = ind2sub(sz,IA);
% row = sort(row);
five = five(row);


[row,~] = ind2sub(sz,IA);


	%create K if it exists & is empty
	if ~isempty(sphK)
		K = sphconvhulln(meshList);
		% 		K = convhulln(meshList);
	end

if strcmp(sampleMethod,'oct_vtx')
	%project back into d-dimensional space
	meshList = padarray(meshList,[0 ndegdim],'post')*v'+mean(pts);
end

% 	load('octvtx_pairmin.mat','five','pts','sphK','usv','avg')
% 	meshList = pts;
	%remove null dimensions
	% 		[pts,V,avg] = proj_down(pts);
	% 		sphK = sphconvhulln(pts);

% if strcmp(sampleMethod,{'5DOF_oct_vtx','5DOF_hsphext')
% 	meshList = pts; %just make sure to output V
% end


		meshList2 = meshList;


switch sampleMethod
	case {'5DOF','5DOF_interior','5DOF_exterior','5DOF_vtx','5DOF_vtx_deleteO',...
			'5DOF_vtx_deleteOz','5DOF_ext_deleteO','5DOF_misFZfeatures',...
			'5DOF_oct_vtx','5DOF_hsphext','5DOF_exterior_hsphext'}
		ctrcuspQ = false;
		[five,sept,o] = mesh5DOF(featureType,ctrcuspQ,res,nint);
		meshList = sept;
end


	switch sampleMethod
		case {'5DOF','5DOF_interior','5DOF_exterior','5DOF_vtx',...
				'5DOF_misFZfeatures','Rohrer2009','5DOF_vtx_deleteO',...
				'5DOF_vtx_deleteOz','5DOF_oct_vtx','5DOF_hsphext',...
				'5DOF_exterior_hsphext'}
			disp('GB5DOF')
			propList = GB5DOF_setup(five);
			
			
		otherwise
			propList = [];
			warning('propList is empty. Ok if sampleType == mesh. Otherwise, did you add to switch case, or is "five" empty?')
	end

% if any([strcmp(sampleMethod,'5DOF_oct_vtx'),contains(sampleMethod,'hsphext')])
%get symmetrized octonions with respect to two points ('O' and
%'interior', both +z)
savename = [sampleMethod(6:end) '_pairmin.mat'];
if size(meshList,2) == 7
	meshList = [meshList zeros(size(meshList,1),1)];
end
o2addQ = true;
[meshList,usv,five,~,~] = get_octpairs(meshList,five,savename,o2addQ);
meshList = proj_down(meshList,1e-6,usv);
% end

% 	if any([strcmp(sampleMethod,'5DOF_oct_vtx'),contains(sampleMethod,'hsphext')])
		%restore null dimensions for oct --> five
		meshList = proj_up(meshList,usv);
% 	end


load_method = 'evalc';
switch load_method
	case 'evalc'
		for i = 1:length(vars)
			var = vars{i};
			val = opts.(var); %#ok<NASGU> %temporary value of vName
			evalc([var '= val']); %assign temp value to a short name
		end
	case 'manual'
		res = opts.res;
		nint = opts.nint;
		sphK = opts.sphK;
		ocuboOpts = opts.ocuboOpts;
end


%package geometry into "five" (e.g. 'A', 'O', 'AC', etc.)
geomQ = true;
if geomQ
	for i = 1:length(five)
		five(i).geometry = findgeometry(disorientation(five(i).q,'cubic'));
	end
end

	opts.octsubdiv = 1
	opts.res = []; %double
	opts.nint = []; %double
	opts.sphK = []; % double
	opts.ocuboOpts = []; % struct

	%checking to see if this increases the interpolation accuracy - it did not
	NV = {'o2addQ',false,'method','pairwise','wtol',1e-3};
	meshList = get_octpairs(meshList,savename,NV{:});


	% 	if any([strcmp(sampleMethod,'5DOF_oct_vtx'),contains(sampleMethod,'hsphext')])
	%restore null dimensions for oct --> five

% save('temp.mat')



	% 	if isempty(meshListTemp)
	% 		warning('creating new usv')
	% 		[meshList,usv] = proj_down(meshList,1e-6);
	% 	end


%originally implemented for just '5DOF_oct_vtx' and 'hsphext'
%get symmetrized octonions with respect to two points ('O' and
%'interior', both +z)
if contains(sampleMethod,'5DOF_')
	savename = [sampleMethod(6:end) '_pairmin.mat'];
else
	savename = [sampleMethod '_pairmin.mat'];
end




pseudoQ = contains(sampleMethod,'pseudo');

if pseudoQ
	str_id = strfind(sampleMethod,'_pseudo');
	sampleMethod = sampleMethod(1:str_id-1); %assumes _pseudo is at end, and remove it
end


	elseif pseudoQ
		tricollapseQ = false;
		[~,K,meshList] = hypersphere_subdiv(meshList,[],opts.octsubdiv,tricollapseQ);
		Ktr = [];



if ~pseudoQ
	else
	try
		load(savename,'octvtx','usv')
		% 	meshListTemp = proj_down(meshList,1e-6,usv);
		[meshList,usv] = proj_down(octvtx,1e-6,usv,'zeroQ',true);
	catch
		warning("error with meshList = proj_down(octvtx,1e-6,usv,'zeroQ',true); or load")
		[meshList,usv] = proj_down(meshList,1e-6,struct.empty,'zeroQ',true);
	end
end

if contains(sampleMethod,'hsphext') || opts.octsubdiv > 1
	%renormalize each quaternion (i.e. bring back into space of rotations)
	% 	meshList(:,1:4) = normr(meshList(:,1:4));
	% 	meshList(:,5:8) = normr(meshList(:,5:8));
	
		% 	meshListTmp = get_octpairs(meshListTmp);

	
	meshList = get_octpairs(meshList);
	five = GBoct2five(meshList,false);
end


%}
