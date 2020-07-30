function [nndistList,databary] = get_interp(mesh,data,inttol,barytol)
arguments
	mesh struct {mustContainFields(mesh,{'pts','props','fname'})}
	data struct {mustContainFields(data,{'pts','props','fname'})}
	inttol(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 1e-6
	barytol(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 1e-6
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: get barycentric interpolation values
% 
% Inputs:
%		a -	a
%
% Outputs:
%		b -	b
%
% Usage:
%		a = b(a);
%
% Dependencies:
%		*
%
% Notes:
%		*
%--------------------------------------------------------------------------

if size(mesh.pts,2) ~= size(data.pts,2)
	errmsg = ['mesh.pts and data.pts dims should be equal, but mesh dim == ' ...
		int2str(size(mesh.pts,2)) ' and data dim == ' int2str(size(data.pts,2))];
	error(errmsg)
end

datapts = data.pts;
meshpts = mesh.pts;

normQ = true;
if normQ
	%normalize points
	datapts = normr(datapts);
	meshpts = normr(meshpts);
end

mesh.K = sphconvhulln(mesh.pts,true);

%get intersecting facets
maxnormQ = false;
disp('get intersecting facets')
intfacetIDs = intersect_facet(meshpts,mesh.K,datapts,inttol,maxnormQ);

%save name
meshdata.fname = ['mesh_' mesh.fname(1:end-4) '_data_' data.fname];

%nearest neighbor list
nnList = dsearchn(meshpts,datapts);

%initialize
databary = NaN(size(datapts));
facetprops = databary;
datainterp = NaN(size(datapts,1),1);
nndistList = datainterp;
nonintDists = datainterp;
nnID = [];
ilist = [];

ndatapts = size(datapts,1);
% databaryTemp = cell(1,size(datapts,1));
disp('loop through datapoints')
for i = 1:ndatapts
	
	datapt = datapts(i,:); %use down-projected data (and mesh)
	baryOK = false; %initialize
	
	if ~isempty(intfacetIDs{i})
		%setup
		intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
		vtxIDs = mesh.K(intfacetID,:);
		facet = meshpts(vtxIDs,:); %vertices of facet
		facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
		prop = data.props(i,:);
		
		baryType = 'spherical'; %'spherical', 'planar'
		%% barycentric coordinates
		switch baryType
			case 'spherical'
				databary(i,:) = sphbary(datapt,facet); %need to save for inference input
				nonNegQ = all(databary(i,:) >= -barytol);
				greaterThanOneQ = sum(databary(i,:)) >= 1-barytol;
				numcheck = all(~isnan(databary(i,:)) & ~isinf(databary(i,:)));
				baryOK = nonNegQ && greaterThanOneQ && numcheck;
				
			case 'planar'
				[~,databaryTemp] = intersect_facet(facet,1:7,datapt,inttol,true);
				if ~isempty(databaryTemp{1})
					databary(i,:) = databaryTemp{1};
					nonNegQ = all(databary(i,:) >= -barytol);
					equalToOneQ = abs(sum(databary(i,:)) - 1) < barytol;
					numcheck = all(~isnan(databary(i,:)) & ~isinf(databary(i,:)));
					baryOK = nonNegQ && equalToOneQ && numcheck;
				end
		end
		
		if baryOK
			%% interpolate using bary coords
			datainterp(i) = dot(databary(i,:),facetprops(i,:));
		else
			disp([num2str(databary(i,:),2) ' ... sum == ' num2str(sum(databary(i,:)),10)]);
		end
	end
	if ~baryOK
		disp(['i == ' int2str(i) ...
			'; no valid intersection, taking NN with dist = ' num2str(nndistList(i))])
		nonintDists(i) = nndistList(i);
		nndistList(i) = NaN; %to distinguish interp vs. NN distances in plotting
		nnID = [nnID nnList(i)]; %nearest neighbor indices
		ilist = [ilist i]; % possible to separate out making baryOK a logical array & using 2 for loops
		% 			datainterp(i) = mesh.props(k(i));
	end
end
fpath = fullfile('data',meshdata.fname);
save(fpath)


disp(['# non-intersections: ' int2str(sum(~isnan((nnID)))) '/' int2str(ndatapts)])

end