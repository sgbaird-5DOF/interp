function [datainterp,databary,savename,varargout] = ...
	get_interp(mesh,data,intfacetIDs,barytype,barytol)
arguments
	mesh struct {mustContainFields(mesh,{'pts','props','sphK','fname'})}
	data struct {mustContainFields(data,{'pts','props','fname'})}
	intfacetIDs cell
	barytype char {mustBeMember(barytype,{'planar','spherical'})} = 'planar'
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
%	a - a
%
% Outputs:
%	b - b
%
% Usage:
%		a = b(a);
%
% Dependencies:
%		*
%
% Notes:
%		consider adding option for 'bary' vs. 'knn' interpolation
%--------------------------------------------------------------------------
meshpts = mesh.ppts;
datapts = data.ppts;

if size(mesh.ppts,2) ~= size(data.ppts,2)
	errmsg = ['mesh.pts and data.pts dims should be equal, but mesh dim == ' ...
		int2str(size(mesh.ppts,2)) ' and data dim == ' int2str(size(data.ppts,2))];
	error(errmsg)
end

%save name
savename = ['mesh_' mesh.fname(1:end-4) '_data_' data.fname];

%nearest neighbor list
nnList = dsearchn(meshpts,datapts);
nndistList = get_omega(sqrt2norm(mesh.pts(nnList,:)),sqrt2norm(data.pts));

%initialize
databary = NaN(size(datapts));
facetprops = databary;
datainterp = NaN(size(datapts,1),1);
% nndistList = datainterp;
nonintDists = datainterp;
nnID = [];
ilist = [];

ndatapts = size(datapts,1);
disp('loop through datapoints')
for i = 1:ndatapts
	datapt = datapts(i,:); %use down-projected data (and mesh)
	baryOK = false; %initialize
	if ~isempty(intfacetIDs{i})
		%setup
		intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
		vtxIDs = mesh.sphK(intfacetID,:);
		facet = meshpts(vtxIDs,:); %vertices of facet
		facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
		
		%% barycentric coordinates
		switch barytype
			case 'spherical'
				databary(i,:) = sphbary(datapt,facet); %need to save for inference input
				nonNegQ = all(databary(i,:) >= -barytol);
				greaterThanOneQ = sum(databary(i,:)) >= 1-barytol;
				numcheck = all(~isnan(databary(i,:)) & ~isinf(databary(i,:)));
				baryOK = nonNegQ && greaterThanOneQ && numcheck;
				
			case 'planar'
				[~,~,databaryTemp,~,~] = projray2hypersphere(facet,1:7,datapt,barytol,true);
				if ~isempty(databaryTemp)
					databary(i,:) = databaryTemp;
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
		nnID = [nnID nnList(i)]; %#ok<AGROW> %nearest neighbor indices
		ilist = [ilist i]; %#ok<AGROW> % possible to separate out making baryOK a logical array & using 2 for loops
		% 			datainterp(i) = mesh.props(k(i));
	end
end

%calculate SE and RMSE values
ids = ~isnan(datainterp);
interpSE = (data.props(ids)-datainterp(ids)).^2;
nnSE = (data.props(ilist)-mesh.props(nnID)).^2;
totSE = [interpSE;nnSE];
allnnSE = (data.props-mesh.props(nnList)).^2;

interpRMSE = sqrt(mean(interpSE));
nnRMSE = sqrt(mean(nnSE));
totRMSE = sqrt(mean(totSE));
allnnRMSE = sqrt(mean(allnnSE));

fpath = fullfile('data',savename);
save(fpath)

disp(' ')
disp(['# non-intersections: ' int2str(sum(~isnan((nnID)))) '/' int2str(ndatapts)])
disp(' ')
disp(['RMSE (J/m^2): interp == ' num2str(interpRMSE,'%3.4f') ', NN == ' num2str(nnRMSE,'%3.4f')])
disp(' ')
disp(['total RMSE: ' num2str(totRMSE,'%3.4f') ', all NN RMSE comparison: ' num2str(allnnRMSE,'%3.4f')])
varargout = {nndistList,nonintDists,intfacetIDs};

end

%-----------------------------CODE GRAVEYARD-------------------------------
%{

datapts = data.pts;
meshpts = mesh.pts;

normQ = true;
if normQ
	%normalize points
	datapts = normr(datapts);
	meshpts = normr(meshpts);
end

%get intersecting facets
maxnormQ = false;
disp('get intersecting facets')
intfacetIDs = intersect_facet(meshpts,mesh.K,datapts,inttol,maxnormQ);

% 		prop = data.props(i,:);

%}
