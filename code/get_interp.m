function [datainterp,databary,facetprops,facetIDs,barypars] = ...
	get_interp(mesh,data,intfacetIDs,barytype,barytol,NV)
arguments
	mesh struct {mustContainFields(mesh,{'pts','ppts','props','sphK'})}
	data struct {mustContainFields(data,{'pts','ppts','props'})}
	intfacetIDs cell
	barytype char {mustBeMember(barytype,{'planar','spherical'})} = 'planar'
	barytol double = double.empty
    NV.saveQ logical = false
    NV.savename char = 'temp.mat'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: Interpolate query point property values based on their
% spherical or planar barycentric coordinates relative to a mesh.
%
% Inputs:
%	mesh - struct containing octonions (pts), down-projected octonions
%	(ppts), property values (props), a triangulation (sphK), and a filename
%	(fname) for the "predictors" and "predictor responses"
%
%	data - struct containing octonions (pts), down-projected octonions
%	(ppts), property values which can be all NaNs (props), and a (fname)
%	for the query points and interpolated values
%
%   intfacetIDs - intersecting facet IDs of data.ppts relative to mesh.ppts
%   and mesh.sphK
%
%   barytype - type of barycentric interpolation, 'spherical' or 'planar'
%
%   barytol - tolerance for barycentric coordinates, defaults to 0.2 for
%   spherical and 1e-6 for planar. *
%
% Outputs:
%	datainterp - interpolated property values at the query points
%
%   databary - barycentric coordinates of query points relative to mesh
%
%   savename - catenatation of mesh.fname and data.fname where the
%   workspace is saved
%
%   varargout{1} - nndistList, NN euclidean distances
%
%   varargout{2} - nonintDists, NN euclidean distances for non-intersecting
%   points
%
% Usage:
%   [datainterp,databary,savename,varargout] =
%   get_interp(mesh,data,intfacetIDs,barytype,barytol);
%
% Dependencies:
%   get_omega.m
%
%   projray2hypersphere.m
%
%   sphbary.m
%    -projfacet2hyperplane.m
%     --projray2hyperplane.m
%
%   sqrt2norm.m
%
% Notes:
%	consider adding option for 'bary' vs. 'knn' interpolation
%
%   *0.2 was chosen as the spherical barytol default because it decreased
%   the interpolation error relative to planar, whereas lower tolerances
%   produced more non-intersections.
%--------------------------------------------------------------------------
if isempty(barytol)
    %assign default barytol
    switch barytype
        case 'spherical'
            barytol = 0.2;
        case 'planar'
            barytol = 1e-6;
    end
end

%unpack
meshpts = mesh.ppts;
datapts = data.ppts;

%check dimensions
if size(mesh.ppts,2) ~= size(data.ppts,2)
	errmsg = ['mesh.ppts and data.ppts dims should be equal, but mesh dim == ' ...
		int2str(size(mesh.ppts,2)) ' and data dim == ' int2str(size(data.ppts,2))];
	error(errmsg)
end

%save name
% savename = ['mesh_' mesh.fname(1:end-4) '_data_' data.fname];

%nearest neighbor list
nnList = dsearchn(meshpts,datapts);
nndistList = get_omega(sqrt2norm(mesh.pts(nnList,:)),sqrt2norm(data.pts));

%initialize
[databary,facetIDs,facetprops] = deal(nan(size(datapts)));
[datainterp,nonintDists] = deal(nan(size(datapts,1),1));
[nnID,ilist] = deal([]);

meannorm = mean(vecnorm(mesh.ppts,2,2));

ndatapts = size(datapts,1);
disp('loop through datapoints')
for i = 1:ndatapts
	datapt = datapts(i,:); %use down-projected data (and mesh)
	baryOK = false; %initialize
	if ~isempty(intfacetIDs{i})
		%% setup
		intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
		vtxIDs = mesh.sphK(intfacetID,:);
		facet = meshpts(vtxIDs,:); %vertices of facet
        facetIDs(i,:) = vtxIDs;
		facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
		
		%% barycentric coordinates
		switch barytype
			case 'spherical'
				databary(i,:) = sphbary(datapt,facet); %need to save for inference input
				nonNegQ = all(databary(i,:) >= -barytol);
				greaterThanOneQ = sum(databary(i,:)) >= meannorm-barytol;
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
			% interpolate using bary coords
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
		datainterp(i) = mesh.props(nnList(i)); %assign NN value
	end
end

%% Error Metrics
ids = setdiff(1:ndatapts,ilist); %ids = ~isnan(datainterp);

interp_errmetrics = get_errmetrics(datainterp(ids),data.props(ids));
nn_errmetrics = get_errmetrics(mesh.props(nnID),data.props(ilist));

errmetrics = get_errmetrics(...
    [datainterp(ids);mesh.props(nnID)],...
    [data.props(ids);data.props(ilist)]);

allnn_errmetrics = get_errmetrics(mesh.props(nnList),data.props);

interpRMSE = interp_errmetrics.rmse;
nnRMSE = nn_errmetrics.rmse;
totRMSE = errmetrics.rmse;
allnnRMSE = allnn_errmetrics.rmse;

nints = length(ids);
numnonints = length(ilist);
int_fraction = nints/(nints+numnonints);

barypars = var_names(errmetrics,interp_errmetrics,nn_errmetrics,allnn_errmetrics,...
    ilist,ids,nnList,nndistList,nonintDists,nints,numnonints,int_fraction);

% fpath = fullfile('data',savename);
% save(fpath,'-v7.3')

disp(' ')
disp(['# non-intersections: ' int2str(sum(~isnan((nnID)))) '/' int2str(ndatapts)])
disp(' ')
disp(['RMSE (J/m^2): interp == ' num2str(interpRMSE,'%3.4f') ', NN == ' num2str(nnRMSE,'%3.4f')])
disp(' ')
disp(['total RMSE: ' num2str(totRMSE,'%3.4f') ', all NN RMSE comparison: ' num2str(allnnRMSE,'%3.4f')])

%% Saving
if NV.saveQ
    save(NV.savename)
end

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

% % initialize
% databary = NaN(size(datapts));
% facetprops = databary;
% datainterp = NaN(size(datapts,1),1);
% nndistList = datainterp;
% nonintDists = datainterp;
% nnID = [];
% ilist = [];


% interpSE = (data.props(ids)-datainterp(ids)).^2;
% nnSE = (data.props(ilist)-mesh.props(nnID)).^2;
% totSE = [interpSE;nnSE];
% allnnSE = (data.props-mesh.props(nnList)).^2;
% interpRMSE = sqrt(mean(interpSE));
% nnRMSE = sqrt(mean(nnSE));
% totRMSE = sqrt(mean(totSE));
% allnnRMSE = sqrt(mean(allnnSE));

%calculate SE and RMSE values
%}
