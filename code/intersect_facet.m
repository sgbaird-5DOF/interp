function [intfacetIDs,dataBary] = intersect_facet(pts,K,datalist,tol,maxnormQ,baryMethod)
arguments
	pts double {mustBeFinite,mustBeReal}
	K double {mustBeFinite,mustBeReal}
	datalist double {mustBeFinite,mustBeReal}
	tol(1,1) double {mustBeFinite,mustBeReal} = 1e-6
	maxnormQ(1,1) logical = false
	baryMethod char {mustBeMember(baryMethod,{'planar','spherical'})} = 'planar'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-06
%
% Description: Project a ray (each row of pts) onto each facet connected to
% the nearest neighbor of that ray, and compute the barycentric coordinates
% of the projected datapoint. If all coordinates of a facet are positive,
% mark that facet as an intersecting facet. If no intersecting is found
% with the first nearest neighbor search, continue looking at next nearest
% neighbors until an intersecting facet has been found or all facets have
% been looped through.
%
% Inputs:
%		pts			===	rows of points that fall on a unit sphere, must
%								match indices in K
%
%		K				===	convex hull triangulation of pts.
%
%		datalist		===	rows of datapoints that represent rays from origin
%								to the datapoint to be projected onto the facet.
%
% Outputs:
%		intfacetIDs ===	intersecting facet IDs (as in row index of K).
%								Returns NaN if no intersecting facet is found, even
%								after looping through all facets instead of just
%								ones connected to first NN. Facet with
%								largest norm returned if multiple facets found.
%
% Dependencies:
%		projray2hypersphere.m
%			numStabBary.m (optional)
%
% Notes:
%		Consider changing it to only look at the k-nearest neighbors before
%		just projecting onto all remaining facets since there are a lot of
%		repeats with using the next nearest neighbor approach for all facets.
%		Alternatively, keep track of which facets have already been
%		considered and remove them from the list of intersecting facets
%		before doing the projection.
%
%		To relax the requirement that pts need to close to on the unit
%		sphere, then there's a set of lines in projray2hypersphere.m that can
%		be changed or removed.
%
% 		if (t(j) <= 0.1) || (t(j) >= 2.1)
% 			posQ(j) = 0;
% 			continue
% 		end
%
%		Something wrong with 'spherical' still I think (too many
%		intersections) 2020-07-21
%
%--------------------------------------------------------------------------

%% find nearest vertex for each datapoint
nnList = dsearchn(pts,datalist);

nmeshpts = size(pts,1);
ndatapts = size(datalist,1);

%% initalize
dataProj = cell(1,ndatapts); %projected data
facetPts = dataProj; %facet points
dataBary = dataProj; %barycentric data
subfacetIDs = dataProj; %sub IDs (i.e. IDs from sublist of facets)
intfacetIDs = dataProj; % facet IDs that can be used to index into K
t = dataProj;

%textwaitbar setup
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateProgress);
N=ndatapts;
p=1;
reverseStr = '';
nreps = floor(N/20);
nreps2 = nreps;

%% loop through datapts
%disp('intersect_facet ')
parfor i  = 1:ndatapts % for parallelized, use parfor
	%% first NN projection
	data = datalist(i,:);
	nn = nnList(i);
	
	%find vertices of facets attached to NN vertex (or use all facets)
	[row,col]=find(K==nn);
	facetPtIDs= K(row,:);
	
	%compute projections
	switch baryMethod
		case 'planar'
            [dataProj{i},facetPts{i},dataBary{i},subfacetIDs{i},t{i}] = ...
				projray2hypersphere(pts,facetPtIDs,data,tol,maxnormQ);
		case 'spherical'
			[dataBary{i},subfacetIDs{i}] = sphbary_setup(pts,facetPtIDs,data,tol); %I think this is buggy 2020-07-16
				
	end
	
	%% keep using next NNs if facet not found
	ptsTemp = pts; %dummy variable to be able to sift through new NN's
	k = 0;
	oldrow = row;
	nnMax = 5; %nnMax = nmeshpts;
	while isempty(subfacetIDs{i}) && k < nnMax
		k = k+1;
		%remove previous NN
		ptsTemp(nn,:) = NaN(1,size(pts,2));
		
		%find next NN
		nn = dsearchn(ptsTemp,data);
		
		%find facets attached to next NN
		[row,~]=find(K==nn);
		
		rownext = setdiff(row,oldrow);
		oldrow = [row;oldrow];
		
		if ~isempty(rownext)
			facetPtIDsNext= K(rownext,:);
			
			%compute projections
			switch baryMethod
				case 'planar'
                    [dataProj{i},facetPts{i},dataBary{i},subfacetIDs{i},t{i}] = ...
						projray2hypersphere(pts,facetPtIDsNext,data,tol,maxnormQ);
				case 'spherical'
					[dataBary{i},subfacetIDs{i}] = sphbary_setup(pts,facetPtIDs,data,tol);
			end
		end
	end
	
	if ~isempty(subfacetIDs{i})
		intfacetIDs{i} = row(subfacetIDs{i}); %convert from facetPtIDs or facetPtIDsNext index to K index
	else
		intfacetIDs{i} = [];
	end
	
	% 	%info about next NN search
	% 	if k > 1 && k < nmeshpts-1
	% 		disp(['---datapoint ',int2str(i)])
	% 		disp(['intersection repeated up to ' int2str(k) '-th nearest neighhor'])
	% 		disp(' ')
	% 	elseif k == nmeshpts
	% 		disp('Looped through all facets.')
	% 		if isempty(dataProj{i})
	% 			disp('no intersecting facet found')
	% 		end
	% 		disp(' ')
	% 	end
	
	% 	disp(t{i})
	if mod(i,nreps2) == 0
		send(D,i);
	end
end

	function nUpdateProgress(~)
		percentDone = 100*p/N;
		msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
		fprintf([reverseStr, msg]);
		reverseStr = repmat(sprintf('\b'), 1, length(msg));
		p = p + nreps;
	end

end

%-------------------------------CODE GRAVEYARD-----------------------------
%{
		facetPtIDsNext= setdiff(K(row,:),facetPtIDsOld);

		[row2,~] = find(ismember(K(row,:),unique(facetPtIDsOld)));
		row2 = row(~ismember(1:length(row),row2));

% 	facetPtIDsOld = facetPtIDs;

% 			facetPtIDsOld = [facetPtIDsOld;facetPtIDsNext]; %#ok<AGROW>


%}
