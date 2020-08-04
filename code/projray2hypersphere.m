function [dataProj,facetPts,dataBary,facetIDs,tvals] = ...
	projray2hypersphere(meshpts,facetPtIDs,datanorm,tol,maxnormQ,invmethod)
arguments
	meshpts double
	facetPtIDs double
	datanorm double
	tol(1,1) double = 1e-6
	maxnormQ(1,1) logical = false
	invmethod char {mustBeMember(invmethod,{'mldivide','pinv','extendedCross'})} = 'mldivide'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-24
%
% Description: project ray to hypersphere, compute barycentric coordinates,
% compute intersecting facet
%
% Inputs:
%
% Outputs:
%		facetIDs		===	list of facet IDs (i.e. rows of K triangulation)
%								for which an intersection was found.
%
% Dependencies:
%		numStabBary.m (if using invMethod == 'extendedCross')
%
% References:
%		https://math.stackexchange.com/q/1256236/798661
%
% Notes:
%		To relax the requirement that pts need to be close to on the unit
%		sphere, then there's a set of lines in projray2hypersphere.m that can
%		be changed or removed.
%
% 		if (t(j) <= 0.1) || (t(j) >= 2.1)
% 			posQ(j) = 0;
% 			continue
% 		end
%--------------------------------------------------------------------------

% invmethod = 'mldivide'; %'pinv', 'mldivide', 'extendedCross'

adjSize = size(facetPtIDs);
%%
%setup vectors and matrices for ray projection
p = cell(adjSize);
nmatTemp = cell(adjSize);
nmat = p;
nvec = p;
ddet = zeros(adjSize(1),1);
dmat = cell(adjSize(1),1);

for j = 1:adjSize(1) %loop through facets
	for k = 1:adjSize(2) %loop through vertices of facet
		%package vertices
		p{j,k} = meshpts(facetPtIDs(j,k),:);
	end
	%package vertex matrix
	nmatTemp{j} = vertcat(p{j,:});
	dmat{j} = nmatTemp{j};
	for k = 1:adjSize(2)
		nmat{j,k} = nmatTemp{j};
		nmat{j,k}(:,k) = 1;
	end
	for k = 1:adjSize(2) %loop through dimensions
		nvec{j}(k) = det(nmat{j,k});
	end
	ddet(j) = det(dmat{j});
end


%%
%project ray onto each facet
a = cell([adjSize(1),1]);
lambda = a; %initialize
posQ = zeros([adjSize(1),1]);

warnID1 = 'MATLAB:nearlySingularMatrix';
warnID2 = 'MATLAB:singularMatrix';
warnID3 = 'symbolic:mldivide:InconsistentSystem';

warning('off',warnID1)
warning('off',warnID2)
warning('off',warnID3)

for j = 1:adjSize(1)
	if ddet(j) ~= 0
		t(j) = ddet(j)/dot(datanorm,nvec{j});
		
		tol2 = 1e-12;
		if (t(j) <= tol2) || (t(j) >= 1/tol2)
			posQ(j) = 0;
			continue
		end
		
		a{j} = -datanorm*-t(j);
		
		switch invmethod
			case 'mldivide'
				lambda{j} = (dmat{j}'\a{j}')'; %barycentric coordinates
				posQ(j) = all(lambda{j} >= -tol) && ~any(lambda{j} == Inf) && (sum(lambda{j}) >= 1-tol); % && any(c{j} > 0)
			case 'extendedCross'
				%compute numerically stable barycentric coordinates
				lambda{j} = numStabBary(nmatTemp{j},a{j});
% 				disp(lambda{j})
				posQ(j) = all(lambda{j} >= -tol) && (sum(lambda{j}) >= 1-tol);
		end
	end
end

warning('on',warnID1)
warning('on',warnID2)
warning('on',warnID3)

posQtot = sum(posQ);

%%
% find facet that intersects ray and compile relevant data for datapoint
if posQtot > 0
	
	idList = find(posQ > 0);
	
	if (length(idList) >= 2) && maxnormQ
		%disp('taking datapoint with largest norm')
		[~, maxpos] = max(vecnorm(vertcat(a{idList})'));
		id = idList(maxpos);
	else
		id = idList;
	end
	
	dataProj = a{id}; %datapoint projected onto intersecting facet
	facetPts = vertcat(p{id,:}); %each row is a vertex
	dataBary = lambda{id}; %barycentric datapoints (non-generalized)
	facetIDs = id;
	tvals = t(logical(posQ));
	%disp(dataBary{i});
	
else
	dataProj = [];
	facetPts = [];
	dataBary = [];
	facetIDs = [];
	tvals = [];
end
end %projray2hypersphere


%--------------------------CODE GRAVEYARD----------------------------------
%{
vpaQ = false;
if vpaQ
	precision = 32;
end
waitbarQ = false;

if waitbarQ
	disp([int2str(adjSize(1)) ' neighboring facets']);
	f = waitbar(0,['looping through setup for ' int2str(adjSize(1)) ' neighboring facets']);
end

		if vpaQ
			t(j) = double(ddet(j)/dot(datanorm,nvec{j}));
		else
			t(j) = ddet(j)/dot(datanorm,nvec{j});
		end

		if ~strcmp(invmethod,'extendedCross')
			warning('off',warnID1)
			warning('off',warnID2)
			warning('off',warnID3)
		end

			case 'pinv'
				lambda{j} = (pinv(dmat{j}')*a{j}')';
				posQ(j) = all(lambda{j} > 0) && (sum(lambda{j}) >= 1-tol);


				%c{j} = (vpa(dmat{j})'\vpa(a{j})')'; %barycentric coordinates

		if waitbarQ && (mod(j,100) == 0)
			waitbar(j/adjSize(1),f);
		end


		if posQ(j) == 1
			%test for warnings on the intersecting facet
			
			warning('on',warnID1)
			warning('on',warnID2)
			warning('on',warnID3)
			
			lastwarn('')
			
			(dmat{j}'\a{j}')'; %barycentric coordinates
			
			lstwarn = warning('query','last');
			if strcmp(lstwarn,{warnID1,warnID2,warnID3})
				posQ(j) = 0; %don't consider a point as intersecting a facet if one of the warnings triggers
			end
		end

if waitbarQ
	close(f)
end
%turn warnings back on
warning('on',warnID1)
warning('on',warnID2)
warning('on',warnID3)


if posQtot == 0
	%warning('no intersecting facet detected')
end

	% info about multiple facets
	% 	if sum(posQ) == 2
	% 		%warning('two facets detected')
	% 		disp(['common facets found: ' int2str(sum(posQ))])
	% 	elseif sum(posQ) >= 3
	% 		%warning('3+ facets detected')
	% 		disp(['common facets found: ' int2str(sum(posQ))])
	% 	end
	
	%take max norm
% 	if length(idList) >= 2
% 		%disp('taking datapoint with largest norm')
% 		[~, maxpos] = max(vecnorm(vertcat(a{idList})'));
% 		id = idList(maxpos);
% 	else
%		id = idList;
% 	end



for j = 1:adjSize(1) %loop through facets
	for k = 1:adjSize(2) %loop through vertices of facet
		p{j,k} = meshpts(facetPtIDs(j,k),:);
	end
	
	nmatTemp{j} = vertcat(p{j,:});
	
	dmat{j} = nmatTemp{j};
	for k = 1:adjSize(2)
		nmat{j,k} = nmatTemp{j};
		nmat{j,k}(:,k) = 1;
	end
	for k = 1:adjSize(2) %loop through dimensions
		if vpaQ
			nvec{j}(k) = det(vpa(nmat{j,k},precision));
		else
			nvec{j}(k) = det(nmat{j,k});
		end
	end
	if vpaQ
		ddet(j) = det(vpa(dmat{j},precision));
	else
		ddet(j) = det(dmat{j});
	end
	if waitbarQ && (mod(j,100) == 0)
		waitbar(j/adjSize(1),f);
	end
end
if waitbarQ
	close(f);
end


if waitbarQ
	f = waitbar(0,['looping through projection for ' int2str(adjSize(1)) ' neighboring facets']);
end



	%take max norm
	if length(idList) >= 2
		%disp('taking datapoint with largest norm')
		[~, maxpos] = max(vecnorm(vertcat(a{idList})'));
		id = idList(maxpos);
	else
		id = idList;
	end
	

%}