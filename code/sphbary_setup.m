function [c,intfacetIDs] = sphbary_setup(pts,facetPtIDs,data,tol)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
% 
% Inputs:
%
% Outputs:
%
% Dependencies:
%
%--------------------------------------------------------------------------

%initialize
nfacets = size(facetPtIDs,1);
posQ = false(nfacets,1);
c = zeros(nfacets,size(data,2));

%disable warnings
warnID1 = 'MATLAB:nearlySingularMatrix';
warnID2 = 'MATLAB:singularMatrix';
warnID3 = 'symbolic:mldivide:InconsistentSystem';

warning('off',warnID1)
warning('off',warnID2)
warning('off',warnID3)

%loop through facets
for i = 1:nfacets
	idtemp = facetPtIDs(i,:);
	vertices = pts(idtemp,:);
	[c(i,:),Pnew] = sphbary(data,vertices);
	
	%check if spherical barycentric coordinate reqs are met
	posbaryQ = (sum(c(i,:)) >= 1-tol) && all(c(i,:) >= -tol);
	if posbaryQ
		
		%enable warnings
		warning('on',warnID1)
		warning('on',warnID2)
		warning('on',warnID3)
		
		sphbary(data,vertices,Pnew);
		lstwarn = warning('query','last');
		
		if ~strcmp(lstwarn,{warnID1,warnID2,warnID3})
			posQ(i) = true; %position Q, true == intersecting facet
		end
	end
	
end

%enable warnings
warning('on',warnID1)
warning('on',warnID2)
warning('on',warnID3)

%find intersecting facet IDs
if sum(posQ) > 0
	intfacetIDs = find(posQ);
	c = c(posQ,:);
else
	intfacetIDs = [];
	c = [];
end

end