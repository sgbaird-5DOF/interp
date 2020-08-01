function [Ktr,K,meshpts] = hsphext_subdiv(pts,nint,tricollapseQ)
arguments
	pts double
	nint(1,1) double = 1
	tricollapseQ(1,1) logical = true
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: Hypersphere exterior hull subdivision.
%
% Inputs:
%
% Outputs:
%
% Dependencies:
%	hypersphere_subdiv.m		
%		-facet_subdiv.m
%
%		-tricollapse.m
%
%		-normr.m
%
%	sphere_stereograph.m
%
%	proj_down.m
%
% Notes:
%		Assumes maximum arc length is less than pi (2020-07-20, probably just
%		a bug I can fix) and also less than a hemisphere (?) (since
%		sphere_stereograph.m is used)
%--------------------------------------------------------------------------

% projpts = sphere_stereograph(pts);
projpts = projfacet2hyperplane(mean(pts),pts); %valid for max arc length < pi

projpts = proj_down(projpts);


% get exterior
%joggle input so that coplanar facets aren't merged (because of stereographic projection)
if size(pts,2) < 5
	opts = {'QJ'};
else
	opts = {'Qx','QJ'};
end
K = convhulln(projpts,opts);

%reformat pts and K
IDs = unique(K);
pts = pts(IDs,:);

for i = 1:length(IDs)
	ID = IDs(i);
	K(K==ID) = i;
end

if nint > 1
	% subdivide the "edges"
	
	[Ktr, K, meshpts] = hypersphere_subdiv(pts,K,nint,tricollapseQ);

	
else
	Ktr.main = K;
	%only output exterior points (update: 2020-07-30, realized this didn't preserve indexing, so fixing that)
	Ktr.pts = K;
	meshpts = pts;
end

%------------------------------CODE GRAVEYARD------------------------------
%{
% for i = 1:size(K2,1);
% 	ptIDs = K2(i,:);
% 	extPts = pts(ptIDs,:); %exterior points
% % 	[extprojpts,ext_usv] = proj_down(extPts,1e-6); %projection of exterior facet (i.e. edge)
% 	
% 	delaunayQ = true;
% 	
% 	%subdivide "edge"
% 	[ptstmp,TRItmp] = facet_subdiv(extPts,nint,delaunayQ);
% 		
% 	%renormalize to unit hypersphere	
%  	ptstmp2 = normr(ptstemp);
% 	
% end


% meshpts2 = sphere_stereograph_inverse(meshpts);

% meshpts2 = proj_up(meshpts,usv);
% 
% meshpts3 = sphere_stereograph_inverse(meshpts2);
%}
