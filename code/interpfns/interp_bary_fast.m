function [dataprops,facetprops,NNextrapID,nnList] = ...
    interp_bary_fast(ppts,ppts2,meshprops,databary,facetIDs)
arguments
    ppts double
    ppts2 double
    meshprops double
    databary double
    facetIDs double
end
<<<<<<< HEAD
% INTERP_BARY_FAST  short-circuit barycentric interpolation (same input/prediction points, new mesh property values)
% Cannot be used with new prediction points.
=======
% INTERP_BARY_FAST  short-circuit barycentric interpolation (same input/query points, new mesh property values)
% Cannot be used with new predict points.
>>>>>>> master

%find NaN values & replace with NN values (NN extrapolation)
[NNextrapID,~] = isnan(databary);
nnList = dsearchn(ppts2(NNextrapID),ppts);
d = size(databary,2);

%properties of each vertex of each intersecting facet
facetprops = get_facetprops(ppts,meshprops,facetIDs);

% adjust NaN values to do NN extrapolation
% e.g. databary == [NaN NaN NaN], facetprops == [NaN NaN NaN]
% --> [1 0 0], [1.213 0 0], such that dot([1 0 0],[1.213 0 0])
% == 1.213, where 1.213 is the extrapolated NN value
databary(NNextrapID,1) = 1;
databary(NNextrapID,2:d) = 0;
facetprops(NNextrapID,1) = meshprops(nnList);
facetprops(NNextrapID,2:d) = 0;

dataprops = dot(databary,facetprops,2);

end

%% HELPER FUNCTION(S)
function facetprops = get_facetprops(ppts,props,facetIDs)
arguments
    ppts double
    props double
    facetIDs double
end

ndatapts = size(facetIDs,1);
facetprops = nan([ndatapts size(ppts,2)]);
for i = 1:ndatapts
    vtxIDs = facetIDs(i,:);
    facetprops(i,:) = props(vtxIDs).'; %properties of vertices of facet
end

end

%% CODE GRAVEYARD
%{

%unpack
% databary = NV.databary;
% ppts = mesh.ppts;


%replace mesh.props with mesh.propList if specified
if ~isempty(NV.propList)
    mesh.props = NV.propList;
end

% pts2 = get_pts(qm2,nA2);
% ppts2 = get_ppts(pts2,projtol,usv,zeroQ);

%}