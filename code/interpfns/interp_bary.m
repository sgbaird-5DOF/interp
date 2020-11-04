function props = interp_bary(mesh,intfacetIDs,qm2,nA2,usv,zeroQ,barytype,barytol,projtol,nnMax,brkQ,NV)
arguments
    mesh(1,1) struct {mustContainFields(mesh,{'pts','ppts','props','sphK'})}
    intfacetIDs cell
    qm2(:,4) double
    nA2(:,3) double
    usv(1,1) struct
    zeroQ(1,1) logical = false
    barytype char {mustBeMember(barytype,{'spherical','planar'})} = 'spherical'
    barytol(1,1) double = 0.2
    projtol(1,1) double = 1e-4
    nnMax(1,1) double = 10
    brkQ(1,1) logical = false
    NV.databary = []
    NV.propList = []
end
% INTERP_BARY  interpolate using spherical or planar barycentric coordinates
pts2 = get_pts(qm2,nA2);
ppts2 = get_ppts(pts2,projtol,usv,zeroQ);
%% when barycentric coordinates aren't specified (new predict points)
if isempty(intfacetIDs)
    intfacetIDs = intersect_facet(mesh.ppts,mesh.sphK,ppts2,inttol,'inttype','planar','nnMax',nnMax);
end

data = struct('pts',pts2,'ppts',ppts2);
if brkQ
    five = struct('q',qm2,'nA',nA2);
    data.props = GB5DOF_setup(five);
else
    data.props = nan(size(pts2));
end

props = get_interp(mesh,data,intfacetIDs,barytype,barytol);

end

%% CODE GRAVEYARD
%{

mdlcmd = @(mesh,data,intfacetIDs,barytol) ...
    get_interp(mesh,data,intfacetIDs,'spherical',barytol);

mdlcmd = @(mesh,data,intfacetIDs,barytol) ...
    get_interp(mesh,data,intfacetIDs,'planar',barytol);

interpfn = @(qm2,nA2) mdlcmd(mesh,...
    struct('pts',get_pts(qm2,nA2),'ppts',get_ppts(qm2,nA2),'props',GB5DOF_setup(GBoct2five(get_pts(qm2,nA2))),'fname','temp.mat'),...
    intersect_facet(mesh.ppts,mesh.sphK,get_ppts(qm2,nA2),inttol,'inttype','planar','nnMax',nnMax),...
    barytol);


interpfn = @(qm2,nA2) get_interp(mesh,...
    ... ---data---
    structcat(rmfield(data,{'pts','ppts'}),...
    struct('pts',get_pts(qm2,nA2),'ppts',get_ppts(qm2,nA2))),...
    ... ^^^data^^^
    intersect_facet(mesh.ppts,mesh.K,get_ppts(qm2,nA2),...
    NV.modelparsspec.inttol,'inttype','planar','nnMax',NV.modelparsspec.nnMax),...
    getinterpmethod,NV.modelparsspec.barytol);


function facetprops = get_facetprops(mesh,intfacetIDs)
arguments
    mesh(1,1) struct {mustContainFields(mesh,{'ppts','sphK','props'})}
    intfacetIDs cell
end

ndatapts = length(intfacetIDs);
facetprops = nan([ndatapts size(mesh.ppts,2)]);
for i = 1:ndatapts
    intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
    vtxIDs = mesh.sphK(intfacetID,:);
    facetprops(i,:) = mesh.props(vtxIDs).'; %properties of vertices of facet
end

end


%     intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
%     vtxIDs = mesh.sphK(intfacetID,:);



% if ~isempty(NV.databary)
%     %% short-circuit barycentric coordinate calculation (same predict points, new mesh property values)
%     %unpack
%     databary = NV.databary;
%     ppts = mesh.ppts;
%     
%     %replace mesh.props with mesh.propList if specified
%     if ~isempty(NV.propList)
%         mesh.props = NV.propList;
%     end
%     
%     %find NaN values & replace with NN values (NN extrapolation)
%     [NNextrapID,~] = isnan(databary);
%     nnList = dsearchn(ppts2(NNextrapID),ppts);
%     d = size(databary,2);
%     
%     facetprops = get_facetprops(ppts,props,facetIDs);
%     
%     % e.g. databary == [NaN NaN NaN], facetprops == [NaN NaN NaN]
%     % --> [1 0 0], [1.213 0 0], such that dot([1 0 0],[1.213 0 0])
%     % == 1.213, where 1.213 is the extrapolated NN value
%     databary(NNextrapID,1) = 1;
%     databary(NNextrapID,2:d) = 0;
%     facetprops(NNextrapID,1) = propList(nnList);
%     facetprops(NNextrapID,2:d) = 0;
%     
%     props = dot(databary,facetprops,2);
%     
% else
% 
% 
% 
% %% HELPER FUNCTION(S)
% function facetprops = get_facetprops(ppts,props,facetIDs)
% arguments
%     ppts double
%     props double
%     facetIDs double
% end
% 
% ndatapts = size(facetIDs,1);
% facetprops = nan([ndatapts size(ppts,2)]);
% for i = 1:ndatapts
%     vtxIDs = facetIDs(i,:);
%     facetprops(i,:) = props(vtxIDs).'; %properties of vertices of facet
% end
% 
% end

%}