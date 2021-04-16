function [ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,method,epsijk,nv)
arguments
    qm %input misorientation quaternions
    nA %input BP normals
    y(:,1) %property values
    qm2 %query misorientations
    nA2 %query BP normals
    K(1,1) double = 1 % number of ensembles
    method char {mustBeMember(method,{'gpr','sphgpr','pbary','sphbary','idw','nn','avg'})} = 'gpr'
    epsijk(1,1) double = 1
    nv.pgnum(1,1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
    nv.databary = [] %for use with bary methods
    nv.facetIDs = [] %for use with bary methods
    nv.ytrue = [] %user-specified "true" values for error calculations
    nv.modelparsspec = struct()
    nv.brkQ(1,1) logical = false %whether to compute BRK values as ytrue
    nv.mygpropts = struct.empty %for use with gpr methods 'gpr' or 'sphgpr'
    nv.r double = [] %for use with 'idw' method, alternatively set to [] for automatic estimation
    nv.uuid(1,8) char = get_uuid() %unique ID associated with this interpolation run
    nv.o = [] %input octonions, specify these or qm/nA pairs
    nv.o2 = [] %query octonions, specify these or qm2/nA2 pairs
    nv.sigma(1,1) double = 0 %synthetic input noise
    nv.dispQ(1,1) logical = false
    nv.KdispQ(1,1) logical = true
    nv.IncludeTies(1,1) {mustBeLogical} = true
    nv.nNN(1,1) double = 1
end
% ENSEMBLEVFZO  take the ensemble average of K VFZs for same set of GBs
KdispQ = nv.KdispQ;
[ypredlist,interpfnlist,mdllist,mdlparslist] = deal(cell(K,1)); %initialize
for k = 1:K
    if KdispQ
        disp(' ')
        disp([int2str(k) '-th ensemble component (ensembleVFZO.m)'])
    end
    %random reference octonion
    nv.oref = get_ocubo();
    nvsub = rmfield(nv,{'KdispQ'});
    nvpairs = namedargs2cell(nvsub);
    %interpolate
    [ypredlist{k},interpfnlist{k},mdllist{k},mdlparslist{k}] = ...
        interp5DOF(qm,nA,y,qm2,nA2,method,epsijk,nvpairs{:});
end
%concatenate
ypredtmp = [ypredlist{:}];
%parse values
% ypred = min(ypredtmp,[],2);
% ypred = median(ypredtmp,2);
ypred = mean(ypredtmp,2);
