function [ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO(qm,nA,y,qm2,nA2,K,method,NV)
arguments
    qm %input misorientation quaternions
    nA %input BP normals
    y(:,1) %property values
    qm2 %query misorientations
    nA2 %query BP normals
    K(1,1) double = 1 % number of ensembles
    method char {mustBeMember(method,{'gpr','sphgpr','pbary','sphbary','idw','nn','avg'})} = 'gpr'
    NV.pgnum(1,1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
    NV.databary = [] %for use with bary methods
    NV.facetIDs = [] %for use with bary methods
    NV.ytrue = [] %user-specified "true" values for error calculations
    NV.modelparsspec = struct()
    NV.brkQ(1,1) logical = true %whether to compute BRK values as ytrue
    NV.mygpropts = struct.empty %for use with gpr methods 'gpr' or 'sphgpr'
    NV.r double = [] %for use with 'idw' method, alternatively set to [] for automatic estimation
    NV.uuid(1,8) char = get_uuid() %unique ID associated with this interpolation run
    NV.o = [] %input octonions, specify these or qm/nA pairs
    NV.o2 = [] %query octonions, specify these or qm2/nA2 pairs
end
[ypredlist,interpfnlist,mdllist,mdlparslist] = deal(cell(K,1)); %initialize
for k = 1:K
    %random reference octonion
    NV.oref = get_ocubo();
    NVpairs = namedargs2cell(NV);
    %interpolate
    [ypredlist{k},interpfnlist{k},mdllist{k},mdlparslist{k}] = ...
        interp5DOF(qm,nA,y,qm2,nA2,method,NVpairs{:});
end
%concatenate
ypredtmp = [ypredlist{:}];
%parse values
% ypred = min(ypredtmp,[],2);
ypred = median(ypredtmp,2);