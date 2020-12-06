function o = qmA2oct(pA,pB,mA,epsijk)
arguments
    pA(:,4) double {mustBeFinite,mustBeReal}
    pB(:,4) double {mustBeFinite,mustBeReal}
    mA(:,3) double {mustBeFinite,mustBeReal}
    epsijk(1,1) double = 1
end
% QMA2OCT Convert lab coordinates (pA,pB,mA) to octonions (o).
% Sample frame quaternions of grain A and grain B and sample frame boundary
% plane normal pointing outward from grain A towards grain B to octonion
% with BP normal = [0 0 1];
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-12-05
% 
% Inputs:
%  pA, pB - Quaternions grains A and B in sample reference frame,
%    resp.
%  mA - boundary plane normal (sample reference frame) pointing from grain
%  A to grain B
%
% Outputs:
%   o - octonion, with BP normal = [0 0 1]
%
% Usage:
%  o = qmA2oct(pA,pB,mA)
%
% Dependencies:
%  vecpair2rmat.m
%  om2qu.m
%  qu2om.m
%  qmult.m
%  qinv.m
%
% References:
%  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
%  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
%  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

%--------------------------------------------------------------------------
npts = size(pB,1);
qR = zeros(npts,4);
for i = 1:npts
    mAtmp = mA(i,:);
    %rotation matrix to go from mAtmp to [0 0 1]
    R = vecpair2rmat(mAtmp,[0 0 1],1); %not sure why this has to stay in active interpretation for epsijk==1 *and* epsijk==-1
    %convert to quaternion
    qR(i,:) = om2qu(R,epsijk);
end

%apply rotation to pA and pB
qA1 = qmult(qR,pA,epsijk);
qB1 = qmult(qR,pB,epsijk);

o = [qA1 qB1];

end
%----------------------------CODE GRAVEYARD--------------------------------
%{
% %convert to orientation matrices
% omA = qu2om(qA1,epsijk);
% omB = qu2om(qB1,epsijk);

% R = zeros(3,3,npts);

if isempty(pA)
    pA = repmat([1 0 0 0],npts,1);
end
if epsijk == -1
    pA = qinv(pA);
    pB = qinv(pB);
end
%}