function o = om2oct(omA,omB,mA,epsijk)
arguments
    omA(3,3,:) double {mustBeFinite,mustBeReal}
    omB(3,3,:) double {mustBeFinite,mustBeReal}
    mA(:,3) double {mustBeFinite,mustBeReal}
    epsijk(1,1) double = 1
end
% OM2OCT Convert orientation matrices (omA,omB) to octonions (o).
% Orientation matrices of grain A and grain B and sample frame boundary
% plane normal pointing outward from grain A towards grain B to octonion
% with BP normal = [0 0 1]; epsijk == 1 (active rotation matrix), epsijk ==
% -1 (passive rotation matrix)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-12-05
% 
% Inputs:
%  omA, omB - Orientation matrices of grains A and B in sample reference
%  frame, resp.
%
% Outputs:
%   o - octonion, with BP normal = [0 0 1]
%
% Usage:
%  o = om2oct(omA,omB)
%
% Dependencies:
%  om2qu.m
%  qmA2oct.m
%
% References:
%  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
%  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
%  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

%--------------------------------------------------------------------------
pA = om2qu(omA,epsijk);
pB = om2qu(omB,epsijk);

o = qmA2oct(pA,pB,mA,epsijk);

end
%----------------------------CODE GRAVEYARD--------------------------------
%{
%convert to orientation matrices
omA = qu2om(qA1,epsijk);
omB = qu2om(qB1,epsijk);

% rotation matrix to go from mA to [0 0 1]
R = vecpair2rmat(mA(1,:),[0 0 1],1); %not sure why this has to stay in active interpretation for epsijk==1 *and* epsijk==-1

%convert to quaternion
qR = om2qu(R,epsijk);

%apply rotation to pA and pB
qA1 = qmult(qR,pA,epsijk);
qB1 = qmult(qR,pB,epsijk);

o = [qA1 qB1];
%}