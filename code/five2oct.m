function o = five2oct(qm,nA,epsijk)
arguments
    qm(:,4) double {mustBeFinite,mustBeReal}
    nA(:,3) double {mustBeFinite,mustBeReal}
    epsijk(1,1) double = 1
end
% FIVE2OCT Convert 5DOF coordinates (qm,nA) to octonions (o).
% Misorientation quaternion and grain A crystal frame boundary plane normal
% pointing outward from grain A towards grain B to octonion with BP normal
% = [0 0 1];
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-12-05
% 
% Inputs:
%  qm - Misorientation quaternion
%    resp.
%  nA - boundary plane normal (grain A crystal frame) pointing from grain
%  A to grain B
%
% Outputs:
%   o - octonion, with BP normal = [0 0 1]
%
% Usage:
%  o = five2oct(qm,nA)
%
% Dependencies:
%  qinv.m
%  qmA2oct.m
%
% References:
%  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
%  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
%  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

%--------------------------------------------------------------------------

npts = size(qm,1);
% assign identity quaternion to grain A
pA = repmat([1 0 0 0],npts,1);
% assign misorientation quaternion to grain B
pB = qm;
% because pA is identity quaternion, mA == nA
mA = nA;

if epsijk == -1
    %correct for misorientation convention
    pA = qinv(pA);
    pB = qinv(pB);
end

% conversion
o = qmA2oct(pA,pB,nA,epsijk);

end
%----------------------------CODE GRAVEYARD--------------------------------
%{
% %convert to orientation matrices
% omA = qu2om(qA1,epsijk);
% omB = qu2om(qB1,epsijk);

% R = zeros(3,3,npts);
%}