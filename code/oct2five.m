function [qm,nA] = oct2five(o,epsijk)
arguments
    o(:,8) double
    epsijk(1,1) double = 1
end
% OCT2FIVE Convert octonions (o) to 5DOF coordinates (qm,nA).
% octonion with BP normal = [0 0 1]. Misorientation quaternion and grain A
% crystal frame boundary plane normal pointing outward from grain A towards grain B.
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-01-27
% 
% Inputs:
%   o - octonion, with BP normal = [0 0 1]
%
% Outputs:
%  qm - Misorientation quaternion
%  nA - boundary plane normal (grain A crystal frame) pointing from grain
%  A to grain B
%
% Usage:
%  [qm,nA] = oct2five(o)
%
% Dependencies:
%  qinv.m
%  qmA2five.m
%
% References:
%  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
%  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
%  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.
%--------------------------------------------------------------------------
npts = size(o,1);
qA = normr(o(:,1:4));
qB = normr(o(:,5:8));
mA = repmat([0 0 1],npts,1);
[qm,nA] = qmA2five(qA,qB,mA,epsijk);