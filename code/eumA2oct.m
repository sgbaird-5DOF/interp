function o = eumA2oct(eA,eB,mA,epsijk)
arguments
    eA(:,3) double {mustBeFinite,mustBeReal}
    eB(:,3) double {mustBeFinite,mustBeReal}
    mA(:,3) double {mustBeFinite,mustBeReal}
    epsijk(1,1) double = 1
end
% EUMA2OCT Convert lab coordinates (eA,eB,mA) to octonions (o).
% Sample frame euler angles of grain A and grain B and sample frame
% boundary plane normal pointing outward from grain A towards grain B to
% misorientation quaternion and boundary plane normal (crystal reference
% frame of grain A).
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-08-22
% 
% Inputs:
%  eA, eB - Euler angles of grains A and B in sample reference frame,
%    resp.
%  mA - boundary plane normal (sample reference frame) pointing from grain
%  A to grain B
%
% Outputs:
%   o - octonion, with BP normal = [0 0 1]
%
% Usage:
%  o = eumA2oct(euA,euB,mA)
%
% Dependencies:
%  eu2qu.m
%  qmA2oct.m
%
% References:
%  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
%  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
%  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.


%--------------------------------------------------------------------------
qA = eu2qu(eA,epsijk); % grain A, sample frame
qB = eu2qu(eB,epsijk); % grain B, sample frame

o = qmA2oct(pA,pB,mA,epsijk);

%----------------------------CODE GRAVEYARD--------------------------------
%{
%  [2]
%  Seita, M.; Volpi, M.; Patala, S.; McCue, I.; Schuh, C. A.; Diamanti, M.
%  V.; Erlebacher, J.; Demkowicz, M. J. A High-Throughput Technique for
%  Determining Grain Boundary Character Non-Destructively in
%  Microstructures with through-Thickness Grains. Npj Computational
%  Materials 2016, 2, 16016.
%}