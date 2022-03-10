function qm = qlab2qm(qAlab,qBlab,epsijk)
arguments
    qAlab(:,4) double
    qBlab(:,4) double
    epsijk(1,1) double = 1
end
% QLAB2QM  Convert lab/sample frame quaternions of grain A and grain B to misorientation quaternion
% active (epsijk==1) or passive (epsijk==-1) convention
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-08-22
%
% Inputs:
%  qA, qB - Orientations of grains A and B in sample reference frame, resp.
%  convention - francis or johnson convention
%
% Outputs:
%  qm - misorientation quaternion (francis: qB qA*, johnson: qA* qB)
%
% Usage:
%  qm = qlab2qm(qAlab,qBlab,convention)
%
% Dependencies:
%  qmult.m
%  qinv.m
%
% References:
%  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
%
%  [2] Rowenhorst, D.; Rollett, A. D.; Rohrer, G. S.; Groeber, M.; Jackson,
%  M.; Konijnenberg, P. J.; De Graef, M. Consistent Representations of and
%  Conversions between 3D Rotations. Modelling Simul. Mater. Sci. Eng.
%  2015, 23 (8), 083501. https://doi.org/10.1088/0965-0393/23/8/083501.
%--------------------------------------------------------------------------

% qm = qmult(qBlab,qinv(qAlab),epsijk); %issue1: produces consistent results internally within GBdist(), but not during conversion, (passive convention)
qm = qmult(qinv(qAlab),qBlab,epsijk); %issue2: produces consistent results in 5DOF-->octonion conversion, but not within GBdist(), (active convention)