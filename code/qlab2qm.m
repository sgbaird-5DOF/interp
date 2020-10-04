function qm = qlab2qm(qAlab,qBlab,convention)
arguments
   qAlab(:,4) double
   qBlab(:,4) double
   convention char {mustBeMember(convention,{'francis','johnson'})} = 'johnson'
end
% QLAB2QM  Convert lab/sample frame quaternions of grain A and grain B and
% and compute misorientation quaternion according to Toby Francis's or
% Oliver Johnson's convention
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
%  qinv_francis.m
%
% References:
%  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
%--------------------------------------------------------------------------
switch convention
    case 'francis'
        qm = qmult(qBlab,qinv_francis(qAlab)); %ref [1], eqn. S-4
    case 'johnson'
        qm = qmult(qinv_francis(qAlab),qBlab);
end