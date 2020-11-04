function five = get_five(npts)
arguments
    npts(1,1) double = 1
end
% GET_FIVE  generate random cubochoric misorientation and random boundary plane normal pairs
%--------------------------------------------------------------------------
% Inputs:
%  npts - number of 5DOF points to generate
%
% Outputs:
%  five - struct containing at minimum misorientations and boundary plane
%  normals.
%
% Usage:
%  five = get_five(npts);
%
% Dependencies:
%  get_cubo.m
%
% Notes:
%  *
%
% Author(s): Sterling Baird
%
% Date: 2020-10-06
%--------------------------------------------------------------------------
five.q = get_cubo(npts);
five.nA = normr(rand(npts,3)-0.5);
end