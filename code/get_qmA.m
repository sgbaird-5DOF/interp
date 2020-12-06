function [pA,pB,mA] = get_qmA(npts)
arguments
    npts(1,1) double = 1
end
% GET_QMA  generate random lab coordinates (pA,pB,mA)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-12-05
%
% Inputs:
%  npts - number of 5DOF points to generate
%
% Outputs:
%  pA, pB - lab frame quaternions of grains A and B, resp.
%
%  mA - lab frame boundary plane normal
%
% Usage:
%  [pA,pB,mA] = get_qmA(npts);
%
% Dependencies:
%  get_cubo.m
%--------------------------------------------------------------------------
pA = get_cubo(npts);
pB = get_cubo(npts);
mA = normr(rand(npts,3)-0.5);
end