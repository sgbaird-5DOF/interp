%-------------------------------------------------------------------------%
%Filename:  nneq.m
%Author:    Oliver Johnson
%Date:      7/15/2014
%
% NNEQ performs a numerical not equal to relational comparison with
% user specified precision.
%
% Inputs:
%   A,B - Two arrays of identical, but arbitrary, size and dimensions.
%   prec - A scalar indicating the number of digits after the decimal point
%          to consider in the comparison.
%
% Outputs:
%   tf - A logical array of the same size and dimensions as A and B, with
%        logical value true wherever A ~= B to within the given precision.
%
% Note: nneq(x,x-10^(-prec),prec) == true        (greater-than)
%       nneq(x,x-10^(-(prec+1)),prec) == false   (equal)
%       nneq(x,x+10^(-(prec+1)),prec) == false   (equal)
%       nneq(x,x+10^(-prec),prec) == true        (less-than)
%       nneq(x,x+10^(-(prec-1)),prec) == true    (less-than)
%-------------------------------------------------------------------------%
function tf = nneq(A,B,prec)

ep = roundp(A-B,prec);
tf = ep ~= 0;