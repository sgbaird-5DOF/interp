%-------------------------------------------------------------------------%
%Filename:  nle.m
%Author:    Oliver Johnson
%Date:      7/15/2014
%
% NLE performs a numerical less than or equal relational comparison with
% user specified precision.
%
% Inputs:
%   A,B - Two arrays of identical, but arbitrary, size and dimensions.
%   prec - A scalar indicating the number of digits after the decimal point
%          to consider in the comparison.
%
% Outputs:
%   tf - A logical array of the same size and dimensions as A and B, with
%        logical value true wherever A <= B to within the given precision.
%
% Note: nle(x,x-10^(-prec),prec) == false       (greater-than)
%       nle(x,x-10^(-(prec+1)),prec) == true    (equal)
%       nle(x,x+10^(-(prec+1)),prec) == true    (equal)
%       nle(x,x+10^(-prec),prec) == true        (less-than)
%       nle(x,x+10^(-(prec-1)),prec) == true    (less-than)
%-------------------------------------------------------------------------%
function tf = nle(A,B,prec)

ep = roundp(A-B,prec);
tf = ep <= 0;