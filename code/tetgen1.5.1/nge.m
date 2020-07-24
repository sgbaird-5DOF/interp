%-------------------------------------------------------------------------%
%Filename:  nge.m
%Author:    Oliver Johnson
%Date:      7/15/2014
%
% NGE performs a numerical greater than or equal relational comparison with
% user specified precision.
%
% Inputs:
%   A,B - Two arrays of identical, but arbitrary, size and dimensions.
%   prec - A scalar indicating the number of digits after the decimal point
%          to consider in the comparison.
%
% Outputs:
%   tf - A logical array of the same size and dimensions as A and B, with
%        logical value true wherever A >= B to within the given precision.
%
% Note: nge(x,x-10^(-prec),prec) == true        (greater-than)
%       nge(x,x-10^(-(prec+1)),prec) == true    (equal)
%       nge(x,x+10^(-(prec+1)),prec) == true    (equal)
%       nge(x,x+10^(-prec),prec) == false       (less-than)
%       nge(x,x+10^(-(prec-1)),prec) == false   (less-than)
%-------------------------------------------------------------------------%
function tf = nge(A,B,prec)

ep = roundp(A-B,prec);
tf = ep >= 0;