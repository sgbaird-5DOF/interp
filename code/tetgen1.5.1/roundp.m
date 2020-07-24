%-------------------------------------------------------------------------%
%Filename:  roundp.m
%Author:    Oliver Johnson
%Date:      7/16/2014
%
% ROUNDP rounds values in x to n digits past the decimal sign.
%
% Inputs:
%   x - An array of real numbers.
%   n - A non-negative scalar integer specifying the desired precision.
%
% Outputs:
%   X - An array of size(x) giving the values in x rounded to n digits
%       past the decimal sign.
%
% Example:
%   Round the value 8.3498 to 2 digits past the decimal sign.
%   x = 8.3498;
%   n = 2;
%   X = roundp(x,n)
%   X = 8.35   
%-------------------------------------------------------------------------%

function X = roundp(x,n)

prec = (10^n);
X = round(prec*x)/prec;