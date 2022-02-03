function [argnames,argvals] = get_args(S,fn)
arguments
    S struct
    fn function_handle
end
% GET_ARGS get argument names and values for a function from a structure
%--------------------------------------------------------------------------
% Inputs:
%  S - struct array containing e.g. parameter combinations with
%  fieldnames that match fn
%
%  fn - function handle, e.g. @(a,b,c) a + b + c, where a,b,c are the
%  arguments to extract names and values for*
%
% Outputs:
%  argnames - names of arguments to fn
%
%  argvals - values of arguments from S
%
% Usage:
%  [argnames,argvals] = get_args(S,fn)
%
% Dependencies:
%  get_argnames.m
%  get_argvals.m
%
% Notes:
%  the arguments to the function handle need to have matches in the fields
%  to parcombs.
%
% Author(s): Sterling Baird
%
% Date: 2020-09-09
%--------------------------------------------------------------------------

argnames = get_argnames(fn); %argument names
argvals = get_argvals(S,fn,argnames); %argument values
end %get_args