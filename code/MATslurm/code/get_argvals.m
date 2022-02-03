function argvals = get_argvals(S,fn,argnames)
arguments
    S
    fn
    argnames = []
end
% GET_ARGVALS get argument values of a function by parsing values from a struct
%--------------------------------------------------------------------------
% Inputs:
%  S - struct array containing e.g. parameter combinations with
%  fieldnames that match fn
%
%  fn - function handle, e.g. @(a,b,c) a + b + c, where a,b,c are the
%  arguments to extract names and values for*
%
% Outputs:
%  argvals - values of arguments from S
%
% Usage:
%  argvals = get_args(S,fn)
%
% Author(s): Sterling Baird
%
% Date: 2020-09-09
%--------------------------------------------------------------------------
if isempty(argnames)
    argnames = get_argnames(fn); %argument names
end
nargs = length(argnames); %number of arguments

%% argument values
ncombs = length(S);
%initialize
argvals = cell(ncombs,nargs);

%loop through variables
for i = 1:nargs
    %unpack argument name
    argname = argnames{i};
    %concatenate and convert to cell
    if ~isfield(S,argname)
        error(['double check that "pars" variable in randOctParityData contains ',argname])
    end
        argvals(:,i) = num2cell(vertcat(S.(argname)),2);
end
