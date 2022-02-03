function out = exec_argfn(fn,S,argoutnames)
arguments
   fn function_handle
   S struct
   argoutnames = []
end
% EXEC_ARGFN execute function_handle by automatically retrieving argument values from a structure
%  and return single output (1 output argument) or struct with fields from
%  argoutnames (multiple output arguments)
%--------------------------------------------------------------------------
% Inputs:
%  fn - function to execute
%
%  S - struct to sample argument values from
%
%  argoutnames (optional) - names of output arguments to use as fields for
%  packaging of output structure
%
% Outputs:
%  out - first output argument from fn if argoutnames is empty. If
%  argoutnames is supplied by the user, then out is a struct containing the
%  output variables as a struct with fieldnames from argoutnames.
%
% Usage:
%  out = exec_argfn(fn,S,argoutnames)
%
% Dependencies:
%  get_argvals.m
%
% Notes:
%  *
%
% Author(s): Sterling Baird
%
% Date: 2020-09-09
%--------------------------------------------------------------------------
%% input/output setup for fn
%input
argvals = get_argvals(S,fn);

%output
if isempty(argoutnames)
    n = 1; %number of output arguments
else
    n = length(argoutnames);
end
argout = cell(1,n); %row vector

%% execute function
[argout{:}] = fn(argvals{:});

%% output
if n == 1
    out = argout{1};
else
    %package output
    out = cell2struct(argout,argoutnames,2);
%     NVpairs = reshape([argoutnames;argout],[],1); %interleave names and values
%     out = struct(NVpairs{:}); %package into struct with supplied names as fields
end

end


%% CODE GRAVEYARD
%{
%     NVpairs = reshape([argoutnames;argout],[],1); %interleave names and values
%     out = struct(NVpairs{:}); %package into struct with supplied names as fields
%}