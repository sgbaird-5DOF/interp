%GCFPOS  get current figure position and copy to clipboard
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-08-18
%
% Usage:
%  gcfpos()
%--------------------------------------------------------------------------
fig = gcf;
pos = fig.Position;
clipboard('copy',pos);