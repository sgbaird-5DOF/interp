%gcfpos
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-08-18
%
% Description: get current figure position and copy to clipboard
%
% Usage:
%  gcfpos
%--------------------------------------------------------------------------
fig = gcf;
pos = fig.Position;
clipboard('copy',pos);