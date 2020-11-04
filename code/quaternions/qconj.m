function q_star = qconj(q)
% QCONJ takes the conjugate of quaternions
%-------------------------------------------------------------------------%
%Filename:  qconj.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%
% Inputs:
%   q - An npts-by-4 array containing quaternion components in 4-vector 
%       format.
%
% Outputs:
%   q_star - An npts-by-4 array containing the conjugate of each of the 
%            npts quaternions in q, i.e., q_star(i,:) = qconj(q(i,:)).
%-------------------------------------------------------------------------%
%---check inputs---%
assert((size(q,2) == 4),'q must be an npts-by-4 array.')

%---compute conjugates of quaternions---%
q_star = [q(:,1), -q(:,2), -q(:,3), -q(:,4)];