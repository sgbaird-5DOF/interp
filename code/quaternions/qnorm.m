%-------------------------------------------------------------------------%
%Filename:  qnorm.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%
% QNORM returns the magnitude (2-norm) of user supplied quaternions.
%
% Inputs:
%   q - An npts-by-4 array containing quaternion components in 4-vector 
%        format.
%
% Outputs:
%   qn - An npts-by-1 vector containing the norm of each of the npts
%        quaternions in q, i.e., qn(i) = norm(q(i,:)).
%-------------------------------------------------------------------------%

function qn = qnorm(q)

%---check inputs---%
assert((size(q,2) == 4),'q must be an npts-by-4 array.')

%---compute quaternion magnitudes---%
qn = sqrt(sum(q.^2,2));