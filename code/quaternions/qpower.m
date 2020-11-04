% QPOWER computes the power of a quaternion, defined as: q^k = q*q*q*... (k times)
%-------------------------------------------------------------------------%
%Filename:  qpower.m
%Author:    Oliver Johnson
%Date:      5/25/2013
%
% Inputs:
%   q - An nq-by-4 array of quaternions. Each row is a quaternion.
%   k - A real scalar defining the exponent.
%
% Outputs:
%   q - An nq-by-4 array of quaternions.
%-------------------------------------------------------------------------%

function q = qpower(q,k)

%---parse inputs---%
assert(isreal(k) && isscalar(k),'k must be a real scalar.')

%---compute power---%
[w,t,p] = q2rot(q);
w = k*w;
q = rot2q(w,t,p);

% Old Way: Recursive definition scales with k
% if k == 0 %define the zero-th power as the identity quaternion
%     nq = size(q,1);
%     q = [ones(nq,1), zeros(nq,3)];
% elseif k == 1 %q^1 = q
%     return
% else %call recursively
%     q = qmultiply(q,qpower(q,k-1));
% end