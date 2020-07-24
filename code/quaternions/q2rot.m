%-------------------------------------------------------------------------%
%Filename:  q2rot.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%
% Q2ROT converts quaternion components to their corresponding three 
% independent rotation angles, w (omega), t (theta), p (phi). Q2ROT
% uses the following convention in defining the quaterion representation:
%
% q(1) = cos(w/2);
% q(2) = sin(w/2).*sin(t).*cos(p);
% q(3) = sin(w/2).*sin(t).*sin(p);
% q(4) = sin(w/2).*cos(t);
%
% Inputs:
%   q - An npts-by-4 array containing the user supplied quaternions in 
%       4-vector format.
%
% Outputs:
%   w - An npts-by-1 vector containing hyperspherical angles in the range 
%       [0,pi].
%   t - An npts-by-1 vector containing polar angles in the range [0,pi].
%   p - An npts-by-1 vector containing azimuthal angles in the range 
%       [0,2*pi];
%
% Change log:
%   (1) 7/29/2011 - Changed definition of t to just take the real part.
%       Numerical errors sometimes result in an argument to acos that is
%       greater than 1 (or less than -1), giving an imaginary part. Made 
%       second assertion more strict. Becuase of sqrt in qnorm it sometimes
%       gives a norm of 1 when it is really 1+/-eps, which causes some 
%       errors.
%-------------------------------------------------------------------------%

function [w,t,p] = q2rot(q)

%---check inputs---%
assert((size(q,2) == 4),'q must be an npts-by-4 array.')

%---normalize---%
q = bsxfun(@rdivide,q,qnorm(q));

%---convert negative quaternions to equivalent positive quaternions---%
q(q(:,1) < 0,:) = -q(q(:,1) < 0,:);

%---correct almost zero values---%
q(abs(q) < eps) = 0;

%---compute angles---%
w = 2*acos(q(:,1));
t = real(acos(q(:,4)./sqrt(1-q(:,1).^2)));
p = atan2(q(:,3),q(:,2));
p(p < 0) = p(p < 0)+2*pi;

%---correct results for any identity quaternions---%
isidentity = (w == 0 | w == 2*pi);
if any(isidentity)
    w(isidentity) = 0;
    t(isidentity) = 0;
    p(isidentity) = 0;
end