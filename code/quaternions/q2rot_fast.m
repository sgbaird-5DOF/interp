%-------------------------------------------------------------------------%
%Filename:  q2rot_fast.m
%Author:    Oliver Johnson
%Date:      5/27/2013
%
% Same thing as q2rot, but does not perform the same checks in order to
% achieve a speedup of 3x. As a result the output is not guaranteed to lie
% in the same intervals as q2rot. Instead we have w = [0,2*pi], t =
% [0,pi]?, p = [-pi,pi].
%-------------------------------------------------------------------------%

function [w,t,p] = q2rot_fast(q)

% %---check inputs---%
% assert((size(q,2) == 4),'q must be an npts-by-4 array.')

% %---normalize---%
% q = bsxfun(@rdivide,q,qnorm(q));

% %---convert negative quaternions to equivalent positive quaternions---%
% q(q(:,1) < 0,:) = -q(q(:,1) < 0,:);

% %---correct almost zero values---%
% q(abs(q) < eps) = 0;

%---compute angles---%
w = 2*acos(q(:,1));
t = real(acos(q(:,4)./sqrt(1-q(:,1).^2)));
p = atan2(q(:,3),q(:,2));
% p(p < 0) = p(p < 0)+2*pi;

%---correct results for any identity quaternions---%
isidentity = (w == 0 | w == 2*pi);
if any(isidentity)
    w(isidentity) = 0;
    t(isidentity) = 0;
    p(isidentity) = 0;
end