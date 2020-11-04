% ROT2Q  convert rotation angles to quaternions
%  converts the three independent rotation angles, w (omega), t
%  (theta), p (phi), to their corresponding quaternion representation. ROT2Q
%  uses the following convention in defining the quaterion representation:
%-------------------------------------------------------------------------%
%Filename:  rot2q.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%
% q(1) = cos(w/2);
% q(2) = sin(w/2).*sin(t).*cos(p);
% q(3) = sin(w/2).*sin(t).*sin(p);
% q(4) = sin(w/2).*cos(t);
%
% Inputs:
%   w - An scalar or npts-by-1 vector containing hyperspherical angles in 
%       the range [0,pi].
%   t - An scalar or npts-by-1 vector containing polar angles in the range 
%       [0,pi].
%   p - An scalar or npts-by-1 vector containing azimuthal angles in the 
%       range [0,2*pi];
%
% Outputs:
%   q - An npts-by-4 array containing the quaternions corresponding to the
%       user supplied angles in 4-vector format.
%-------------------------------------------------------------------------%

function q = rot2q(w,t,p)

%---check inputs---%
assert(iscolumn(w),'w must be an npts-by-1 vector.')
assert(iscolumn(t),'t must be an npts-by-1 vector.')
assert(iscolumn(p),'p must be an npts-by-1 vector.')

npts = max([size(w,1), size(t,1), size(p,1)]);
assert((size(w,1) == npts | isscalar(w)) &...
       (size(t,1) == npts | isscalar(t)) &...
       (size(p,1) == npts | isscalar(p)),'w, t, and p must all be of the same size or scalars.')

%---extend any scalar inputs to match sizes---%
if isscalar(w)
    w = w*ones(npts,1);
end
if isscalar(t)
    t = t*ones(npts,1);
end
if isscalar(p)
    p = p*ones(npts,1);
end
   
%---compute quaternion components---%
q(:,1) = cos(w/2);
q(:,2) = sin(w/2).*sin(t).*cos(p);
q(:,3) = sin(w/2).*sin(t).*sin(p);
q(:,4) = sin(w/2).*cos(t);