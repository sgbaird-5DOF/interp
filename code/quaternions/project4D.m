% PROJECT4D uses an isovolumetric projection to project points on the unit hypersphere into the unit ball in 3D
%-------------------------------------------------------------------------%
%Filename:  project4D.m
%Author:    Oliver Johnson
%Date:      6/21/2011
%
% Inputs:
%   w - An array giving the rotation angles in the range [0,pi]
%   t - An array giving the polar angles, measured from the +z axis and in 
%       the range [0,pi].
%   p - An array giving the azimuthal angles, measured from the +x axis and
%       in the range [0,2*pi].
%
% Outputs:
%   x,y,z - Arrays of size(w) containing the cartesian coordinates of the
%           projected points.
%
% [1] Mason, J, and C Schuh. “Hyperspherical harmonics for the
%     representation of crystallographic texture.” Acta Materialia 56.20
%     (2008): 6141-6155.
%-------------------------------------------------------------------------%

function [x,y,z] = project4D(w,t,p,projection)

assert(all(w(:) >= 0) && all(w(:) <= pi),'w must be in the range [0,pi].')
assert(all(t(:) >= 0) && all(t(:) <= pi),'t must be in the range [0,pi].')
assert(all(p(:) >= 0) && all(p(:) <= 2*pi),'p must be in the range [0,2*pi].')

if nargin < 4 || strcmpi(projection,'isovolumetric')
    %---use isovolumetric projection to 3D---%
    r = (3/2)^(1/3).*(w/2-sin(w/2).*cos(w/2)).^(1/3);
    x = sin(t).*cos(p).*r;
    y = sin(t).*sin(p).*r;
    z = cos(t).*r;
elseif nargin == 4 && strcmpi(projection,'rodrigues')
    %---use Rodriguez projection---%
    r = tan(w/2);
    x = sin(t).*cos(p).*r;
    y = sin(t).*sin(p).*r;
    z = cos(t).*r;
end