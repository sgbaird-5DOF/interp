function q = gmat2q(g)
arguments
    g(3,3,:) double
end
% GMAT2Q converts rotation matrices to their equivalent quaternion
% representation.
%-------------------------------------------------------------------------%
%Filename:  gmat2q.m
%Author:    Oliver Johnson
%Date:      3/28/2020
%
% Inputs:
%   g - A 3-by-3-by-npts array of rotation matrices.
%
% Outputs:
%   q - An npts-by-4 array of equiavalent quaternions.
%
% [1] Method from https://doi.org/10.1007/978-3-319-93188-3_5
%
% OLD METHOD (4/1/2013):
% [1] Code adapted from http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
% This algorithm only works when the rotation matrix is stricly skew
% symmetric, in the case that it is symmetric the difference of the off
% diagonal terms is zero so the sign function returns zero. I need to get a
% better method to assign the signs consistently. Take a look at MTEX's
% mat2quat function.
%-------------------------------------------------------------------------%

%---pre-allocate---%
npts = size(g,3);
q = zeros(npts,4);

%---compute transformation---%

%q0
t = g(1,1,:)+g(2,2,:)+g(3,3,:); %traces
t = t(:);
isPos = t > 0; 
q(isPos,1) = 0.5*(sqrt(1+t(isPos)));
q(~isPos,1) = 0.5*(sqrt(...
    squeeze((g(3,2,~isPos)-g(2,3,~isPos)).^2+(g(1,3,~isPos)-g(3,1,~isPos)).^2+(g(2,1,~isPos)-g(1,2,~isPos)).^2)./...
    squeeze(3-t(~isPos))));

%q1
t = g(1,1,:)-g(2,2,:)-g(3,3,:);
t = t(:);
isPos = t > 0; 
q(isPos,2) = 0.5*(sqrt(1+t(isPos)));
q(~isPos,2) = 0.5*(sqrt(...
    squeeze((g(3,2,~isPos)-g(2,3,~isPos)).^2+(g(1,2,~isPos)+g(2,1,~isPos)).^2+(g(3,1,~isPos)+g(1,3,~isPos)).^2)./...
    squeeze(3-t(~isPos))));

%q2
t = -g(1,1,:)+g(2,2,:)-g(3,3,:);
t = t(:);
isPos = t > 0; 
q(isPos,3) = 0.5*(sqrt(1+t(isPos)));
q(~isPos,3) = 0.5*(sqrt(...
    squeeze((g(1,3,~isPos)-g(3,1,~isPos)).^2+(g(1,2,~isPos)+g(2,1,~isPos)).^2+(g(2,3,~isPos)+g(3,2,~isPos)).^2)./...
    squeeze(3-t(~isPos))));

%q3
t = -g(1,1,:)-g(2,2,:)+g(3,3,:);
t = t(:);
isPos = t > 0; 
q(isPos,4) = 0.5*(sqrt(1+t(isPos)));
q(~isPos,4) = 0.5*(sqrt(...
    squeeze((g(2,1,~isPos)-g(1,2,~isPos)).^2+(g(3,1,~isPos)+g(1,3,~isPos)).^2+(g(3,2,~isPos)+g(2,3,~isPos)).^2)./...
    squeeze(3-t(~isPos))));

%signs
q(:,2) = copysign(q(:,2), squeeze(g(3,2,:)-g(2,3,:)));
q(:,3) = copysign(q(:,3), squeeze(g(1,3,:)-g(3,1,:)));
q(:,4) = copysign(q(:,4), squeeze(g(2,1,:)-g(1,2,:)));

% %---compute transformation---%
% t = g(1,1,:)+g(2,2,:)+g(3,3,:); %traces
% t = t(:);
% r = sqrt(1+t);
% q(:,1) = 0.5*r;
% q(:,2) = copysign(0.5*sqrt(1+g(1,1,:)-g(2,2,:)-g(3,3,:)), g(3,2,:)-g(2,3,:));
% q(:,3) = copysign(0.5*sqrt(1-g(1,1,:)+g(2,2,:)-g(3,3,:)), g(1,3,:)-g(3,1,:));
% q(:,4) = copysign(0.5*sqrt(1-g(1,1,:)-g(2,2,:)+g(3,3,:)), g(2,1,:)-g(1,2,:));

end

function x = copysign(x,y)

sgn = ones(size(y));
sgn(y < 0) = -1;

x = sgn.*x; % note we don't need abs(x) because the equations above guarantee x >= 0 to begin with

end