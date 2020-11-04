% Q2GMAT converts quaternions to their canonical 3x3 rotation matrix representation
%-------------------------------------------------------------------------%
%Filename:  q2gmat.m
%Author:    Oliver Johnson
%Date:      6/2/2011
%
% Inputs:
%   q - An npts-by-4 matrix of quaternion components.
%
% Outputs:
%   g - A 3-by-3-by-npts array of rotation matrices, where g(:,:,i) gives
%       the rotation matrix for q(i,:).
%-------------------------------------------------------------------------%

function g = q2gmat(q)

assert(size(q,2) == 4,'q must be an npts-by-4 array.')

%---determine the number of quaternions given---%
npts = size(q,1);

%---extract quaternion components---%
q0 = q(:,1);
q1 = q(:,2);
q2 = q(:,3);
q3 = q(:,4);

%---initialize g-matrix array---%
g = zeros(3,3,npts);

%---compute elements of g-matrices---% %consider using instead, the formula from http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
g(1,1,:) = reshape(q0.^2+q1.^2-q2.^2-q3.^2,[1 1 npts]);
g(1,2,:) = reshape(2.*(q1.*q2-q0.*q3),[1 1 npts]);
g(1,3,:) = reshape(2.*(q1.*q3+q0.*q2),[1 1 npts]);
g(2,1,:) = reshape(2.*(q1.*q2+q0.*q3),[1 1 npts]);
g(2,2,:) = reshape(q0.^2-q1.^2+q2.^2-q3.^2,[1 1 npts]);
g(2,3,:) = reshape(2.*(q2.*q3-q0.*q1),[1 1 npts]);
g(3,1,:) = reshape(2.*(q1.*q3-q0.*q2),[1 1 npts]);
g(3,2,:) = reshape(2.*(q2.*q3+q0.*q1),[1 1 npts]);
g(3,3,:) = reshape(q0.^2-q1.^2-q2.^2+q3.^2,[1 1 npts]);