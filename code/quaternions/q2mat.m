% Q2MAT  transforms quaternion components to their matrix representation.
%-------------------------------------------------------------------------%
%Filename:  q2mat.m
%Author:    Oliver Johnson
%Date:      2/23/2011
%
%Updated 6/8/2011 to include isocliny flag option.
%
% Inputs:
%   q - An npts-by-4 array containing quaternion components in 4-vector 
%       format.
%   flag - (optional) A string indicating the isocliny of the 4-by-4
%          rotation matrix representation.  Valid strings are 'right' for a
%          right isoclinic rotation, or 'left' for a left isoclinic
%          rotation. The default is 'left'.
%
% Outputs:
%   Q - A 4-by-4-by-npts array containing the matrix representations of
%       quaternions in q, where Q(:,:,i) is the matrix representation of
%       q(i,:).
%-------------------------------------------------------------------------%

function Q = q2mat(q,flag)

%---check inputs---%
assert((size(q,2) == 4),'q must be an npts-by-4 array.')
if nargin == 1
    flag = 'left';
elseif nargin == 2
    assert(ischar(flag),'flag must be a string.')
    assert(strcmp('left',flag) || strcmp('right',flag),sprintf('Invalid string contained in flag.\nValid strings are:\n(1) ''left''\n(2) ''right'''))
else
    error('Wrong number of input arguments.')
end

%---find number of points---%
npts = size(q,1);

%---perform transformation---%
Q(1,1,:) = reshape(q(:,1),[1 1 npts]);
Q(2,1,:) = reshape(q(:,2),[1 1 npts]);
Q(3,1,:) = reshape(q(:,3),[1 1 npts]);
Q(4,1,:) = reshape(q(:,4),[1 1 npts]);
Q(1,2,:) = reshape(-q(:,2),[1 1 npts]);
Q(2,2,:) = reshape(q(:,1),[1 1 npts]);
Q(3,2,:) = reshape(-q(:,4),[1 1 npts]);
Q(4,2,:) = reshape(q(:,3),[1 1 npts]);
Q(1,3,:) = reshape(-q(:,3),[1 1 npts]);
Q(2,3,:) = reshape(q(:,4),[1 1 npts]);
Q(3,3,:) = reshape(q(:,1),[1 1 npts]);
Q(4,3,:) = reshape(-q(:,2),[1 1 npts]);
Q(1,4,:) = reshape(-q(:,4),[1 1 npts]);
Q(2,4,:) = reshape(-q(:,3),[1 1 npts]);
Q(3,4,:) = reshape(q(:,2),[1 1 npts]);
Q(4,4,:) = reshape(q(:,1),[1 1 npts]);

if strcmp('right',flag)
    Q(2:4,2:4) = Q(2:4,2:4).';
end