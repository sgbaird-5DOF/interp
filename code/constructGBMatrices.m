function [gA_R,gB_R] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab,convention)
arguments
    qA_Lab(:,4)
    qB_Lab(:,4)
    nA_Lab(:,3)
    convention char {mustBeMember(convention,{'olmsted','livermore'})} = 'livermore'
end
%CONSTRUCTGBMATRICES  Make Olmsted GB matrices for GB5DOF using labframe quaternion/normal inputs
% Convert the crystal orientations (expressed as quaternions) of grains
% meeting at a grain boundary (GB) together with the corresponding GB
% normals to the pair of GB matrices defined by Olmsted [1] and required by
% the GB5DOF function [2].
%
%-------------------------------------------------------------------------%
% Filename: constructGBMatrices.m
% Author: Oliver Johnson
% Date: 3/11/2020
%
% NOTE: Olmsted's GB matrices are defined such that the GB plane lies in
% the z = 0 plane of the reference frame. However, the GB5DOF matrices are
% defined such that the GB plane lies in the x = 0 plane of the reference
% frame. Consequently, we provide for calculation of the matrices using
% either convention.
%
% Inputs:
%   qA_Lab - An nGB-by-4 array of quaternions representing the orientations 
%            of one of the grains at the GB (grain A). qA_Lab(i,:) 
%            represents the quaternion defining the orientation of grain A 
%            for the i-th GB. With a reference coordinate system defined 
%            such that z = 0 is the GB plane, grain A is the region z < 0.
%   qB_Lab - An nGB-by-4 array of quaternions representing the orientations 
%            of one of the grains at the GB (grain B). qB_Lab(i,:) 
%            represents the quaternion defining the orientatino of grain B 
%            for the i-th GB. With a reference coordinate system defined 
%            such that z = 0 is the GB plane, grain B is the region z > 0.
%   nA_Lab - A 3-by-nGB array of 3D vectors representing the GB normals
%            expressed in the Lab coordinate system. As a convention, we
%            define nA_Lab such that it is pointing away from grain A (i.e.
%            it is the outward-pointing normal). nA_Lab(:,i) give the
%            vector representing the GB normal for the i-th GB.
%
% Outputs:
%   gA_R - A 3-by-3-by-nGB array of rotation matrices representing the
%          orientation of grain A expressed in the reference coordinate
%          system which has the z = 0 plane aligned with the GB plane.
%          gA_R(:,:,i) represents the orientation of grain A for the i-th
%          GB.
%   gB_R - A 3-by-3-by-nGB array of rotation matrices representing the
%          orientation of grain B expressed in the reference coordinate
%          system which has the z = 0 plane aligned with the GB plane.
%          gB_R(:,:,i) represents the orientation of grain B for the i-th
%          GB.
%
%   NOTE: When using the 'livermore' convention, gA_R is defined such that 
%   gA_R(1,:,i) = nA_A(:,i).' (i.e. the 1st row of gA_R corresponds to the 
%   GB normal pointing away from grain A, expressed in the coordinate 
%   system of grain A, with nA_A(:,i) = gA_Lab(:,:,i).'*nA_Lab(:,i)).
%   Meanwhile, gA_B is defined such that gB_R(1,:,i) = -nB_B(:,i).' (i.e. 
%   the 1st row of gB_R corresponds to the GB normal pointing away from 
%   grain B, expressed in the coordinate system of grain B, with 
%   nB_B(:,i) = gB_Lab(:,:,i).'*nB_Lab(:,i)).
%
%   NOTE: The inputs P and Q for GB5DOF correspond, respectively to gA_R
%   and gB_R. It appears that the convention used in the GB5DOF code is to
%   define a single GB normal, which points away from grain B. This is why
%   the documentation for GB5DOF says that Nq = Q(1,:), whereas by our
%   definition gB_R(1,:,i) = -nB_B(:,i).' (Nq = -Q(1,:)).
%
% [1] Olmsted, D. L. (2009). A new class of metrics for the macroscopic 
%     crystallographic space of grain boundaries. Acta Materialia, 57(9), 
%     2793–2799. https://doi.org/10.1016/j.actamat.2009.02.030
% [2] Bulatov, V. V, Reed, B. W., & Kumar, M. (2014). Grain boundary energy 
%     function for fcc metals. Acta Materialia, 65, 161–175. 
%     https://doi.org/10.1016/j.actamat.2013.10.057
%-------------------------------------------------------------------------%

nA_Lab = nA_Lab.';
% number of GBs
Ngb = size(qA_Lab,1);

% ensure GB normals are normalized
nA_Lab = nA_Lab./sqrt(sum(nA_Lab.^2,1));

% normals pointing away from grain B in the Lab frame
nB_Lab = -nA_Lab;

%% Construct reference frame quaternions (qR)

switch convention
    case 'olmsted' % olmsted convention has z-axis aligned to -nB_Lab
        zR = -nB_Lab;
        
        xR = cross(zR,rand(3,Ngb)); % cross with random vector to get an arbitrary vector in the GB plane to use as the x-axis
        xR = xR./sqrt(sum(xR.^2,1));
    case 'livermore' % livermore convention has x-axis aligned to -nB_Lab
        xR = -nB_Lab;
        
        zR = cross(xR,rand(3,Ngb)); % cross with random vector to get an arbitrary vector in the GB plane to use as the x-axis
        zR = zR./sqrt(sum(zR.^2,1));
end

yR = cross(zR,xR); % cross z with x to define y-axis
yR = yR./sqrt(sum(yR.^2,1));

% make rotation matrices
gR = [permute(xR,[1,3,2]),permute(yR,[1,3,2]),permute(zR,[1,3,2])];

% convert to quaternions
qR = gmat2q(gR);

%% Construct qA_R and qB_R (quaternions in the reference frame)

qA_R = qmultiply(qinv_johnson(qR),qA_Lab);
qB_R = qmultiply(qinv_johnson(qR),qB_Lab);

%% Convert to rotation matrices

gA_R = q2gmat(qA_R);
gB_R = q2gmat(qB_R);

% NOTE: P = gA_R; Q = gB_R;

end

%% testing stuff
% 
% gA_Lab = q2gmat(qA_Lab);
% gB_Lab = q2gmat(qB_Lab);
% 
% nA_A = (gA_Lab.'*nA_Lab).'
% gA_R
% 
% nB_B = (gB_Lab.'*nB_Lab).'
% gB_R
% 
% nB_R = gR.'*nB_Lab;
% 
% figure; % orientation A and orientation R in Lab frame
% quiver3(0,0,0,gA_Lab(1,1),gA_Lab(2,1),gA_Lab(3,1),0,'r'); 
% hold on; 
% quiver3(0,0,0,gA_Lab(1,2),gA_Lab(2,2),gA_Lab(3,2),0,'g'); 
% quiver3(0,0,0,gA_Lab(1,3),gA_Lab(2,3),gA_Lab(3,3),0,'b'); 
% axis equal tight vis3d
% 
% quiver3(0,0,0,gR(1,1),gR(2,1),gR(3,1),0,'r--'); 
% quiver3(0,0,0,gR(1,2),gR(2,2),gR(3,2),0,'g--'); 
% quiver3(0,0,0,gR(1,3),gR(2,3),gR(3,3),0,'b--');
% 
% quiver3(0,0,0,nB_Lab(1),nB_Lab(2),nB_Lab(3),0,'k','linewidth',2);
% 
% figure; % orientation A and orientation R in R frame
% quiver3(0,0,0,gA_R(1,1),gA_R(2,1),gA_R(3,1),0,'r'); 
% hold on; 
% quiver3(0,0,0,gA_R(1,2),gA_R(2,2),gA_R(3,2),0,'g'); 
% quiver3(0,0,0,gA_R(1,3),gA_R(2,3),gA_R(3,3),0,'b'); 
% axis equal tight vis3d
% 
% quiver3(0,0,0,1,0,0,0,'r--'); 
% quiver3(0,0,0,0,1,0,0,'g--'); 
% quiver3(0,0,0,0,0,1,0,'b--');
% 
% quiver3(0,0,0,1,0,0,0,'k','linewidth',2);

%% CODE GRAVEYARD
%{
%% Ensure proper formatting of inputs
assert(size(qA_Lab,2) == 4 && size(qB_Lab,2) == 4,'qA_Lab and qB_Lab must be n-by-4 arrays of quaternions.')
assert(size(nA_Lab,1) == 3,'nA_Lab must be a 3-by-n array of vectors.')
%}