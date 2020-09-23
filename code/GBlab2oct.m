%-------------------------------------------------------------------------%
% Filename: GBlab2oct.m
% Author: Oliver Johnson
% Date: 9/22/2020
%
% GBLAB2OCT will convert the crystal orientations (expressed as 
% quaternions) of grains meeting at a grain boundary (GB) together with the
% corresponding GB normals to the GB octonion defined by Francis [1] and 
% required by the GBdist function [2].
%
% NOTE: The convention used by [1] and required for proper use of the 
% GBdist function [2] is that the GB plane reference frame is defined such 
% that the GB plane lies in the zR = 0 plane, with grain A occupying zR > 0
% and grain B occupying zR < 0, which implies that the normal pointing away 
% from grain B (which we call nB, using the outward-pointing normal 
% convention) is parallel to zR. Note that this is in contrast to Olmsted's 
% convention which is that the GB plane lies in the zR = 0 plane with grain
% A occupying zR < 0 and grain B occupying zR > 0, which implies that the
% normal pointing away from grain A (which we call nA, using the
% outward-pointing normal convention) is parallel to zR. Moreover, the 
% GB5DOF function has a completely different convention, which is that the 
% GB plane reference frame is defined such that the GB plane lies in the 
% xR = 0 plane with grain A occupying xR < 0 and grain B occupying xR > 0,
% which implies that the normal pointing away from grain A (which we call
% nA, using the outward-pointing normal convention) is parallel to xR.
%
% Inputs:
%   qA_Lab - An nGB-by-4 array of quaternions representing the orientations 
%            of one of the grains at the GB (grain A). qA_Lab(i,:) 
%            represents the quaternion defining the orientation of grain A 
%            for the i-th GB. With a reference coordinate system defined 
%            such that z = 0 is the GB plane, grain A is the region z < 0.
%   qB_Lab - An nGB-by-4 array of quaternions representing the orientations 
%            of one of the grains at the GB (grain B). qB_Lab(i,:) 
%            represents the quaternion defining the orientation of grain B 
%            for the i-th GB. With a reference coordinate system defined 
%            such that z = 0 is the GB plane, grain B is the region z > 0.
%   nA_Lab - A 3-by-nGB array of 3D vectors representing the GB normals
%            expressed in the Lab coordinate system. As a convention, we
%            define nA_Lab such that it is pointing away from grain A (i.e.
%            it is the outward-pointing normal). nA_Lab(:,i) give the
%            vector representing the GB normal for the i-th GB.
%
% Outputs:
%   oAB - An nGB-by-8 array of octonions where oAB(i,:) represents the i-th
%         GB. Note that these are NOT normalized octonions, rather they are
%         composed of two normalized quaternions, so that the magnitude of
%         each octonion is sqrt(2). This is the expected format of
%         octonions used as input to GBdist.
%
% [1] T. Francis, I. Chesser, S. Singh, E.A. Holm, M. De Graef, A geodesic 
%     octonion metric for grain boundaries, Acta Mater. 166 (2019) 135–147.
%     doi:10.1016/j.actamat.2018.12.034.
% [2] https://github.com/ichesser/GB_octonion_code
%-------------------------------------------------------------------------%

function oAB = GBlab2oct(qA_Lab,qB_Lab,nA_Lab,convention)

%% Ensure proper formatting of inputs

assert(size(qA_Lab,2) == 4 && size(qB_Lab,2) == 4,'qA_Lab and qB_Lab must be n-by-4 arrays of quaternions.')
assert(size(nA_Lab,1) == 3,'nA_Lab must be a 3-by-n array of vectors.')

% number of GBs
Ngb = size(qA_Lab,1);

% ensure GB normals are normalized
nA_Lab = nA_Lab./sqrt(sum(nA_Lab.^2,1));

% normals pointing away from grain B in the Lab frame
nB_Lab = -nA_Lab;

%% Construct reference frame quaternions (qR)

switch convention
    case 'francis' % francis convention has z-axis aligned to nB_Lab
        zR = nB_Lab;
        
        xR = cross(zR,rand(3,Ngb)); % cross with random vector to get an arbitrary vector in the GB plane to use as the x-axis
        xR = xR./sqrt(sum(xR.^2,1));
%     case 'olmsted' % olmsted convention has z-axis aligned to -nB_Lab
%         zR = -nB_Lab;
%         
%         xR = cross(zR,rand(3,Ngb)); % cross with random vector to get an arbitrary vector in the GB plane to use as the x-axis
%         xR = xR./sqrt(sum(xR.^2,1));
%     case 'livermore' % livermore convention has x-axis aligned to -nB_Lab
%         xR = -nB_Lab;
%         
%         zR = cross(xR,rand(3,Ngb)); % cross with random vector to get an arbitrary vector in the GB plane to use as the x-axis
%         zR = zR./sqrt(sum(zR.^2,1));
end

yR = cross(zR,xR); % cross z with x to define y-axis
yR = yR./sqrt(sum(yR.^2,1));

% make rotation matrices
gR = [permute(xR,[1,3,2]),permute(yR,[1,3,2]),permute(zR,[1,3,2])];

% convert to quaternions
qR = gmat2q(gR);

%% Construct qA_R and qB_R (quaternions in the GB reference frame)

qA_R = qmultiply(qA_Lab,qinv(qR)); % note this is using Toby's misorientation convention, these should not be used with our code directly
qB_R = qmultiply(qB_Lab,qinv(qR));

%% Convert to octonion

oAB = [qA_R, qB_R];

end