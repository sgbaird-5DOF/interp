function vq = get_slerp(v0,v1,dt)
arguments
    v0(:,4) double = [1 0 0 0]
    v1(:,4) double = [0 0 0 1]
    dt(1,1) double = 0.1
end
% import numpy as np
% 
% def slerp(v0, v1, t_array):
%     """Spherical linear interpolation."""
%     # >>> slerp([1,0,0,0], [0,0,0,1], np.arange(0, 1, 0.001))
%     t_array = np.array(t_array)
%     v0 = np.array(v0)
%     v1 = np.array(v1)
%     dot = np.sum(v0 * v1)
%
%spaced t values
tarray = 0:dt:1;
npts = length(tarray);

%     Only unit quaternions are valid rotations.
%     Normalize to avoid undefined behavior.
v0 = normr(v0);
v1 = normr(v1);

% Compute the cosine of the angle between the two vectors.
x = dot(v0,v1,2);
if x < 0
%      If the dot product is negative, slerp won't take
%     the shorter path. Note that v1 and -v1 are equivalent when
%     the negation is applied to all four components. Fix by 
%     reversing one quaternion.
    v1 = -v1;
    x = -x;
end
xthresh = 0.9995;
if x > xthresh
         % If the inputs are too close for comfort, linearly interpolate
        % and normalize the result.
    vq = v0+tarray*(v1-v0);

else
%     Since dot is in range [0, xthresh], acos is safe
    theta0 = acos(x); % theta_0 = angle between input vectors
    stheta0 = sin(theta0); % theta = angle between v0 and result
    theta = theta0 .* tarray; % 
    stheta = sin(theta);
    s0 = cos(theta) - x .* stheta ./ stheta0; % == sin(theta_0 - theta) / sin(theta_0)
    s1 = stheta ./ stheta0;
    
    v0rep = repmat(v0,npts,1);
    v1rep = repmat(v1,npts,1);
    s0rep = repmat(s0.',1,4);
    s1rep = repmat(s1.',1,4);
    vq = v0rep.*s0rep+v1rep.*s1rep;
end

%     if dot < 0.0:
%         v1 = -v1
%         dot = -dot
%     
%     DOT_THRESHOLD = 0.9995
%     if dot > DOT_THRESHOLD:
%         result = v0[np.newaxis,:] + t_array[:,np.newaxis] * (v1 - v0)[np.newaxis,:]
%         return (result.T / np.linalg.norm(result, axis=1)).T
%     
%     theta_0 = np.arccos(dot)
%     sin_theta_0 = np.sin(theta_0)
% 
%     theta = theta_0 * t_array
%     sin_theta = np.sin(theta)
%     
%     s0 = np.cos(theta) - dot * sin_theta / sin_theta_0
%     s1 = sin_theta / sin_theta_0
%     return (s0[:,np.newaxis] * v0[np.newaxis,:]) + (s1[:,np.newaxis] * v1[np.newaxis,:])