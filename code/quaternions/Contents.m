% QUATERNIONS
%
% Files
%   cell600          - generates the vertices of a 600-cell (the 4D analog of the
%   gmat2q           - converts rotation matrices to their equivalent quaternion
%   newquaternion    - Construct quaternion 
%   project4D        - uses an isovolumetric projection to project points on the unit hypersphere into the unit ball in 3D
%   q2gmat           - converts quaternions to their canonical 3x3 rotation matrix representation
%   q2hsv            - quaternion to continuous, one-to-one coloring scheme (Patala2011)
%   q2mat            - transforms quaternion components to their matrix representation.
%   q2rgb            - convert quaternion to RGB values, adapted from colormap432.m of MTEX 4.5.0
%   q2rod            - convert quaternion to rodriguez vector
%   q2rot            - convert quaternions to rotation angles.
%   q2rot_fast       - Same as q2rot, but does not perform the same checks (3x speedup)
%   qconj            - takes the conjugate of quaternions
%   qinv_johnson     - QINV takes the inverse of quaternions.
%   qmultiply        - performs quaternion multiplication.
%   qnorm            - returns the magnitude (2-norm) of user supplied quaternions.
%   qpower           - computes the power of a quaternion, defined as: q^k = q*q*q*... (k times)
%   quaternionrotate - Rotates a 3D vector by a quaternion 
%   r2q              - convert r (Euclidean BP coordinates?) to quaternion
%   randq            - generates uniformly distributed positive quaternions.
%   randq_validate   - RANDQVALIDATE  validate randq.m
%   rod2q            - convert Rodrigues vector to quaternion
%   rod2rgb          - convert from Rodrigues vectors to RGB values, adapted from colormap432.m of MTEX 4.5.0
%   rot2q            - convert rotation angles to quaternions





