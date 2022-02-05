function o = TJ2oct(EAs,norms,epsijk)
arguments
   EAs(:,3,3) double {mustBeReal,mustBeFinite}
   norms(:,3,3) double {mustBeReal,mustBeFinite}
   epsijk(1,1) double = 1
end
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-08-22
%
% Description: to understand format of B matrix as in B*X = 0
%
% Usage:
%  o = TJ2oct(EAs,norms);
%
% Input:
%  EAs - euler angles in sample frame in in triple junction sets
%
%  norms - boundary plane normals in sample frame in triple junction sets
%
% Output:
%  o - rows of octonions (non-symmetrized)
%
% Dependencies:
%  eu2qu.m
%  eunA2five.m
%  GBfive2oct.m
%
% Notes:
%  GB1 (grain 2 --> grain 3) Boundary normal points from grain 2 towards grain 3
%  GB2 (grain 3 --> grain 1)
%  GB3 (grain 1 --> grain 2)
%
%  example: GB1 of TJ1 is defined by EAs(1,1,:) and norms(1,1,:)
%  example: GB2 of TJ1 is defined by EAs(1,2,:) and norms(1,2,:)
%
%  TJs - triple line direction (sample frame, Cartesian unit vector)
%  e1 - Grain 1 Euler Angles (sample frame, in degrees)
%  e2 - Grain 2 Euler Angles
%  e3 - Grain 3 Euler Angles
%
%  m1 - BP normal from 2 --> 3 (Cartesian unit vector, sample frame)
%  m2 - BP normal from 3 --> 1
%  m3 - BP normal from 1 --> 2
%
% References:
%  [1] Shen, Y. F., Zhong, X., Liu, H., Suter, R. M., Morawiec, A., &
%  Rohrer, G. S. (2019). Determining grain boundary energies from triple
%  junction geometries without discretizing the five-parameter space. Acta
%  Materialia, 166, 126–134. https://doi.org/10.1016/j.actamat.2018.12.022
%
%  [2] https://github.com/sgbaird/TJ2GBE
%
%  [3] https://github.com/Yufeng-shen/TJ2GBE
%
% see also Ni_0131_21520_Cub.mat
%--------------------------------------------------------------------------

[e1,e2,e3,m1,m2,m3,intIDs] = TJ2em(EAs,norms);

q1 = eu2qu(e1,epsijk);
q2 = eu2qu(e2,epsijk);
q3 = eu2qu(e3,epsijk);

%convert euler pairs and BP normal (sample frame) to octonions
o1 = qmA2oct(q2,q3,m1,epsijk); %changed from GBlab2oct to qmA2oct (2020-12-05)
o2 = qmA2oct(q3,q1,m2,epsijk);
o3 = qmA2oct(q1,q2,m3,epsijk);

%interleave so that e.g. o(1,:) corresponds with resE(1)
o(intIDs,:) = [o1; o2; o3];

% %convert from triple junctions to 5DOF
% [qm,nA] = TJ2five(EAs,norms,qmconvention);
% 
% %convert from 5DOF to octonions
% o = GBfive2oct(qm,nA);
