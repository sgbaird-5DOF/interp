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

%-----------------------CODE GRAVEYARD-------------------------------------
%{
% %interleave quaternions (see note)
% qm(id1,:) = qmult(q2,qinv(q3)); % GB1 (grains 2 and 3)
% qm(id2,:) = qmult(q3,qinv(q1)); % GB2 (grains 3 and 1)
% qm(id3,:) = qmult(q1,qinv(q2)); % GB3 (grains 1 and 2)
%
% 
% %rotate boundary plane normals into grain A crystal frame
% Zero = zeros(nGB,1);
% n1tmp = qmult(q1,qmult([Zero m1],qinv(q1)));
% n2tmp = qmult(q2,qmult([Zero m2],qinv(q2)));
% n3tmp = qmult(q2,qmult([Zero m2],qinv(q2)));
% 
% %interleave BP normals (see note)
% nA(id1,:) = n1tmp(:,2:4);
% nA(id2,:) = n2tmp(:,2:4);
% nA(id3,:) = n3tmp(:,2:4);
% 
% %consider packaging above into function, and calling once each for 1, 2,
% %and 3
% 
% %package
% five.q = qm;
% five.nA = nA;

%  I think that it takes a triple junction, and stacks the grain boundaries.
%Question is, how are the grain boundaries defined? Is GB1 grain 1 and
%grain 2 for example?

% M = [];
% for i = 1:nTJ
%     t=TJs(i,:);
%     xls=[0,-t(3),t(2);t(3),0,-t(1);-t(2),t(1),0];
%     M = horzcat(eu2om(EAs(i,2,:)).',eu2om(EAs(i,3,:)).',eu2om(EAs(i,1,:))).';
%     M=dot(xls,M);
% end

% e1 = squeeze(vertcat(EAs(:,1,:)));
% e2 = squeeze(vertcat(EAs(:,2,:)));
% e3 = squeeze(vertcat(EAs(:,3,:)));

% m1 = vertcat(norms(:,1,:));
% m2 = vertcat(norms(:,2,:));
% m3 = vertcat(norms(:,3,:));

% %--initialize
% five = struct();
% five.q = zeros(nGB,4);
% five.nA = zeros(nGB,3);

%interleave quaternions and BP normals
fivetmp.q([id1 id2 id3],:) = vertcat(five1.q, five2.q, five3.q);
fivetmp.nA([id1 id2 id3],:) = vertcat(five1.nA, five2.nA, five3.nA);

% qm = zeros(nGB,4);

% GB1 (i.e. line 1) has grains 2 and 3: EAs(1,2,:) and EAs(1,3,:). Then
% grain 2 == grain A, and grain 3 == grain B since the normal (nA) points
% from grain 2 --> grain 3.

%convert to quaternions
q1 = eu2qu(e1);
q2 = eu2qu(e2);
q3 = eu2qu(e3);


%    qmconvention char {mustBeMember(qmconvention,{'francis','johnson'})} = 'francis'
%}
