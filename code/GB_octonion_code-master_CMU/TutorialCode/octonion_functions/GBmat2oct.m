function o = GBmat2oct(OA1,OB1)
%% INPUT DATA 

% OA, OB are normalized orientation matrices in the crystal frames of grain 1 and 2 with miller indices along the COLUMNS. 
% IMPORTANT CONVENTION: third column of orientation matrix is aligned with boundary plane (z direction)

%% OUTPUT

% o, grain boundary octonion expressed in boundary plane reference frame
OA = normr(OA1')'; OB = normr(OB1')';
% disp(OA)
% disp(OB)
%quaternions representing grain orientations in crystal frame, BP along z
pA = om2qu(OA); nA = OA(:,3)';
pB = om2qu(OB); nB = OB(:,3)';
% disp(pA)
% disp(pB)

%transform BP vectors from crystal frames into sample frame 
mA0 = qmult(qinv(pA),qmult([0 nA],(pA))); mA = mA0(2:4);
%mB0 = qmult(qinv(pB),qmult([0 nB],(pB))); mB = mB0(2:4);

%transform BP vector in sample frame to GB frame (aligned with z = [0 0 1])
phiA = acos(mA(3)); axisA = [mA(2) -mA(1) 0];
%phiB = acos(mB(3)); axisB = [mB(2) -mB(1) 0];

%it suffices to only consider rotation of mA into z direction,
%since mA || mB 

Pab = [cos(phiA/2) sin(phiA/2)*axisA];

qA = qmult(pA,Pab);
qB = qmult(pB,Pab);

o = [qA qB];

end


