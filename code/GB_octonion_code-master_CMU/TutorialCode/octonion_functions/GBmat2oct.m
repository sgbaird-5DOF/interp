function o = GBmat2oct(OA1,OB1,epsijk)
arguments
    OA1
    OB1
    epsijk(1,1) double = 1
end
%% INPUT DATA 

% OA, OB are normalized orientation matrices in the crystal frames of grain 1 and 2 with miller indices along the COLUMNS. 
% IMPORTANT CONVENTION: third column of orientation matrix is aligned with boundary plane (z direction)

%% OUTPUT

% o, grain boundary octonion expressed in boundary plane reference frame
OA = normr(OA1')'; OB = normr(OB1')';
% disp(OA)
% disp(OB)
%quaternions representing grain orientations in crystal frame, BP along z
pA = om2qu(OA,epsijk); nA = OA(:,3)';
pB = om2qu(OB,epsijk); nB = OB(:,3)';
% disp(pA)
% disp(pB)

%transform BP vectors from crystal frames into sample frame 
% mA0 = qmult(qinv(pA),qmult([0 nA],(pA),epsijk),epsijk); mA = mA0(2:4);
mA = Lpr(qinv(pA),nA,epsijk);
thr = 1e-8;
mA(abs(mA) < thr) = 0;
%mB0 = qmult(qinv(pB),qmult([0 nB],(pB))); mB = mB0(2:4);

%transform BP vector in sample frame to GB frame (aligned with z = [0 0 1])
phiA = acos(mA(3)); axisA = [mA(2) -mA(1) 0];
%phiB = acos(mB(3)); axisB = [mB(2) -mB(1) 0];

%it suffices to only consider rotation of mA into z direction,
%since mA || mB 

Pab = [cos(phiA/2) -epsijk*sin(phiA/2)*axisA];

qA = qmult(pA,Pab,epsijk);
qB = qmult(pB,Pab,epsijk);

o = [qA qB];

end


