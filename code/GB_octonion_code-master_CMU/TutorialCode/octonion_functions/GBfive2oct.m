function o = GBfive2oct(qmis,nA)
%% INPUT DATA 

% qmis: misorientation quaternion
% nA: boundary planes in respective crystal frames of upper or lower grain
% (be consistent with choice across GB's)
% orientations of GB 1 and GB 2
% IMPORTANT CONVENTION: third column of orientation matrix is aligned with boundary plane (z direction)

%% OUTPUT
% o, octonion in GB plane reference frame

mA0 = qmult((qmis),qmult([0 nA],qinv(qmis)));
mA = mA0(2:4);

phiA = acos((mA(3)));

axisA = normr([mA(2) -mA(1) 0]);

pA = ax2qu([axisA phiA]);

qA = qmult(qinv(qmis),pA);
qB = pA;

o = [qA qB];

end


