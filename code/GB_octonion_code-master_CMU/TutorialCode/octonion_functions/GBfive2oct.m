function o = GBfive2oct(qmis,nA)
arguments
	qmis(:,4) double {mustBeReal,mustBeFinite}
	nA(:,3) double {mustBeReal,mustBeFinite}
end
%% INPUT DATA 

% qmis: misorientation quaternion
% nA: boundary planes in respective crystal frames of upper or lower grain
% (be consistent with choice across GB's)
% orientations of GB 1 and GB 2
% IMPORTANT CONVENTION: third column of orientation matrix is aligned with boundary plane (z direction)

%% OUTPUT
% o, octonion in GB plane reference frame
npts = size(qmis,1);
assert(npts == size(nA,1),['# quaternions: ' int2str(npts) ', # normals: ' int2str(size(nA,1))]);
Zero = zeros(npts,1);
mA0 = qmult(qmis,qmult([Zero nA],qinv(qmis)));
mA = mA0(:,2:4);

%set mA(:,3) values that are close to -1 or 1 to -1 or 1, respectively
tol = 1e-6;
ids1 = abs(mA(:,3) + 1) < tol; %close to -1 ids
ids2 = abs(mA(:,3) - 1) < tol; %close to 1 ids
mA(ids1,3) = -1;
mA(ids2,3) = 1;

phiA = acos((mA(:,3)));

axisA = normr([mA(:,2) -mA(:,1) Zero]);

pA = ax2qu([axisA phiA]);

qA = qmult(qinv(qmis),pA);
qB = pA;

o = [qA qB];

if ~isreal(o);
	1+1;
end

end


