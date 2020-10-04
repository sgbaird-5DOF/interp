function o = GBfive2oct(qmis,nA,qmconvention)
arguments
    qmis(:,4) double {mustBeReal,mustBeFinite}
    nA(:,3) double {mustBeReal,mustBeFinite}
    qmconvention char {mustBeMember(qmconvention,{'francis','johnson'})} = 'francis'
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

mA = qmA2nA(qmis,nA,qmconvention);

%set mA(:,3) values that are close to -1 or 1 to -1 or 1, respectively
tol = 1e-4; %reduced tolerance 2020-08-22 b.c. a single value was complex (mA(:,3) had a value greater than 1)
ids1 = abs(mA(:,3) + 1) < tol; %close to -1 ids
ids2 = abs(mA(:,3) - 1) < tol; %close to 1 ids
mA(ids1,3) = -1;
mA(ids2,3) = 1;

phiA = acos((mA(:,3)));

axisA = normr([mA(:,2) -mA(:,1) Zero]);

pA = ax2qu([axisA phiA]);

qA = qlab2qm(qmis,pA,qmconvention);
qB = pA;

o = [qA qB];


%% Code Graveyard
%{

% qA = qmult(qinv_francis(qmis),pA);


% npts = size(qmis,1);
% qA = repmat([1 0 0 0],npts,1); %identity octonion
% qB = qmis;
% o = GBlab2oct(qA,qB,nA,'francis');
% end


% mA0 = qmult(qmis,qmult([Zero nA],qinv(qmis)));
% mA0 = qmult(qmis,[0 Lpr(qmis,nA)]); %this attempt still produce different mA values

% 2020-09-22, switched with wrapper for GBlab2oct.m
% %% INPUT DATA
%
% % qmis: misorientation quaternion
% % nA: boundary planes in respective crystal frames of upper or lower grain
% % (be consistent with choice across GB's)
% % orientations of GB 1 and GB 2
% % IMPORTANT CONVENTION: third column of orientation matrix is aligned with boundary plane (z direction)
%
% %% OUTPUT
% % o, octonion in GB plane reference frame
% npts = size(qmis,1);
% assert(npts == size(nA,1),['# quaternions: ' int2str(npts) ', # normals: ' int2str(size(nA,1))]);
% Zero = zeros(npts,1);
%
% mA0 = qmult(qmis,qmult([Zero nA],qinv(qmis)));
% % mA0 = qmult(qmis,[0 Lpr(qmis,nA)]); %this attempt still produce different mA values
%
% mA = mA0(:,2:4);
%
% % disp(['mA = ' num2str(mA)])
%
% %set mA(:,3) values that are close to -1 or 1 to -1 or 1, respectively
% tol = 1e-4; %reduced tolerance 2020-08-22 b.c. a single value was complex (mA(:,3) had a value greater than 1)
% ids1 = abs(mA(:,3) + 1) < tol; %close to -1 ids
% ids2 = abs(mA(:,3) - 1) < tol; %close to 1 ids
% mA(ids1,3) = -1;
% mA(ids2,3) = 1;
%
% phiA = acos((mA(:,3)));
%
% axisA = normr([mA(:,2) -mA(:,1) Zero]);
%
% pA = ax2qu([axisA phiA]);
%
% qA = qmult(qinv(qmis),pA);
% qB = pA;
%
% o = [qA qB];
%
% % if ~isreal(o)
% % 	1+1;
% % end
%}