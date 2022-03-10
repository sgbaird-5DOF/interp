function propList = GB5DOF_setup(pA,pB,mA,mat,epsijk,nv)
arguments
    pA(:,4) = []
    pB(:,4) = []
    mA(:,3) = [0 0 1] %default is octonion convention
    mat char = 'Ni'
    epsijk(1,1) double = 1
    nv.o double = []
end
%GB5DOF_SETUP  Compute 5DOF GB energy from BRK function
%--------------------------------------------------------------------------
% Author(s): Oliver Johnson, Sterling Baird
% Date: 2020-07-27
% 
% Inputs:
%  five - struct containing at minimum misorientation quaternions (q) and
%  boundary plane normals in crystal reference frame of grain A (nA)
%  pointing towards grain B.
%
% Outputs:
%  propList - grain boundary energies computed at each grain boundary in
%  five
%
% Usage:
%  propList = GB5DOF_setup(five)
%
% Dependencies:
%  constructGBMatrices.m
%  GB5DOF.m
%
% Notes:
%  Find a good verification/ground truth (i.e. here is a GBE for a GB with
%  this representation.)
%--------------------------------------------------------------------------

% Compute GB matrices
% if ~isempty(five)
%     pB = vertcat(five.q);
%     pA = repmat([1 0 0 0],nGB,1);
%     mA = vertcat(five.nA).';
% end

o = nv.o;
if isempty([pA;pB])
    assert(~isempty(o),'specify either pA and pB or o')
    pA = o(:,1:4);
    pB = o(:,5:8);
end
pA = normr(pA);
pB = normr(pB);

npts = size(pB,1);

if isempty(mA)
    mA = [0 0 1];
end

if size(mA,1) == 1
    mA = repmat(mA,npts,1);
end

if isempty(pA)
    % case that pB is actually misorientation
    pA = repmat([1 0 0 0],npts,1);
end

if epsijk == -1
    pA = qinv(pA);
    pB = qinv(pB);
end

% [gA_R,gB_R] = constructGBMatrices(pA,pB,mA,'livermore');

[omA,omB] = deal(zeros(3,3,npts));
parfor i = 1:npts
    mAtmp = mA(i,:);
    R = vecpair2rmat(mAtmp,[1 0 0],1);
    qR = om2qu(R,1);
    pAtmp = pA(i,:);
    pBtmp = pB(i,:);
    qA = qmult(qR,pAtmp,1);
    qB = qmult(qR,pBtmp,1);
    omA(:,:,i) = qu2om(qA,1);
    omB(:,:,i) = qu2om(qB,1);
end

%Calculate GB Energies
% mat = 'Ni'; %'Al', 'Au', 'Cu'
E(npts) = struct;
E(1).(mat) = [];
parfor k = 1:npts
	E(k).(mat) = GB5DOF(omA(:,:,k),omB(:,:,k),mat);
end
propList = vertcat(E.(mat));
end