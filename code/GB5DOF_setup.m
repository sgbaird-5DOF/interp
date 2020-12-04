function propList = GB5DOF_setup(pA,pB,mA,epsijk)
arguments
    pA(:,4)
    pB(:,4)
    mA(:,3)
    epsijk(1,1) double = 1
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

npts = size(pB,1);

% [gA_R,gB_R] = constructGBMatrices(pA,pB,mA,'livermore');

[omA,omB] = deal(zeros(3,3,npts));
for i = 1:npts
    mAtmp = mA(i,:);
    R = vecpair2rmat(mAtmp,[1 0 0]);
    qR = om2qu(R,epsijk);
    qA = qmult(qR,pA,epsijk);
    qB = qmult(qR,pB,epsijk);
    omA(:,:,i) = qu2om(qA,epsijk);
    omB(:,:,i) = qu2om(qB,epsijk);
end

%Calculate GB Energies
element = 'Ni'; %'Al', 'Au', 'Cu'
E(npts) = struct;
E(1).(element) = [];
parfor k = 1:npts
	E(k).(element) = GB5DOF(omA(:,:,k),omB(:,:,k),element);
end
propList = vertcat(E.Ni);
end

%% CODE GRAVEYARD
%{

% f = waitbar(0,['calculating GB energies for ',int2str(nGB),' points.']);	

    %     E.Al(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Al');
	%     E.Au(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Au');
	%     E.Cu(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Cu');

% close(f);
%E.Ni = reshape(E.Ni,size(x)); %x was an output from sphere() in other code

% 	waitbar(k/nGB,f)
%}