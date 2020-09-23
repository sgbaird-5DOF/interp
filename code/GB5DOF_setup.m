function propList = GB5DOF_setup(five)
% GB5DOF_SETUP  Compute 5DOF GB energy from BRK function
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
qB_Lab = vertcat(five.q);
nGB = size(qB_Lab,1);
qA_Lab = repmat([1 0 0 0],nGB,1);
nA_Lab = vertcat(five.nA).';

[gA_R,gB_R] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab,'livermore');

%Calculate GB Energies
element = 'Ni'; %'Al', 'Au', 'Cu'
E(nGB) = struct;
E(1).(element) = [];
parfor k = 1:nGB
	E(k).(element) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),element);
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