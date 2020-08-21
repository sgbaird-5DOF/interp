function propList = GB5DOF_setup(five)
%Compute 5DOF GB energy from BRK function

% Compute GB matrices
qB_Lab = vertcat(five.q);
nGB = size(qB_Lab,1);
qA_Lab = repmat([1 0 0 0],nGB,1);
nA_Lab = vertcat(five.nA).';

[gA_R,gB_R] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab,'livermore');

%Calculate GB Energies
% f = waitbar(0,['calculating GB energies for ',int2str(nGB),' points.']);
E(nGB) = struct;
E(1).Ni = [];
parfor k = 1:nGB
% 	waitbar(k/nGB,f)
	E(k).Ni = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Ni');
	%     E.Al(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Al');
	%     E.Au(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Au');
	%     E.Cu(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Cu');
end
% close(f);
%E.Ni = reshape(E.Ni,size(x)); %x was an output from sphere() in other code

propList = vertcat(E.Ni);
end