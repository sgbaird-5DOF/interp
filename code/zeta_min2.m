function zm = zeta_min2(o1,o2)
arguments
	o1(:,8) double {mustBeFinite,mustBeReal}
	o2(:,8) double {mustBeFinite,mustBeReal}
end
%--------------------------------------------------------------------------
% Date: 2020-07-27
%
% Description: Alternative version of CMU group function zeta_min(),
% vectorized by Sterling Baird
% 
% Inputs:
%		(o1,o2)	-	lists of octonions
%
% Outputs:
%		zm			-	list of minimized zeta angles
%
% Usage:
%		zm = zeta_min2(o1,o2);
%
% Dependencies:
%		*
%
% Notes:
%		*
%--------------------------------------------------------------------------

%quaternion dot products = cos(omega/2), omega is misorientation angle

%unpack quaternions
qA = o1(:,1:4);
qB = o1(:,5:8);
qC = o2(:,1:4);
qD = o2(:,5:8);

%dot() applies conj, so be careful that inputs are real. Alternative: sum(qA.*qC);
qdot_AC = dot(qA,qC,2);
qdot_BD = dot(qB,qD,2);

mu_num1 = qA(:,4).*qC(:,1)-qC(:,4).*qA(:,1)+qB(:,4).*qD(:,1)-qD(:,4).*qB(:,1);
crossAC = cross(qA(:,2:4),qC(:,2:4),2);
crossBD = cross(qB(:,2:4),qD(:,2:4),2);

mu_arg = (mu_num1 + crossAC(:,3) + crossBD(:,3))./(qdot_AC+qdot_BD);
mu = 2*atan(mu_arg);

% shift negative values
zm = mu;
negIDs = mu < 0;
zm(negIDs) = zm(negIDs)+2*pi;

end


%----------------------original code from CMU group------------------------
%{
function zm = zeta_min(qA,qB,qC,qD)
%%% zeta is twist angle of U(1) symmetry (6 --> 5 DOF)
%%% GBOM angle can be analytically minimized w.r.t. zeta (EQN 25, octonion paper)

% [cA,sA,aA,~] = q2ax(qA);
% [cB,sB,aB,~] = q2ax(qB);
% [cC,sC,aC,~] = q2ax(qC);
% [cD,sD,aD,~] = q2ax(qD);

%quaternion dot products = cos(omega/2), omega is misorientation angle

qdot_AC = sum(qA.*qC); % dot(qA,qC);%qdot(cA,cC,sA,sC,aA,aC);
qdot_BD = sum(qB.*qD); %dot(qB,qD);%qdot(cB,cD,sB,sD,aB,aD);

mu_num1 = qA(4)*qC(1)-qC(4)*qA(1)+qB(4)*qD(1)-qD(4)*qB(1);
crossAC = crossp(qA(2:4),qC(2:4));
crossBD = crossp(qB(2:4),qD(2:4));

mu_arg = (mu_num1 + crossAC(3) + crossBD(3))/(qdot_AC+qdot_BD);
mu = 2*atan(mu_arg);

if mu >= 0
	zm = mu;
else
	zm = mu + 2*pi;
end

end
%}

%-----------------------------CODE GRAVEYARD-------------------------------
%{

%}