clear; close all

frametype = 'sample'; %'sample', 'lab'

filepathgen = fullfile('**','gmat2q.m');
file = dir(filepathgen);
if ~strcmp(pwd,file.folder)
	addpath(file.folder);
end

%% Test for symmetry in GB FZ

dO = [0, 0, 0];
dA = [sqrt(2)-1, 0, 0];
dB = [sqrt(2)-1, sqrt(2)-1, 0];
dC = [sqrt(2)-1, sqrt(2)-1, 3-2*sqrt(2)];
dD = [sqrt(2)-1, 1-(1/sqrt(2)), 1-(1/sqrt(2))];
dE = [1/3, 1/3, 1/3];

% d = (dO+dA+dB)/3;
d = dC;
q(1) = sqrt(1/(1+sum(d.^2)));
q(2) = q(1)*d(1);
q(3) = q(1)*d(2);
q(4) = q(1)*d(3);

% Choose points in misorietation FZ
% q = rot2q(30*degree,90*degree,22.5*degree); %(1) OAB
% q = rot2q(30*degree,90*degree,45*degree); %(2) OB
% q = rot2q(30*degree,90*degree,45*degree); %(3) OE
q = disorientation(q,'cubic');

figure
plotFZrodriguez();
hold on
d = q(2:4)./q(1);
plot3(d(1),d(2),d(3),'ko','markerfacecolor','r')

% Define symmetry axes of GB FZ
% % (8) OE
% ax = (1/sqrt(2))*[q(1)+q(2),q(2)-q(1),-2*q(2)].';
% ay = (1/sqrt(6))*[-3*q(2)+q(1),3*q(2)+q(1),-2*q(1)].';
% az = (1/sqrt(3))*[1,1,1].';

% % (9) OA
% ax = [0,q(1),q(2)].';
% ay = [0,-q(2),q(1)].';
% az = [1,0,0].';

% % (12) E
% ax = (1/sqrt(6))*[2,-1,-1].';
% ay = (1/sqrt(2))*[0,1,-1].';
% az = (1/sqrt(3))*[1,1,1].';

% (14) C
ax = [1,0,0].';
ay = (1/sqrt(2))*[0,1,1].';
az = (1/sqrt(2))*[0,-1,1].';

% Generate GB normals


[x,y,z] = sphere(40);
% x = x(21:end,:);
% y = y(21:end,:);
% z = z(21:end,:);

nA_Lab = [x(:),y(:),z(:)].';
switch frametype
	case 'lab'
		nA_Lab = [ax,ay,az]*nA_Lab;
	case 'sample'
		nA_Lab = eye(3)*nA_Lab;
end

% Compute GB matrices

qA_Lab = repmat([1 0 0 0],size(nA_Lab,2),1);
qB_Lab = repmat(q,size(nA_Lab,2),1);

[gA_R,gB_R] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab,'livermore');

%% Calculate GB Energies

f = waitbar(0,['calculating GB energies for ',int2str(length(x)),' points.']);
nGB = size(qA_Lab,1);
for k = 1:nGB
	waitbar(k/nGB,f)
	E.Ni(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Ni');
	%     E.Al(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Al');
	%     E.Au(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Au');
	%     E.Cu(k) = GB5DOF(gA_R(:,:,k),gB_R(:,:,k),'Cu');
end
close(f);
E.Ni = reshape(E.Ni,size(x));

%% Plot

figure;
surf(x,y,z,E.Ni);
shading interp;
axis equal tight vis3d
xlabel('x')
ylabel('y')
zlabel('z')
hold on

switch frametype
	case 'lab'
		%use standard x,y,z cartesian axes
		quiver3(0,0,0,1,0,0,1.5,'r','linewidth',2)
		quiver3(0,0,0,0,1,0,1.5,'g','linewidth',2)
		quiver3(0,0,0,0,0,1,1.5,'b','linewidth',2)
	case 'sample'
		%use symmetry axis
		quiver3(0,0,0,ax(1),ax(2),ax(3),2,'r','linewidth',2)
		quiver3(0,0,0,ay(1),ay(2),ay(3),2,'g','linewidth',2)
		quiver3(0,0,0,az(1),az(2),az(3),2,'b','linewidth',2)
end
