%test meshBP

% dependency: disorientation.m for 'interior' option
clear; close all

%assumption is 'FCC' or 'BCC' crystal type (2020-06-30)

% geomchecktypeList = {...
% 	'OAB','OBCE','OADE','CDE',... surfaces
% 	'OB','CE','ED','OE','OA','AC',... lines
% 	'B','E','A','C','O','interior'}; %points
 
% geomchecktypeList = {...
% 	'A','B','C','D','E','O'};

% geomchecktypeList = {...
% 	'C','B','interior'}; %points

geomchecktypeList = {'twosphere'};

% geomchecktypeList = {...
% 	'interior'}; %points

%add info for quaternions and rodriguez vectors for key features of the
%misorientation FZ
filename = 'misFZfeatures.mat'; %comes from plotFZrodriguez_test.m
file = dir(filename);
addpath(file(1).folder);
%qlist === quaternion list, dlist === rodriguez vec list
load(filename,'qlist','dlist','qnames','lnames','snames');

nint = 3; % nint == 1 subdivides 0 zero times, nint == 2 subdivides once, etc.

k = sqrt(2)-1;

figure(1)
plotFZrodriguez();

scl = 0.1;
rx = [1,0,0]*scl;
ry = [0,1,0]*scl;
rz = [0,0,1]*scl;

w = [0,0,0];
hold on
quiver3(w,w,w,rx,ry,rz,1,'linewidth',1)
text(rx(1),rx(2)+0.01,rx(3),'d_1','FontWeight','bold')
text(ry(1),ry(2),ry(3)+0.01,'d_2','FontWeight','bold')
text(rz(1)+0.01,rz(2),rz(3),'d_3','FontWeight','bold')

%figure(2) options
fig2 = figure(2);
fig2.Position = [162,66,980.5,758.5];

nchecks = length(geomchecktypeList);
sq = sqrt(nchecks);

numsubplotQ = true;

if numsubplotQ
	p = numSubplots(nchecks);
	rows = p(1);
	cols = p(2);
else
	if mod(sq,1) == 0
		rows = sq;
		cols = sq;
	else
		rows = ceil(sq);
		cols = ceil(sq);
	end
end

ctrcuspQ = false;

for i = 1:length(geomchecktypeList)
	geomchecktype = geomchecktypeList{i};
	q = qlist.(geomchecktype);
	%various test cases
	q = q/norm(q); %normalize quaternion
	%q = disorientation(q,'cubic');
	
	figure(1)
	hold on
	d = q(2:4)./q(1);
	plot3(d(1),d(2),d(3),'ko','markerfacecolor','none')
	text(d(1),d(2),d(3)+0.01,geomchecktype)
	
	[newpts,A,R,TRI,alen] = meshBP(q,nint,ctrcuspQ);
	
	figure(2)
	hax(i) = subplot(rows,cols,i); %axis handle
	%hax.View = [rad2deg(alen)+45,45];
	hax(i).View = [135,45];
	
	ax = A(:,1);
	ay = A(:,2);
	az = A(:,3);
	
	hold on
	%sample frame (symmetry axis)
	quiver3(0,0,0,1,0,0,1,'k','linewidth',2)
	quiver3(0,0,0,0,1,0,1,'k','linewidth',2)
	quiver3(0,0,0,0,0,1,1,'k','linewidth',2)
	text(1.01,0,0,'x_s')
	text(0,1.01,0,'y_s')
	text(0,0,1.01,'z_s')
	
	%lab frame (symmetry axis)
	quiver3(0,0,0,ax(1),ax(2),ax(3),1,'r','linewidth',2)
	quiver3(0,0,0,ay(1),ay(2),ay(3),1,'r','linewidth',2)
	quiver3(0,0,0,az(1),az(2),az(3),1,'r','linewidth',2)
	text(ax(1),ax(2),ax(3),'x_l')
	text(ay(1),ay(2),ay(3),'y_l')
	text(az(1),az(2),az(3),'z_l')
	
	%plot sphere
	[x,y,z] = sphere(40);
	scl = 0.5;
	x = scl*x;
	y = scl*y;
	z = scl*z;
	C(:,:,1) = zeros(size(x))+0.5;
	C(:,:,2) = C(:,:,1)+0.5;
	C(:,:,3) = 0.5*ones(size(x)); % blue
	surf(x,y,z,C)
	shading interp;
	
	if isempty(newpts(1).sub)
		field = 'main';
		n = 1;
	else
		field = 'sub';
		n = length(newpts);
	end
	for j = 1:n
		newptVals = newpts(j).(field);
		TRIvals = TRI(j).(field);
		%point(s) in sample frame
		pts = (R\newptVals.').';
		
		x = pts(:,1);
		y = pts(:,2);
		z = pts(:,3);
		
		w = repelem(0,length(x)).';
		
		quiver3(w,w,w,x,y,z,1,'k','linewidth',1)
		plot3(x,y,z,'k*')
		trisurf(TRIvals,x,y,z,'FaceColor','none','EdgeColor','black')
		
		%point(s) in lab frame
		x = newptVals(:,1);
		y = newptVals(:,2);
		z = newptVals(:,3);
		
% 		quiver3(w,w,w,x,y,z,1,'r','linewidth',1)
		plot3(x,y,z,'r*')
		trisurf(TRIvals,x,y,z,'FaceColor','none','EdgeColor','red')
		
		axis equal vis3d
		xlabel('x')
		ylabel('y')
		zlabel('z')
		t = [-1 1];
		xlim(t)
		ylim(t)
		zlim(t)
		axis off tight
	end
	
	%remove whitespace
	outerpos = hax(i).OuterPosition;
	ti = hax(i).TightInset;
	left = outerpos(1) + ti(1);
	bottom = outerpos(2) + ti(2);
	hax_width = outerpos(3) - ti(1) - ti(3);
	hax_height = outerpos(4) - ti(2) - ti(4);
	hax(i).Position = [left bottom hax_width hax_height];
	
	title(geomchecktype);
	
end

%link the camera view angles, etc.
% linknames = {'CameraUpVector', 'CameraPosition','CameraTarget','CameraViewAngle'};
% Link = linkprop(hax, linknames);
% setappdata(gcf, 'StoreTheLink', Link);

%----------------------------CODE GRAVEYARD--------------------------------
%{
switch geomchecktype %surfaces case 'OAB' q = [1,0.6,1,0]; case 'OBCE' q =
[3,0.5,0.5,1]; case 'OADE' q = [3 1 2 2]; case 'CDE' q = [3,0.5,1.5,1];
%lines case 'OB' q = [0.5,1,1,0]; case 'CE' q = [5,1,1,3]; case 'ED' q = [5
1 2 2]; case 'OE' q = [0.5,1,1,1]; case 'OA' q = [0.5,1,0,0]; case 'AC' q =
[1,sqrt(2)-1,2,(sqrt(2)-1)*2]; %points case 'B' q =
[1/sqrt(1+2*k^2),k/sqrt(1+2*k^2),k/sqrt(1+2*k^2),0]; case 'E' q =
[sqrt(3)/2,1/(2*sqrt(3)),1/(2*sqrt(3)),1/(2*sqrt(3))]; case 'A' q =
[cos(pi/8),sin(pi/8),0,0]; case 'C' q = 1/(2*sqrt(2))*[1/k,1,1,k]; case 'O'
q = [1,0,0,0]; case 'interior' q = [0.7489,0.2646,0.2281,0.1336]; %randomly
generated (and then plugged into disorientation(q,'cubic') %q =
[1,0.5,1,0]+[3,0.5,1.5,1]; %OAB + CDE examples end

Even though I handpicked these values based on Patala2013, many of them
didn't seem to be in the place I would have expected in the misorientation
FZ.

%}

