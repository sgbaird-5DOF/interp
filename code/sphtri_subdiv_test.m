clear; close all;

%setup
alen = pi/2; %arclength

sphpts = [...
	0		  	pi/2	%pt1
	0			0  	%pt2
	alen		0		%pt3
	]; %[az; el]

% sphpts = [...
%    1.570796326794897   1.570796326794897
%    1.570796326794897                   0
%    3.141592653589793                   0
% 	];
% 
% sphpts = [...
%    3.141592653589793   1.570796326794897
%    3.141592653589793                   0
%    4.712388980384690                   0
% 	];
% 
% sphpts = [...
%    4.712388980384690   1.570796326794897
%    4.712388980384690                   0
%    6.283185307179586                   0
% 	];

nint = 1; % nint == 1 subdivides 0 zero times, nint == 2 subdivides once, etc.

%subdivision
[newsphpts,newpts,TRI] = sphtri_subdiv(sphpts,nint);

%plotting
ax = gca;
trisurf(TRI,newpts(:,1),newpts(:,2),newpts(:,3));
ax.View = [90+rad2deg(alen/2) 45];
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
% t = [-1 1];
% xlim(t);
% ylim(t);
% zlim(t);

%points that should have same values in mesh (might be better just to use
%unique() after applying symmetry operations when doing mirroring)
axismirrorptIDs1 = find(round(newsphpts(:,1)-0,12) == 0); %points on one symmetry axis
axismirrorptIDs2 = find(round(newsphpts(:,1)-alen,12) == 0); %points on other symmetry axis

axismirrorptIDs = [axismirrorptIDs1;axismirrorptIDs2]; %combine

mirrorptIDs = setdiff(1:length(newsphpts),axismirrorptIDs); %find the points not on a symmetry axis