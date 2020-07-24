%intersect facet test
%%
clear; close all
%%
d = 3;
pts = [eye(d);-eye(d)];
datalist = [1,1,1];
K = convhulln(pts);

intfacetIDs = intersect_facet(pts,K,datalist,1e-6)
intfacetIDs = vertcat(intfacetIDs{:});

%%
clear
%load filename
fname = '5DOF_vtx_octsubdiv1.mat';
%load mesh
load(fname,'pts');
datalist = pts;

%data
d = size(datalist,2);
pts = [eye(d);-eye(d)]; % define orthoplex points
pts = setdiff(pts,datalist,'rows'); % remove orthoplex points from datalist that are already in pts
K = convhulln(pts);

tol = 1e-6;
intfacetIDs = intersect_facet(pts,K,datalist,tol)
intfacetIDs = vertcat(intfacetIDs{:});

%%
clear
%load filename
fname = 'vtxmesh_octsubdiv_nint2.mat';
%load mesh
load(fname,'meshpts');
datalist = meshpts;

%data
d = size(datalist,2);
pts = [eye(d);-eye(d)]; % define orthoplex points
%pts = setdiff(pts,datalist,'rows'); % remove orthoplex points from datalist that are already in pts
K = convhulln(pts);

intfacetIDs = intersect_facet(pts,K,datalist)
intfacetIDs = vertcat(intfacetIDs{:});

%%
clear
%load filename
fname = 'vtxmesh_octsubdiv_nint4.mat';
%load mesh
load(fname,'uniquePts');
datalist = uniquePts;

%data
d = size(datalist,2);
pts = [eye(d);-eye(d)]; % define orthoplex points
pts = setdiff(pts,datalist,'rows'); % remove orthoplex points from datalist that are already in pts
K = convhulln(pts);

intfacetIDs = intersect_facet(pts,K,datalist)
intfacetIDs = vertcat(intfacetIDs{:});

%% plotting
if d == 3
	tmp = num2cell(pts,1);
	trisurf(K,tmp{:},'FaceColor',[0.5 0.5 0.5])
	hold on
	trisurf(K(intfacetIDs,:),tmp{:},'FaceColor','red')
	zero = num2cell(zeros(size(datalist)),1);
	tmp = num2cell(datalist,1);
	quiver3(zero{:},tmp{:},0.5,'black','LineWidth',1,'MaxHeadSize',1)
	ax = gca;
	ax.View = [30,30];
	axis equal off
	hold off
end