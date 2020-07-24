clear; close all;

% fname = 'vtxmesh.mat';
% fname = '5DOF_vtx_deleteO_octsubdiv1.mat';
% fname = '5DOF_vtx_octsubdiv1.mat';
% fname = '5DOF_exterior_octsubdiv1.mat';
% fname = '5DOF_exterior_mis15_nint3_octsubdiv2.mat' %no need to further subdivide
% fname = '5DOF_exterior_res15_nint2_octsubdiv1';
% fname = '5DOF_misFZfeatures_octsubdiv3.mat';
% fname = '5DOF_interior_octsubdiv1.mat';
% fname = '5DOF_vtx_deleteO_octsubdiv5.mat'; %no subdivision needed
% fname = '5DOF_vtx_octsubdiv1.mat';
% fname = '5DOF_oct_vtx_octsubdiv1.mat';
% fname = '5DOF_hsphext_octsubdiv2.mat';
% fname = '5DOF_hsphext_octsubdiv3.mat';
% fname = '5DOF_exterior_hsphext_octsubdiv1.mat';
% fname = '5DOF_exterior_hsphext_res10_nint1_octsubdiv1.mat';

% fname = get_fname('5DOF_exterior_hsphext',10,1,3); %res,nint,octsubdiv

fname = get_fname('5DOF_exterior_hsphext',13,2,2); %res,nint,octsubdiv


% S = load('misFZfeatures.mat');

%load data
load(fname,'pts','sphK','usv');

addpathdir({'qinv.m','qu2ax.m','q2rod.m','PGnames.mat'});

%parameters
octsubdiv = 1; % # of intervals along facet edge (distinct from levels in K-tree)

if isempty(usv)
	[projpts,usv] = proj_down(pts,1e-6);
else
	[projpts,usv] = proj_down(pts,1e-6,usv);
end

if ~isempty(projpts)
	pts = projpts;
	stereoQ = false;
else
% 	projpts = sphere_stereograph(pts);
	projpts = projfacet2hyperplane(mean(pts),pts);
	stereoQ = true;
end

if (octsubdiv > 1) % && (exist('sphK','var') == 0)
	[Ktr, K, pts] = hypersphere_subdiv(pts,sphK,octsubdiv);
elseif (octsubdiv == 1) && (exist('sphK','var') == 0)
	
	if ~contains(fname,'hsphext')
		K = sphconvhulln(pts);
	else
		[Ktr, K, pts] = hsphext_subdiv(pts,sphK,octsubdiv);
	end
	
else
	K = sphK;
end

if stereoQ
% 	pts = sphere_stereograph_inverse(projpts);
	pts = normr(pts);
	pts = proj_up(pts,usv);
end

% if size(pts,2) == 7
% 	pts = [pts zeros(size(pts,1),1)]; %add row of zeros
% end

if size(pts,2) == 7
	pts = proj_up(pts,usv);
end

tic
disp('octlist to five')

% FUNCTION CALL
five = GBoct2five(pts);

npts = length(five);
disp(['npts = ',int2str(npts)])

toc




%% plot setup
fig = figure;
fig.Position = [340.0000  393.5000  560.0000  252.5000];
tiledlayout(1,2)
%---misorientation FZ---
nexttile
%plotting
plotFZrodriguez_vtx();
hold on

%plotting individual facets
% for i = nrange
% 	dnewtemp = dnew{i};
% 	npts2 = size(dnewtemp,1);
% 	tmp = num2cell(dnewtemp,1);
% 	scatter3(tmp{:},10,repelem(i,npts2,1),'filled')
% 	drawnow
% 	pause(0.5)
% end

% q_all = disorientation(vertcat(five.q),'cubic');
% d_all = q2rod(q_all);

% d_all = vertcat(five.d);

% tmp = num2cell(d_all,1);

% ids = strcmp('C',cellstr(vertcat(five.geometry)));
% pts = pts(ids,:);
% 
% 
% pgnum = 32;
% olist = osymsets(pts,pgnum);
% olist2 = vertcat(olist{:});
% five = GBoct2five(olist2);

q_all = vertcat(five.q);
q_all = disorientation(q_all,'cubic');
d_all = q2rod(q_all);
nA_all = vertcat(five.nA);

%set rng for rangesearch
seed = 10;
rng(seed);

%load test point
load('misFZfeatures.mat','dlist')
r = 0.1; %Euclidean distance
% r = Inf;
Y = dlist.A;

%check for points within specified distance
[idx,D] = rangesearch(d_all,Y,r);

d_all = d_all(idx{1},:);
nA_all = nA_all(idx{1},:);

tmp = num2cell(d_all,1);

C = parula(size(nA_all,1));
[~,I] = sort(D{1});
C = C(I,:);

% scatter3(tmp{:},10,(1:size(d_all,1)).','filled')
ax = scatter3(tmp{:},10,C,'filled'); %color by Euclidean distance
alpha(ax,0.5);

zero = num2cell(zeros(size(d_all,1),3),1);
quiver3(tmp{1},zero{1:3},tmp{2},zero{1},'Color','k','ShowArrowhead','off',...
	'LineStyle',':','AutoScale','off')
quiver3(tmp{1:2},zero{:},tmp{3},'Color','k','ShowArrowhead','off',...
	'LineStyle',':','AutoScale','off')

% Kd = convhulln(d_all);
% tmp = num2cell(d_all,1);
% trisurf(Kd,tmp{:},'FaceColor','none')

%---Boundary Plane Normals---
nexttile
hold on

[x,y,z] = sphere(40);
scl = 0.65;
x = x*scl;
y = y*scl;
z = z*scl;

surf(x,y,z,'EdgeColor','none','FaceColor','cyan','FaceAlpha',0.7)
ax = gca;
ax.View = [135 15];

% for i = nrange
% 	nAtemp = vertcat(five2{i}.nA);
% 	npts2 = length(nAtemp);
% 	tmp = num2cell(nAtemp,1);
% 	scatter3(tmp{:},10,repelem(i,npts2,1),'filled')
% 	drawnow
% 	pause(0.5)
% end

% nAfull = vertcat(vertcat(five2{:}).nA);

% nA_all = vertcat(five.nA);

tmp = num2cell(nA_all,1);
% scatter3(tmp{:},10,(1:size(nA_all,1)).','filled')
ax = scatter3(tmp{:},10,C,'filled'); %color by Euclidean distance, larger dist == smaller marker
alpha(ax,0.3);

KnA = sphconvhulln(nA_all);
tmp = num2cell(nA_all,1);
trisurf(KnA,tmp{:},'FaceColor','none')

scl = 0.75;
x = [1,0,0]*scl;
y = [0,1,0]*scl;
z = [0,0,1]*scl;

w = [0,0,0];
hold on
quiver3(w,w,w,x,y,z,1,'linewidth',1)
text(x(1),x(2)+0.01,x(3),'d_1','FontWeight','bold')
text(y(1),y(2),y(3)+0.01,'d_2','FontWeight','bold')
text(z(1)+0.01,z(2),z(3),'d_3','FontWeight','bold')

axis equal tight vis3d off


sgtitle(fname,'interpreter','none')


%--------------------------------CODE GRAVEYARD----------------------------
%{
nptstot = 0;
for i = 1:nmain
	npts2 = size(Ktr.sub{i}.main.pts);
	
	nrange = nptstot+1:nptstot+npts2;
	five2{i} = five(nrange);
	nptstot = nptstot+npts2;
end



load('misFZfeatures.mat','qlist')
qint = qlist.interior;
nAint = normr([0 0 1]);
oint = GBfive2oct(qint,nAint);

% octnew = zeros(size(octlist,1),16);
% for i = 1:length(octlist)
% 	[~, octnew(i,:), ~] = GBdist([octlist(i,:) oint],32,false);
% end
% 
% octnew = octnew(i,1:8);


%e.g. fname == '5DOF_vtx_octsubdiv1.mat' -> fname == '5DOF_vtx_octsubdiv4'
fname = [fname(1:strfind(fname,'_octsubdiv')) 'octsubdiv' int2str(octsubdiv)];
save(fname,'Ktr','K','meshpts')

% [row,col] = find(isnan(nAfull));



%plotting individual facets
% nmain = length(Ktr.sub);
% nrange = 1:nmain;
% nrange = 5;
% nrange = 1:5;
% for i = nrange
% 	five2{i} = GBoct2five(Ktr.sub{i}.main.pts);
% % 	qnew{i} = disorientation(vertcat(five2{i}.q),'cubic');
% % 	dnew{i} = q2rod(qnew{i});
% 	dnew{i} = vertcat(five2{i}.d);
% end


fnamelist = {'qinv.m','qu2ax.m','q2rod.m','PGnames.mat'};
for i = 1:length(fnamelist)
	fname = fnamelist{i};
	file = dir(fullfile('**',fname));
	addpath(file.folder);
end

%}