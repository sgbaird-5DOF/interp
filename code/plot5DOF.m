function plot5DOF(five,lbl,opts,ilist,NV)
arguments
	five struct
	lbl char {mustBeMember(lbl,{'mesh','data','OSLERP','ocubo'})} = 'mesh'
	opts = struct.empty
	ilist = []
	NV.quiverQ(1,1) logical = false
end

%plotting
nexttile
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

d_all = vertcat(five.d);

Q1 = strcmp(lbl,'mesh') || isempty(ilist);
Q2 = strcmp(lbl,'mesh');

if strcmp(lbl,'ocubo')
	alphaval = 1;
else
	alphaval = 0.5;
end

%misorientation
if Q1
	tmp = num2cell(d_all,1);
	scatter3(tmp{:},5,(1:size(d_all,1)).','filled','markerfacealpha',alphaval)
	
	zero = num2cell(zeros(size(d_all,1),3),1);
	tmp = num2cell(d_all,1);
	if NV.quiverQ
		quiver3(tmp{1},zero{1:3},tmp{2},zero{1},'Color','k','ShowArrowhead','off',...
			'LineStyle',':','AutoScale','off')
		quiver3(tmp{1:2},zero{:},tmp{3},'Color','k','ShowArrowhead','off',...
			'LineStyle',':','AutoScale','off')
	end
	
	if strcmp(lbl,'OSLERP')
		plot3(tmp{:})
	end
else
	% ilist specified
	intIDs = setdiff(1:size(d_all,1),ilist);
	tmp = num2cell(d_all(intIDs,:),1);
	ax1 = scatter3(tmp{:},5,'b','filled','markerfacealpha',0.25);
	tmp = num2cell(d_all(ilist,:),1);
	ax2 = scatter3(tmp{:},5,'r','filled','markerfacealpha',0.25);
	legend([ax1 ax2],{'intersecting','non-intersecting'},'location','southoutside');
end

% legend('intersecting','non-intersecting')
title([lbl ' misFZ'])

nexttile
hold on

%boundary plane
[x,y,z] = sphere(100);
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

camlight('right')

% nAfull = vertcat(vertcat(five2{:}).nA);
nAfull = vertcat(five.nA);

if Q1
	tmp = num2cell(nAfull,1);
	scatter3(tmp{:},10,(1:size(nAfull,1)).','filled','markerfacealpha',alphaval)
	if strcmp(lbl,'OSLERP')
		plot3(tmp{:})
	end
	if Q2
		KnA = sphconvhulln(nAfull);
		disp(' ')
		tmp = num2cell(nAfull,1);
% 		trisurf(KnA,tmp{:},'FaceColor','none')
	end
else
	tmp = num2cell(nAfull(intIDs,:),1);
	ax1 = scatter3(tmp{:},5,'b','filled','markerfacealpha',0.2);
	tmp = num2cell(nAfull(ilist,:),1);
	ax2 = scatter3(tmp{:},5,'r','filled','markerfacealpha',0.2);
	legend([ax1 ax2],{'intersecting','non-intersecting'},'location','southoutside');
end

axis equal tight vis3d off

title([lbl ' BP Normal'])

%--------------------------------CODE GRAVEYARD----------------------------
%{

% if nargin == 3
% 	opts = varargin{1};
% end
% 
% if nargin == 4
% 	opts = varargin{1};
% 	ilist = varargin{2};
% end


if ~isempty(opts)
	title([lbl '_res' int2str(opts.res)],'interpreter','none')
else
	title(lbl)
end


if ~isempty(opts)
	title([lbl '_nint' int2str(opts.nint)],'interpreter','none')
else
	title(lbl)
end
%}
