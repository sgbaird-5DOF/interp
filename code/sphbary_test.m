% spherical barycentric coordinates test
% clear; close all;

seed = 11;
rng(seed);

P = rand(3);
a = rand(1,3);

c = sphbary(a,P);

%% setup
d = 3; %dimensionality
vertices = rand(d); %rows of vertices defining facet on hypersphere

vertices = normr(vertices); %normalize
nvec = mean(vertices); % point that will define tangent hyperplane
nvec = normr(nvec);

%% projection
newvertices = projfacet2hyperplane(nvec,vertices);

sphbarydata = sphbary(nvec,newvertices) %#ok<*NOPTS>

propvals = 10*rand(1,3);
sphbaryinterp = dot(sphbarydata,propvals)

K = 1:d;

[dataProj,facetPts,dataBary,facetID] = projray2hypersphere(vertices,K,nvec);

dataBary

classicbaryinterp = dot(dataBary,propvals)

if d == 3
	fig = paperfigure();
	hold on
	
	temp = num2cell(newvertices,1);
% 	ax = trisurf(K,temp{:},'FaceColor',[0.4 0.6 0.7]);
% 	alpha(ax,0.7);
	
	temp = num2cell(newvertices,1);
	sz = [length(temp{1}),length(temp)];
	zero = num2cell(zeros(sz),1);
% 	quiver3(zero{:},temp{:},1.1)
% 	plot3(temp{:},'k*')
% 	text(temp{1}+0.02,temp{2}+0.02,temp{3}+0.05,{'1','2','3'})
	
	[x,y,z] = sphere(30);
    ids = find(z < 0);
    x(ids) = 0;
    y(ids) = 0;
    z(ids) = 0;
	ax = surf(x,y,z,'EdgeColor',0.01*[67.8, 84.7, 90.2]);
	alpha(ax,0.2);
    ax = gca;
	ax.View = [45, 15];
	ax.ZLim(1) = 0;
		
	temp = num2cell(vertices,1);
	sz = [length(temp{1}),length(temp)];
	zero = num2cell(zeros(sz),1);
	ax=trisurf(K,temp{:},'FaceColor','black');
    alpha(ax,0.5);
    quiver3(zero{:},temp{:},1.1)
    text(temp{1}+0.02,temp{2}+0.02,temp{3}+0.05,{'1','2','3'})
	
	temp = num2cell(nvec,1);
%  plot3(temp{:},'r*','MarkerFaceColor','r')
% 	plot3(temp{:},'r*')
% 	quiver3(0,0,0,temp{:},1.1)
	
	temp = num2cell(dataProj);
	plot3(temp{:},'r*')
    quiver3(0,0,0,temp{:},1.2)
	
	xlabel('x');
	ylabel('y');
	zlabel('z');

	axis equal tight off
%     
%     camlight
%     lighting gouraud
%     material dull
    
%     set(gcf,'InvertHardCopy','off')
end



%----------------------------------CODE GRAVEYARD--------------------------
%{
	% 	[x,y,z] = sphere(30);
	% 	ax = surf(x,y,z,'EdgeColor','white');
	% 	alpha(ax,0.5)
	
	% 	axis equal tight
	% 	ax = gca;
	% 	ax.View = [60, 45];

% 	t = [0 1.3];
% 	xlim(t);
% 	ylim(t);

% 	quiver3(0,0,0,temp{:},1)

% 	quiver3(zero{:},temp{:},1)


% 	temp = num2cell([nvec; vertices],1);
% 	plot3(temp{:},'b*')

%}


