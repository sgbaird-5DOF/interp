% project ray to hypersphere test
% dependencies: projray2hypersphere.m, sphconvhulln.m, normr.m

clear; close all;

d = 3;

seed = 12;
rng(seed);

nmeshpts = 20;

meshpts = rand(d,nmeshpts).'; %transpose to preserve random number generation order
meshpts = normr(meshpts);
K = sphconvhulln(meshpts);

data = mean(meshpts,1);
datanorm = normr(data);

if d == 3
	tmp = num2cell(meshpts,1);
	trisurf(K,tmp{:},'FaceColor','white');
	hold on
	[x,y,z] = sphere(40);
	x = 0.5*x;
	y = 0.5*y;
	z = 0.5*z;
	surf(x,y,z,'EdgeColor','none');
	tmp = num2cell(datanorm,1);
	quiver3(0,0,0,tmp{:},1.2,'LineWidth',1,'MaxHeadSize',0.5)
	axis equal tight vis3d
	ax = gca;
	ax.View(1) = 50;
end

[dataProj,facetPts,dataBary,facetIDs] = projray2hypersphere(meshpts,K,datanorm);

if d == 3
	tmp = num2cell(meshpts,1);
	trisurf(K(facetIDs,:),tmp{:},'FaceColor','red')
end


% fname = 'vtxmesh.mat';
% 
% %load data
% load(fname,'K','meshpts');
