clear; close all force

test = 2;
switch test
	case 1
		%%
		%fname = 'hypersphere_subdiv_testdata.mat';
		% fname note: facet 88001 of 'hypersphere_subdiv_testdata.mat' had an extra
		% degenerate dimension. This only happened once on the way to 88001 (out of
		% ~90000), but I may consider looking at the values of the SVD
		% (specifically diag(S) ) to get a sense of how "skew" the facets are.
		% pca(meshpts) would be another option.
		
		fname = '5DOF_oct_vtx_octsubdiv1.mat';
		
		%load data
		load(fname,'pts');
		
		fnamelist = {'q2rod.m'};
		for i = 1:length(fnamelist)
			fname = fnamelist{i};
			file = dir(fullfile('**',fname));
			addpath(file.folder);
		end
		
		%parameters
		octsubdiv = 2; %single subdivision (distinct from levels in K-tree)
		
		K = sphconvhulln(pts,true);
		[Ktr, K3, meshpts] = hypersphere_subdiv(pts,K,octsubdiv);
		
		
	case 2
		%%		
		seed = 10;
		rng(seed);
		
		d = 3;
		n = 3;
		subdivType = 'orthant';
		meshstart = hypersphereSetup(n,d,subdivType);
		K = sphconvhulln(meshstart);
		
		nint = 3;
		
		delaunayQ = true;
		[Ktr, K3, meshpts] = hypersphere_subdiv(meshstart,K,nint);
		
		if d == 3
			figure(2)
			ax = subplot(1,2,1);
			tmp = num2cell(meshstart,1);
			trisurf(K,tmp{:});
			ax.View = [135 20];
			axis equal tight
			
			ax = subplot(1,2,2);
			tmp = num2cell(meshpts,1);
			trisurf(K3,tmp{:})
			ax.View = [135 20];
			axis equal tight
			
		end
		
	case 3
		%%
		seed = 10;
		rng(seed);
		
		d = 3;
		
		% pts = rand(10,3);
		pts = [eye(3); 1 1 1];
		
		pts = normr(pts);
		
		K1 = sphconvhulln(pts);
		
		projpts = sphere_stereograph(pts);
		
		[projpts,usv] = proj_down(projpts,1e-6);
		
		% get exterior
		K2 = convhulln(projpts);
		
		nint = 3;
		
		[Ktr, K3, meshpts] = hypersphere_subdiv(projpts,K2,nint);
		
		% meshpts2 = sphere_stereograph_inverse(meshpts);
		
		meshpts2 = proj_up(meshpts,usv);
		
		meshpts3 = sphere_stereograph_inverse(meshpts2);
		
		
		if d == 3
			close(figure(1)); figure(1)
			ax = nexttile;
			tmp = num2cell(pts,1);
			trisurf(K1,tmp{:});
			ax.View = [135 20];
			axis equal tight
			
			ax = nexttile;
			tmp = num2cell(projpts,1);
			plot(tmp{:},'*')
			ax.View = [135 20];
			axis equal tight
			
			ax = nexttile;
			tmp = num2cell(meshpts3,1);
			plot3(tmp{:},'*')
			ax.View = [135 20];
			axis equal tight
			
		end
		
end
