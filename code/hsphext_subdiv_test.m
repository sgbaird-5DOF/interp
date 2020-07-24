% hypersphere exterior hull subdivision test
clear; close all

test = 2;
switch test
	case 1
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
		octsubdiv = 3; %single subdivision (distinct from levels in K-tree)
		
		[Ktr, K3, meshpts] = hsphext_subdiv(pts,octsubdiv);

		
	case 2
		seed = 10;
		rng(seed);
		
		d = 3;
		
% 		pts = [eye(3); 1 1 1; 0 -1 0];
		
% 		pts = [eye(3); 1 1 1; 0 -1 0; -0.5 0 0.5];
		
% 		pts = [eye(3); 0 -1 0; -1 0 0]; %produces a 2D ring
		
% 		pts = rand(10,d);

% 		pts = [eye(3); rand(100,d)];

		pts = rand(1000,d);
		
		
		pts = normr(pts);
		
		K1 = sphconvhulln(pts);
		
		nint = 3;
		
		[Ktr,K,meshpts] = hsphext_subdiv(pts,nint,true);
		
		
		if d == 3
			close(figure(1));
			fig = figure(1);
			fig.Position = [247.5000  284.0000  822.0000  360.5000];
			tiledlayout(1,2)
			ax = nexttile;
			tmp = num2cell(pts,1);
			trisurf(K1,tmp{:});
			ax.View = [105 20];
			axis equal tight
			
			% 	ax = nexttile;
			% 	tmp = num2cell(projpts,1);
			% 	plot(tmp{:},'*')
			% 	ax.View = [135 20];
			% 	axis equal tight
			
			ax = nexttile;
			tmp = num2cell(meshpts,1);
 			plot3(tmp{:},'*')
			hold on
			trisurf(K,tmp{:})
			ax.View = [105 20];
			axis equal tight
			
		end
end