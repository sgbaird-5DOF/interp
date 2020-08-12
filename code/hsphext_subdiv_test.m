% hypersphere exterior hull subdivision test
clear; close all

test = 2;

addpathdir({'q2rod.m','cu2qu.m','GBfive2oct.m'})

switch test
	case 1
		fname = '5DOF_oct_vtx_octsubdiv1.mat';
		
		%load data
		load(fname,'pts');
		
		%parameters
		octsubdiv = 3; %single subdivision (distinct from levels in K-tree)
		
		[Ktr, K3, meshpts] = hsphext_subdiv(pts,octsubdiv);
		
	case 2
		seed = 10;
		rng(seed);
		
		d = 7;
		
% 		pts = [eye(3); 1 1 1; 0 -1 0];
% 		pts = [eye(3); 1 1 1; 0 -1 0; -0.5 0 0.5];
% 		pts = [eye(3); 0 -1 0; -1 0 0]; %produces a 2D ring
% 		pts = rand(10,d);
% 		pts = [eye(3); rand(10000,d)];
		npts = 2000;
		
		if d == 7
			pts = get_ocubo(npts,'random',[],seed);
			pts = get_octpairs(pts);
			[pts,usv] = proj_down(pts);
		else
			pts = rand(npts,d);
		end
		
		pts = normr(pts);
		
		K1 = sphconvhulln(pts);
		
		nint = 1;
		[Ktr,K,meshpts] = hsphext_subdiv(pts,nint,true);
		
		disp(['d: ' int2str(d) ', points kept: ' int2str(size(meshpts,1)) '/' int2str(npts)])
		
		if d == 3
			close(figure(1));
			fig = figure(1);
			fig.Position = [247.5000  284.0000  822.0000  360.5000];
			tiledlayout(1,3)
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
			
			hold on
			for i = 1:size(K,1)
				t1 = meshpts(K(i,1),:);
				t2 = meshpts(K(i,2),:);
				tmp = n2c([t1;t2]);
				% plot3(tmp{:},'-')
				plot3(tmp{:},'-')
				plot3(tmp{:},'k*')
				drawnow
				pause(0.1)
			end
% 			tmp = num2cell(meshpts,1);
%  			plot3(tmp{:},'*')
			hold on
% 			trisurf(K,tmp{:})
			ax.View = [105 20];
			axis equal tight
			
			nexttile
			avg = normr(mean(meshpts));
			a = projfacet2hyperplane(avg,meshpts);
			b = proj_down([avg;a],1e-4);
% 			b = proj_down(a);
			t=n2c(b);
			triplot(delaunayn(b),t{:})
			axis equal tight
			
		end
		
		if d == 7
			exto = proj_up(meshpts,usv); %exterior octonions
			five = GBoct2five(exto,true);
			five = correctdis(five);
			geometry = findgeometry(vertcat(five.q),2);
			ninterior = sum(ismember(geometry,'interior')) %#ok<NOPTS>
			[five.geometry] = geometry{:};
			fig=figure;
			fig.Position=[162.0000  385.0000  893.0000  384.5000];
			tiledlayout(1,2)
			plot5DOF(five)
		end
		
	case 3
		%probing a convex hull
		seed = 10;
		rng(seed);
		
		d = 3;
		
		%mesh
		nmeshpts = 100;
		orthoPts = orthoplex(d);
		ids = find(all(orthoPts >= 0, 2));
		orthoPts = orthoPts(ids,:);
		meshpts = normr(rand(nmeshpts,d));
		[~,K,pts] = hsphext_subdiv(meshpts);
		% meshprops = (1:nmeshpts)/nmeshpts*10 + rand(1,nmeshpts);
		
		randpoly = rand(d,1);
		meshprops = meshpts*randpoly;
		
		%data
		ndatapts = 10;
		datapts = normr(rand(ndatapts,d));
		dataprops = datapts*randpoly;
end

%---------------------------CODE GRAVEYARD---------------------------------
%{
		fnamelist = {'q2rod.m'};
		for i = 1:length(fnamelist)
			fname = fnamelist{i};
			file = dir(fullfile('**',fname));
			addpath(file.folder);
		end
%}