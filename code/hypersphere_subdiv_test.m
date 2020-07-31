clear; close all force

test = 4;
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
		
	case 4
		%set random number generator
		seed = 10;
		rng(seed);
		
		d = 4;
		%get points along an arc in 3D, plus origin
		endpts = normr(rand(d-1,d));
		nmeshpts = 100;
		
		v = endpts(2:d-1,:) - endpts(1,:);
		
		%get random perturbations
		randvals = rand(nmeshpts,d-2);
		randpts = normr(randvals*v + endpts(1,:));
		meshpts = [endpts; randpts];
		
		% 		meshpts = [endpts; normr(endpts(1,:) + rand(nmeshpts,1)*(endpts(2,:) - endpts(1,:)))];
		
		%get datapoints along same 3D arc
		ndatapts = 100;
		
		%get random perturbations
		randvals = rand(ndatapts,d-2);
		datapts = normr(randvals*v + endpts(1,:));
		% 		datapts = normr(endpts(1,:) + rand(ndatapts,1)*(endpts(2,:) - endpts(1,:)));
		
		if d == 3
			%figure setup
			fig = figure;
			fig.Position = [231.0000  155.5000  789.5000  611.5000];
			tiledlayout(2,2)
			
			%plot 3D points along arc and origin
			nexttile
			t1=n2c(meshpts);
			t2=n2c(datapts);
			plot3(t1{:},'k*',t2{:},'r*')
			axis equal
			hold on
			plot3(0,0,0,'k*')
		end
		
		%remove degenerate dimension (and keep origin point)
		a = [zeros(1,d); meshpts; datapts];
		% 		a = projfacet2hyperplane(normr(mean(meshpts)),a); %turns out this is unnecessary (maybe even a bad idea)
		[a,usv] = proj_down(a,1e-6,struct.empty,1,false,true);
		
		a(2:end,:) = a(2:end,:) - a(1,:);
		a(1,:) = [];
		
		%recenter rest of data relative to origin
		meshpts = a(1:end-ndatapts,:);
		datapts = a(end-ndatapts+1:end,:);
		
		
		% 		[pts,usv] = proj_down(pts,1000);
		% 		pts = normr(proj_up(pts,usv));
		% 		avg = normr(mean(pts));
		
		K = sphconvhulln(meshpts);
		
		if d == 3
			%plot the 2D points
			nexttile
			t1=n2c(round(meshpts,15));
			t2=n2c(round(datapts,15));
			plot(t1{:},'k*',t2{:},'r*')
			
			%plot the 2D triangulation
			hold on
			x = meshpts(:,1);
			y = meshpts(:,2);
			for i = 1:size(K,1)
				plot(x(K(i,:)),y(K(i,:)),'k-')
			end
			ax = gca;
			% 		ax.XLim(min(abs(ax.XLim))) = 0;
			% 		ax.YLim(min(abs(ax.YLim))) = 0;
			plot(0,0,'k*')
			axis equal
		end
		
		%subdivide the 2D arc
		nint = 2;
		[Ktr,K,meshpts] = hypersphere_subdiv(meshpts,K,nint);
		
		if d == 3
			%plot subdivision
			%--2D arc points
			nexttile(4)
			t1=n2c(meshpts);
			t2=n2c(datapts);
			plot(t1{:},'k*',t2{:},'r*')
			mytitle = ['nint == ' int2str(nint)];
			title(mytitle)
			
			%--2D arc triangulation
			hold on
			x = meshpts(:,1);
			y = meshpts(:,2);
			for i = 1:size(K,1)
				plot(x(K(i,:)),y(K(i,:)),'k-')
			end
			plot(0,0,'k*')
			axis equal
		end
		
		intfacetIDs = intersect_facet(meshpts,K,datapts,1e-1,true,'planar','mldivide');
		
		if d == 3
			%uncenter the 2D arc points (using projected origin point) and add
			%projected origin back in
			meshpts = proj_up(meshpts,usv);
			datapts = proj_up(datapts,usv);
			
			%plot the reprojected 3D arc
			nexttile(3)
			t1=n2c(meshpts);
			t2=n2c(datapts);
			plot3(t1{:},'k*',t2{:},'r*')
			title(mytitle)
			hold on
			plot3(0,0,0,'k*')
			axis equal
		end
		nonints = sum(cellfun(@isempty,intfacetIDs));
		disp(['# non-intersections == ' int2str(nonints) '/' int2str(length(intfacetIDs))])
		
end

%----------------------------CODE GRAVEYARD--------------------------------
%{
% 		plot(x(K),y(K)) %only works for closed polygon & corresponding convex hull I suppose
%}
