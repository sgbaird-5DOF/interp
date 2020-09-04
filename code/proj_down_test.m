%proj_down test
clear; close all

seed=10;
rng(seed)
test = 1;
switch test
	case 1
		%data with degenerate dimension
		d = 3;
		npts = 30;
		endpts = eye(d);
		endpts(1,:) = [];
		
		v = endpts(2:d-1,:) - endpts(1,:);
		
		%get random coefficients
		randvals = rand(npts,d-2);
		randvals = rand(npts,1).*(randvals./sum(randvals,2));
		randpts = normr(randvals*v + endpts(1,:));
		pts = [endpts; randpts];
		
		if d == 3
            R = rotationmat3D(deg2rad(20),rand(3,1))*rotationmat3D(deg2rad(20),rand(3,1));
            pts = (R*(pts.')).';
			fig=figure;
			fig.Position = [376.5000  364.0000  560.0000  249.0000];
			tiledlayout(1,2)
			nexttile
			t=n2c(pts);
			plot3(t{:},'k*')
			
			[x,y,z]=sphere(40);
			scl = 0.8;
			x = scl*x;
			y = scl*y;
			z = scl*z;
			hold on
			surf(x,y,z,'EdgeColor','none');
			axis equal
			ax=gca;
			ax.View = [-7.0769e+01   2.8197e+01];
			shading interp
		end
		
		%remove degenerate dimension (outputting origin point automatically
		%recenters projpts relative to origin)
		[downpts,usv] = proj_down(pts,'zeroQ',true);
		
		if d == 3
			nexttile
			t=n2c(downpts);
			plot(t{:},'k*')
			viscircles([0,0],0.8,'Color','k');
			axis equal tight
		end
		
		uppts = proj_up(downpts,usv);
		
	case 2
		%sub-hemisphere data (i.e. within single orthant)
		d = 3;
		npts = 100;
		pts = normr(rand(npts,d));
		
		if d == 3
			fig=figure;
			fig.Position = [147.5000  392.5000  972.5000  292.5000];
			tiledlayout(1,3)
			nexttile
			t=n2c(pts);
			K=sphconvhulln(pts);
			trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
			
			[x,y,z]=sphere(40);
			scl = 0.8;
			x = scl*x;
			y = scl*y;
			z = scl*z;
			hold on
			surf(x,y,z,'EdgeColor','none');
			axis equal
			ax=gca;
			ax.View = [65.5094   21.7505];
		end
		
		projpts = projfacet2hyperplane(mean(pts),pts);
		
		[downpts,usv,zeropt] = proj_down(projpts);
		
		if d == 3
			nexttile(3)
			t=n2c(downpts);
			K=delaunayn(downpts);
			triplot(K,t{:})
			axis equal tight
		end
		
		uppts = proj_up(downpts,usv);
		
		
		if d == 3
			nexttile
			t=n2c(uppts);
			trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
			
			[x,y,z]=sphere(40);
			scl = 0.8;
			x = scl*x;
			y = scl*y;
			z = scl*z;
			hold on
			surf(x,y,z,'EdgeColor','none');
			axis equal
			ax=gca;
			ax.View = [65.5094   21.7505];
		end
		
		%renormalize to hypersphere
		uppts = normr(uppts);
		
		
end

if ~all(ismembertol(pts,uppts,'ByRows',true))
	error('Oops.. proj_down -> proj_up resulted in different points')
else
	disp('Nice! proj_down -> proj_up resulted in same points within default tolerance')
end


