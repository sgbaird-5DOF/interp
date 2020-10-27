function proj_down_test(test)
arguments
    test(1,1) double {mustBeInteger} = 2
end

seed=10;
rng(seed)

txtfn = @(lbl) text(0.025,0.95,['(' lbl ')'],'Units','normalized','FontSize',12);
txtfn2 = @(lbl) text(-0.145187165616493,0.974242424021265,0,['(' lbl ')'],'Units','normalized','FontSize',12);

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
			fig.Position = [443.4 123.4 261.6 588.8];
			tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
			nexttile
			t=n2c(pts);
			plot3(t{:},'k*')
			
			[x,y,z]=sphere(80);
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
            txtfn('a')
		end
		
		%remove degenerate dimension (outputting origin point automatically
		%recenters projpts relative to origin)
		[downpts,usv] = proj_down(pts,'zeroQ',true);
        [downpts2,usv2] = proj_down(pts,'zeroQ',false);
		
		if d == 3
			nexttile
			t=n2c(downpts);
			plot(t{:},'k*')
			viscircles([0,0],0.8,'Color','k');
			axis equal tight
            txtfn('b')
%             nexttile
            t=n2c(downpts2);
            hold on
            plot(t{:},'r*')
            plot(0,0,'k+')
            hold off
            legend('zeroQ=T','zeroQ=F','(0,0)','Location','south')
%             viscircles([0,0],0.8,'Color','k');
% 			axis equal tight
%             txtfn('c')
		end
		
		uppts = proj_up(downpts,usv);
		
	case 2
		%sub-hemisphere data (i.e. within single orthant)
		d = 3;
		npts = 100;
		pts = normr(rand(npts,d));
		scl = 0.8;
		if d == 3
			fig=figure;
			fig.Position = [405.8 109.8 208.8 660];
			tiledlayout(4,1,'TileSpacing','compact','Padding','compact')
			nexttile
			t=n2c(pts);
            plot3(t{:},'k.')
% 			K=sphconvhulln(pts);
% 			trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
			
			[x,y,z]=sphere(40);
			x = scl*x;
			y = scl*y;
			z = scl*z;
			hold on
			surf(x,y,z,'EdgeColor','none');
			axis equal tight
			ax=gca;
			ax.View = [65.5094   21.7505];
            txtfn2('a')
		end
		
		projpts = projfacet2hyperplane(mean(pts),pts);
		
		[downpts,usv] = proj_down(projpts,'zeroQ',false);
		
		if d == 3
			nexttile(4)
            t=n2c(downpts);
			K=delaunayn(downpts);
% 			hold on
            triplot(K,t{:},'k')
			axis equal tight
            txtfn('d')
            
            nexttile(3)
            plot(t{:},'k.')
            hold off
            axis equal tight
            txtfn('c')
		end
		
		uppts = proj_up(downpts,usv);
		
		
		if d == 3
			nexttile
			t=n2c(uppts);
            plot3(t{:},'k.')
% 			trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
			
			[x,y,z]=sphere(40);
% 			scl = 0.8;
			x = scl*x;
			y = scl*y;
			z = scl*z;
			hold on
			surf(x,y,z,'EdgeColor','none');
			axis equal tight
			ax=gca;
			ax.View = [65.5094   21.7505];
            txtfn2('b')
		end
		
		%renormalize to hypersphere
		uppts = normr(uppts);
		
		
end

if ~all(ismembertol(pts,uppts,'ByRows',true))
	error('Oops.. proj_down -> proj_up resulted in different points')
else
	disp('Nice! proj_down -> proj_up resulted in same points within default tolerance')
end


