function proj_down_test(test)
arguments
    test(1,1) double {mustBeInteger} = 2
end

seed=10;
rng(seed)

txtfn = @(lbl) text(0.025,0.95,['(' lbl ')'],'Units','normalized','FontSize',12);
txtfn2 = @(lbl) text(-0.145187165616493,0.974242424021265,0,['(' lbl ')'],'Units','normalized','FontSize',12);

charlist = num2cell('a':'d');

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
            % 			fig=figure;
            % 			fig.Position = [443.4 123.4 261.6 588.8];
            % 			tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
            paperfigure(1,2);
            nexttile % 3D point cloud
            t=n2c(pts);
            plot3(t{:},'k*')
            hold on
            sphplot(200,'axview',[-7.0769e+01   2.8197e+01])
            axis off
            
            papertext(1);
        end
        
        %remove degenerate dimension (outputting origin point automatically
        %recenters projpts relative to origin)
        [downpts,usv] = proj_down(pts,'zeroQ',true);
        [downpts2,usv2] = proj_down(pts,'zeroQ',false);
        
        if d == 3
            nexttile % 2D point cloud
            t=n2c(downpts);
            plot(t{:},'k*')
            viscircles([0,0],0.8,'Color','k');
            scl2 = 1.15;
            axis([-1 1 -1 1]*scl2,'equal','off')
%             txtfn('b')
            papertext(2);
            %             nexttile
            t=n2c(downpts2);
            hold on
            plot(t{:},'r*')
            plot(0,0,'k+')
            hold off
            lgd = legend('zeroQ=T','zeroQ=F','(0,0)','Location','southeast');
            lgd.Position = [0.766477042811879 0.194469508634414 0.140716805329194 0.148823532609379];
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
        
        projpts = projfacet2hyperplane(mean(pts),pts);
        
        [downpts,usv] = proj_down(projpts,'zeroQ',false);
        
        paperfigure(1,3)
        
        scl2 = 1.15;
        
        if d == 3
            tnum = 2;
            nexttile(tnum) %2D triangulation
            t=n2c(downpts);
            K=delaunayn(downpts);
            % 			hold on
            triplot(K,t{:},'k')
            axis(scl2*[-1 1 -1 1],'equal','tight','off')
%             txtfn(charlist{tnum})
            
            papertext(tnum,'xypos',[0.005 1.1])
            
            cloud2D_Q = false;
            if cloud2D_Q
                tnum = 3;
                nexttile(tnum) %2D point cloud
                plot(t{:},'k.')
                hold off
                axis(scl2*[-1 1 -1 1])
                axis equal tight off
                txtfn(charlist{tnum})
            end
        end
        
        uppts = proj_up(downpts,usv);
        
        
        if d == 3
            tnum = 1;
            nexttile(tnum) % 3D point cloud on hyperplane
            t=n2c(uppts);
            plot3(t{:},'k.')
            % 			trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
            
            hold on
            sphplot('axview',[65.5094   21.7505])
            axis off
            
%             [x,y,z]=sphere(40);
%             % 			scl = 0.8;
%             x = scl*x;
%             y = scl*y;
%             z = scl*z;
%             hold on
%             surf(x,y,z,'EdgeColor','none');
%             axis('equal','tight','off')
%             ax=gca;
%             ax.View = [65.5094   21.7505];
%             txtfn2(charlist{tnum})
            papertext(tnum,'xypos',[0.0005 0.95])

        end
        
        if d == 3
            % 			fig=figure;
            % 			fig.Position = [405.8 109.8 208.8 660];
            % 			tiledlayout(4,1,'TileSpacing','compact','Padding','compact')
            tnum = 3;
            nexttile(tnum) % 3D point cloud on hypersphere w/ tri
            t=n2c(pts);
            plot3(t{:},'k.')
            hold on
            
            sphplot()
            
%             K=sphconvhulln(pts);
            trisurf(K,t{:},'FaceColor','none','EdgeColor','k')
            
%             [x,y,z]=sphere(40);
%             x = scl*x;
%             y = scl*y;
%             z = scl*z;
%             surf(x,y,z,'EdgeColor','none');
            axis(scl2*[-1 1 -1 1 -1 1],'equal','tight','off')
            ax=gca;
            ax.View = [65.5094   21.7505];
%             txtfn2(charlist{tnum})
            papertext(tnum);
        end
        
        %renormalize to hypersphere
        uppts = normr(uppts);
        
        
end

if ~all(ismembertol(pts,uppts,'ByRows',true))
    error('Oops.. proj_down -> proj_up resulted in different points')
else
    disp('Nice! proj_down -> proj_up resulted in same points within default tolerance')
end


