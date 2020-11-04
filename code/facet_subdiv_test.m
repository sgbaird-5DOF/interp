%SUBDIV_TEST  e.g. subdivide a d-1 or d-2 simplex (e.g. triangle or line in 3D Cartesian coordinates) and plot
clear; close all;

test = 2;

switch test
	case 1
		% subdivide a d-1 simplex (e.g. triangle in 3D Cartesian coordinates)
		d = 3;
		
		ptType = 'regular'; %'regular', 'thin'
		switch ptType
			case 'regular'
				pts = -eye(d)+1;
			case 'thin'
				pts = eye(d);
				pts(1,1) = 0.1*pts(1,1);
				pts(2,2) = 0.1*pts(2,2);
		end
		
		nint = 4;
		
		[newpts,TRI] = facet_subdiv(pts,nint,true);
		
		switch d
			case 2
				trisurf(TRI,newpts(:,1),newpts(:,2),zeros(1,length(newpts)))
			case 3
				trisurf(TRI,newpts(:,1),newpts(:,2),newpts(:,3))
			case 4
				%depricated (subdivpts is now internal variable in facet_subdiv 2020-06-27)
				%K = convhulln(subdivpts);
				%trisurf(K,subdivpts(:,1),subdivpts(:,2),subdivpts(:,3))
		end
		ax = gca;
		ax.View = [69.9928   18.5834];
		axis equal tight
		
	case 2
		%subdivide a d-2 simplex (e.g. a line segment in 3D Cartesian coordinates)
		d = 3;
		pts = eye(d);
		pts(1,:) = [];
		pts(2,:) = ones(1,d);
		pts = normr(pts);
		
		% 		% subdivide a 3D line segment
		% 		pts = [0;1];
		% 		pts = [zeros(size(pts)) pts pts];
		
		if d == 3
			fig=figure;
			fig.Position=[368.5000  346.0000  560.0000  272.5000];
			tiledlayout(1,2)
			nexttile
			t=n2c(pts);
			plot3(t{:},'k-*')
			axis equal tight
			% 			xlim([-0.5 0.5])
			ax=gca;
			ax.View = [-117.2292   32.0689];
			
		elseif d == 4
			fig=figure;
			fig.Position=[214.5000  309.5000  909.5000  311.5000];
			% 			tiledlayout(1,2)
			nexttile
			[tmp,usv,zeropt] = proj_down(pts); %remove column of zeros
			t=n2c(tmp);
			triplot(1:3,t{:})
			axis equal tight
			% 			ax=gca;
			% 			ax.View = [101.5023   20.2621];
		end
		
		nint = 6;
		[subdivpts,TRI] = facet_subdiv(pts,nint,true);
		
		if d == 3
			nexttile
			t=n2c(subdivpts);
			plot3(t{:},'k*')
			hold on
			plot3(subdivpts(TRI,1),subdivpts(TRI,2),subdivpts(TRI,3)) %might depend on it being a closed polygon..
			axis equal tight
			% 			xlim([-0.5 0.5])
			ax=gca;
			ax.View = [-117.2292   32.0689];
			
		elseif d == 4
			for nint = 2:6
				[subdivpts,TRI] = facet_subdiv(pts,nint,true);
				
				nexttile
				[tmp,usv,zeropt] = proj_down(subdivpts);
				t=n2c(tmp);
				% 			K = delaunayn(tmp);
				triplot(TRI,t{:})
				axis equal tight
				% 			hold on
				% 			plot3(t{:},'*')
				% 			ax=gca;
				% 			ax.View = [101.5023   20.2621];
			end
			print -clipboard -dbitmap
		end
end


%------------------------------CODE GRAVEYARD------------------------------
%{
tmp = proj_down(subdivpts,1e-6,struct.empty,'nforce',1,'nforceQ',true); %remove column of zeros
%}