%slight rotation matrix test
clear; close all

d = 3;

R = slightRot(d);

% npts = 10;
% pts = normr(rand(npts,d));

% pts = eye(d);
pts = axpolytope(d);

if d == 3
	nexttile
	if size(pts,1) == d
		K = 1:d;
	else
		K = convhulln(pts);
	end
	t = num2cell(pts,1);
	trisurf(K,t{:},'FaceColor','none','EdgeColor','black');
	
	hold on
	rotpts = (R*pts.').';
	t = num2cell(rotpts,1);
	trisurf(K,t{:},'FaceColor','none','EdgeColor','red');
	axis equal
	
	alpha(0.5);
end


check = ismembertol(rotpts,axpolytope(d),1e-6,'ByRows',true);
msg = ['number of points that still fall on axpolytope: ' int2str(sum(check))];
disp(msg);
% assert(sum(check) == 0,msg)