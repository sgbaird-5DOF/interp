clear; close all;

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

newpts = facet_subdiv(pts,nint);

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
axis equal tight