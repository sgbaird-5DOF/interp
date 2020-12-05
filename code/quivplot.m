function ax = quivplot(scl,x,y,z,w,mkr,xlbl,ylbl,zlbl)
arguments
   scl(1,1) double = 1
   x(1,3) double = [1,0,0]
   y(1,3) double = [0,1,0]
   z(1,3) double = [0,0,1]
   w(1,3) double = [0,0,0]
   mkr char = 'k'
   xlbl = 'x';
   ylbl = 'y';
   zlbl = 'z';
end
% QUIVPLOT  plot three quivers in the x-hat, y-hat, and z-hat directions
xyz = scl*[x;y;z];
xcomp = xyz(:,1);
ycomp = xyz(:,2);
zcomp = xyz(:,3);
w = w.';

hold on
ax(1) = quiver3(w,w,w,xcomp,ycomp,zcomp,0,mkr,'linewidth',1,'Autoscale','off');
ax(2) = text(x(1),x(2)+0.05,x(3),xlbl,'FontWeight','bold');
ax(3) = text(y(1),y(2),y(3)+0.05,ylbl,'FontWeight','bold');
ax(4) = text(z(1)+0.05,z(2),z(3),zlbl,'FontWeight','bold');

end

%% CODE GRAVEYARD
%{
% x = x.'*scl;
% y = y.'*scl;
% z = z.'*scl;
%}