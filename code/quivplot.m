function quivplot(scl)
arguments
   scl(1,1) double = 1
end
% QUIVPLOT  plot three quivers in the x-hat, y-hat, and z-hat directions
x = [1,0,0]*scl;
y = [0,1,0]*scl;
z = [0,0,1]*scl;

w = [0,0,0];
hold on
quiver3(w,w,w,x,y,z,1,'linewidth',1)
text(x(1),x(2)+0.01,x(3),'x','FontWeight','bold')
text(y(1),y(2),y(3)+0.01,'y','FontWeight','bold')
text(z(1)+0.01,z(2),z(3),'z','FontWeight','bold')

end