function rotquivplot(q0,q1,dt,scl,epsijk,nv)
arguments
    q0(:,4) double = get_cubo()
    q1(:,4) double = get_cubo()
    dt(1,1) double = 0.1
    scl(1,1) double = 1
    epsijk(1,1) double = 1
    nv.view = [45 38];
end
ax=gca;
ax.View = nv.view;
lims = [-1 1];
axis equal
xlim(lims);ylim(lims);zlim(lims);

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];
w = [0 0 0];
[x1,y1,z1] = rotslerpax(scl*x,scl*y,scl*z,q0,q1,dt,epsijk);
npts = size(x1,1);
% t=n2c(x1);
% paperfigure()
quivplot(scl,x1(1,:),y1(1,:),z1(1,:),w,'b');
for i = 1:npts
    t = n2c(scl*x1(i,:));
    scatter3(t{:},'k')
    t = n2c(scl*y1(i,:));
    scatter3(t{:},'r')
    t = n2c(scl*z1(i,:));
    scatter3(t{:},'b')
    ax = quivplot(scl,x1(i,:),y1(i,:),z1(i,:),w,'b');
    drawnow
    pause(0.1)
    set(ax,'Visible','off')
end

quivplot(scl,x1(end,:),y1(end,:),z1(end,:),w,'r',"x'","y'","z'");

%% CODE GRAVEYARD
%{
om = qu2om(q0,epsijk);
x = om(1,:);
y = om(2,:);
z = om(3,:);
%}