function rotvecplot(q0,q1,r,dt,scl,epsijk)
arguments
    q0(:,4) double = get_cubo()
    q1(:,4) double = get_cubo()
    r(:,3) double = rand(1,3);
    dt(1,1) double = 0.1
    scl(1,1) double = 1
    epsijk(1,1) double = 1
end

ax=gca;
ax.View = [45 38];
axis equal
lims = [-1 1];
xlim(lims);ylim(lims);zlim(lims);

r1 = rotslerp(r,q0,q1,dt,epsijk);
npts = size(r1,1);
t=n2c(1);
% paperfigure()
w = [0,0,0];
vecplot(scl,r(1,:),w,'b');
for i = 1:npts
    t = n2c(scl*r1(i,:));
    scatter3(t{:},'k')
    ax = vecplot(scl,r1(i,:),w,'b');
    drawnow
    pause(0.1)
    set(ax,'Visible','off')
end

vecplot(scl,r1(end,:),w,'r',"r'");