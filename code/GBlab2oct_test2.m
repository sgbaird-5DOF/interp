%troubleshooting and visualizing conversion from lab coordinates to qm/nA pairs
clear; close all
% First, to get a good visualization of axes

scl = 1;
x = [1 0 0];
y = [0 1 0];
z = [0 0 1];
w = [0 0 0];
paperfigure()
hold on
mkr = 'k';
quivplot(scl,x,y,z,w,mkr);
lims=[-1 1];
xlim(lims)
ylim(lims)
zlim(lims)
axis equal off
view([45 38])

rng(15)
v0 = get_cubo();
% v0 = [1 0 0 0];
% v1 = [0 0 0 1];
v1 = get_cubo();
dt = 0.1;

rotquivplot(v0,v1,dt,scl);
