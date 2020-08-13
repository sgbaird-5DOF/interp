%misFZcon test
clear; close all
[A,b] = misFZcon();

pts = cprnd(10000,A,b);

plotFZrodriguez_vtx();
t=n2c(pts);
plot3(t{:},'*')