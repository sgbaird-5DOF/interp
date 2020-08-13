%misFZcon test
clear; close all
[A,b] = misFZcon();

dlist = cprnd(10000,A,b);

plotFZrodriguez_vtx();
t=n2c(dlist);
plot3(t{:},'*')

qlist = disorientation(rod2q(dlist),'cubic');

geometry = findgeometry(qlist,1e-2);
sum(~ismember(geometry,{'interior','exterior'}))
dlist = q2rod(qlist);

misFZ = inmisFZ(dlist);

sum(misFZ)