%hypercube test
clear; close all

d = 3;

[pts,K] = hypercube(d);

%plotting
if d == 3
	ptscell = num2cell(pts,1);
	trisurf(K,ptscell{:})
	axis equal
end