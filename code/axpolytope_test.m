% AXPOLYTOPE_TEST generate vertices on the axes of a polytope in n-D

d = 3;
axpPts = axpolytope(d,1:2);

if d == 3
	axpK = convhulln(axpPts);
	t = num2cell(axpPts,1);
	trisurf(axpK,t{:})
	axis equal
end