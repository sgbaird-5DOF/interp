%project facet to hyperplane test
% Dependencies:
%		normr.m (shadows built-in MATLAB version from Computer Vision
%		Toolbox, should work similarly)

%% setup
d = 7; %dimensionality
vertices = rand(d); %rows of vertices defining facet on hypersphere

vertices = normr(vertices); %normalize
nvec = mean(vertices); % point that will define tangent hyperplane
nvec = normr(nvec);

%% projection
newvertices = projfacet2hyperplane(nvec,vertices)

for i = 1:size(vertices,1)
vtx = newvertices(i,:);
	vtxnorm(i) = norm(vtx); %check norm magnitudes (should no longer be 1)
end

disp('new vertex norms = ')
disp(' ')
disp(vtxnorm.')

