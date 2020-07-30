%knninterp test
clear; close all;

seed = 10;
rng(seed);

d = 3;

%mesh
nmeshpts = 100;
orthoPts = orthoplex(d);
ids = find(all(orthoPts >= 0, 2));
orthoPts = orthoPts(ids,:);
meshpts = normr([orthoPts; rand(nmeshpts,d)]);
% meshprops = (1:nmeshpts)/nmeshpts*10 + rand(1,nmeshpts);

randpoly = rand(d,1);
meshprops = meshpts*randpoly;

%data
ndatapts = 10;
datapts = normr(rand(ndatapts,d));
dataprops = datapts*randpoly;

interpvals = knninterp(meshpts,meshprops,datapts);

disp(['RMSE == ' num2str(immse(interpvals,dataprops))])

if d == 3
	t = num2cell(meshpts,1);
	K = sphconvhulln(meshpts);
	trisurf(K,t{:})
end