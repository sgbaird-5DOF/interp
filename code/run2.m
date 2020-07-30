% run test

NV = {'o2addQ',false,'method','pairwise','wtol',1e-3};

mesh.fname = 'meshtest';
data.fname = 'datatest';

disp('SETUP')
mesh.pts = get_octpairs(get_ocubo(50,'random',[],15),mesh.fname,NV{:}); %output seems to be the same
data.pts = get_octpairs(get_ocubo(50,'random',[],20),data.fname,NV{:});

%project mesh and data together
[a,usv] = proj_down([mesh.pts;data.pts],1e-3);

if size(a,2) < size(mesh.pts,2)
	data.pts = proj_down(data.pts,1e-6,usv);
	mesh.pts = proj_down(mesh.pts,1e-6,usv);
end

mesh.K = sphconvhulln(mesh.pts,true); %sphconvhulln output was different

intfacetIDs = intersect_facet(mesh.pts,mesh.K,data.pts,1e-6,true);

disp(['# non-intersecting facets = ' ...
	int2str(sum(cellfun(@isempty,intfacetIDs))) '/' int2str(length(intfacetIDs))])

disp('PROPERTIES')
mesh.props = GB5DOF_setup(GBoct2five(mesh.pts));
data.props = GB5DOF_setup(GBoct2five(data.pts));


disp('INTERPOLATION')
[nndistList,databary] = get_interp(mesh,data);