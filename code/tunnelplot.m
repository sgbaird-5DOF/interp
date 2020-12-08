function tunnelplot(interpfn,a,b,n)
arguments
    interpfn = []
    a = []
    b = []
    n = []
end

%% load data
% Sstr = 'gpr10000_gitID-9a0720b_puuID-51432301_rohrer-Ni-regularization.mat';
% S2str = 'nn10000_gitID-9a0720b_puuID-d2c68f10_rohrer-Ni-regularization.mat';

% Sstr = 'gpr1000_gitID-9a0720b_puuID-a92bb17f_rohrer-Ni-lamb1000.mat';
% S2str = 'nn1000_gitID-9a0720b_puuID-441cf59f_rohrer-Ni-lamb1000.mat';

% Sstr = 'gpr10000_gitID-9a0720b_puuID-2002d1ea_rohrer-Ni-lamb500m.mat';
Sstr = 'gpr10000_gitID-9a0720b_puuID-2002d1ea_rohrer-Ni-lamb300m.mat';

S = load(Sstr);
% S2 = load(S2str);

% unpack
mdl = S.mdl;
gprMdl = mdl.gprMdl;
% ppts = mdl.mesh.ppts;
ppts2 = mdl.data.ppts;

%% pairwise distance
pd = squareform(pdist(ppts2));
[mx,id] = max(pd,[],'all','linear');
npts2 = size(ppts2,1);

%% max dimension
[row,col] = ind2sub([npts2,npts2],id);
a = ppts2(row,:);
b = ppts2(col,:);
n = 200;
% pts = mdl.mesh.pts;
pts2 = mdl.data.pts;
A = pts2(row,:);
B = pts2(col,:);

%% distance calculations/coordinate interpolation
d1 = get_omega(A,B);
t = arrayfun(@(i) linspace(A(i),B(i),n).',1:size(A,2),'UniformOutput',false);
tpts = [t{:}];
t = n2c(tpts);
arcpts = interparc(n,t{:},'spline');
arcpts = 1/sqrt(2)*sqrt2norm(arcpts);
[d2,seg] = arclength(t{:},'spline');
d2 = 2*d2;

% ids = knnsearch(pts2,arcpts);
% nnpts = uniquetol(pts2(ids,:),'ByRows',true);
% t = n2c(nnpts);
% arcpts = interparc(n,t{:},'spline');
% arcpts = 1/sqrt(2)*sqrt2norm(arcpts);

arcppts = proj_down(arcpts,1e-4,mdl.usv,'zeroQ',false);

%% property interpolation
[tpred,tsd,tint] = predict(gprMdl,arcppts);
scl = 1;
tpred = tpred*scl;
tsd = tsd*scl;
% mdl2 = S2.mdl;
% mdlcmd2 = mdl2.mdlcmd;
propList = mdl.mesh.props*scl;

%%
K = 6;
[~,nnd,~,~,ids] = get_knn(ppts2,'norm',K,'Y',arcppts);
tprednn = cellfun(@(id) propList(id),n2c(ids),'UniformOutput',false);
% tprednn = propList(dsearchn(ppts,arcppts));

% errorbar(tpred,tsd)
x = linspace(0,rad2deg(d2),n);
hold on
ax = cell(3,1);
sz2 = zeros(n,K);
for i = 1:K
    sz = -rescale(nnd{i},-50,-2);
%     sz2 = 2*rad2deg(-nnd{i}+2*min(nnd{i}));
    sz2(:,i) = 2*rad2deg(nnd{i});
    scatter(x,tprednn{i},sz,sz2(:,i),'o')
end
ax{1} = plot(x,tprednn{1},'k--','LineWidth',1);

ax{2}=shadedErrorBar(x,tpred,tsd);
ax{2}.mainLine.LineWidth = 2.5;

pA = arcpts(:,1:4);
pB = arcpts(:,5:8);
brk = GB5DOF_setup(pA,pB,[0 0 1],1);

ax{3} = plot(x,brk,'r-','LineWidth',1);

% kstr = num2cell(strcat(string(num2cell(1:K)),'-NN'));

legend([ax{1},ax{2}.mainLine,ax{3}],{'$\overline{AB}$ 1st-NN','$\overline{AB}$ GPR','$\overline{AB}$ BRK'},'Location','northoutside','Interpreter','latex')

xlabel('$\overline{AB}(t)$ $(^\circ)$','Interpreter','latex')
ylabel('$GBE (J m^{-2})$','Interpreter','latex')
cb = colorbar;
cb.Label.String = 'distance to $\overline{AB}$ $(^\circ)$';
cb.Label.Interpreter = 'latex';
caxis([min(sz2,[],'all'),max(sz2,[],'all')])
colormap jet
axis square tight
end
