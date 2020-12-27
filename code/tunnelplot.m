function tunnelplot(mdl,A,B,n,nv)
arguments
    mdl = []
    A = []
    B = []
    n(1,1) double = 300
    nv.nnQ(1,1) logical = true
    nv.nnQ2(1,1) logical = true
    nv.brkQ(1,1) logical = true
end

nnQ = nv.nnQ;
nnQ2 = nv.nnQ2;
brkQ = nv.brkQ;

% unpack
if isfield(mdl,'gprMdl')
    if ~isempty(mdl.gprMdl)
        gprMdl = mdl.gprMdl;
    end
else
    gprMdl = mdl.cgprMdl;
end
ppts = mdl.mesh.ppts;
ppts2 = mdl.data.ppts;

%% pairwise distance
npts = size(ppts,1);
if npts > 20000
    ids = 1:20000;
else
    ids = 1:npts;
end

pts = mdl.mesh.pts;
if isempty(A) && isempty(B)
    pd = squareform(pdist(ppts(ids,:)));
    [mx,id] = max(pd,[],'all','linear');
    
    % max dimension
    indlen = length(ids);
    [row,col] = ind2sub([indlen,indlen],id);
    A = pts(row,:);
    B = pts(col,:);    
end
% a = ppts2(row,:);
% b = ppts2(col,:);
% a = [-0.204792525808597,-0.0645669209001381,-0.107955224738106,0.0381947287488522,0.0993759868054296,-0.0225722216770392,0.0189130134087235];
% b = [0.287494800345075,-0.00106803759647395,0.0779976779290583,-0.0258833561524949,-0.127050434093774,-0.00810139755838871,0.0258975652026086];

A = sqrt(2)*normr(A);
B = sqrt(2)*normr(B);
% pts2 = mdl.data.pts;
% A = pts2(row,:);
% B = pts2(col,:);
% A = proj_up(a,mdl.usv);
% B = proj_up(b,mdl.usv);

%% distance calculations/coordinate interpolation
d1 = get_omega(A,B);

arcpts = normr(OSLERP(A,B,d1,n));
d2 = d1;

% t = arrayfun(@(i) linspace(A(i),B(i),n).',1:size(A,2),'UniformOutput',false); %1D interpolation across each dimension
% tpts = [t{:}];
% tpts = normr(tpts);
% t = n2c(tpts);
% arcpts = interparc(n,t{:},'spline');
% arcpts = 1/sqrt(2)*sqrt2norm(arcpts);
% [d2,seg] = arclength(t{:},'spline');
% d2 = 2*d2; %convert to omega

% ids = knnsearch(pts2,arcpts);
% nnpts = uniquetol(pts2(ids,:),'ByRows',true);
% t = n2c(nnpts);
% arcpts = interparc(n,t{:},'spline');
% arcpts = 1/sqrt(2)*sqrt2norm(arcpts);

arcppts = proj_down(arcpts,1e-4,mdl.usv,'zeroQ',mdl.zeroQ);

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
[~,nnd,~,~,ids] = get_knn(ppts,'norm',K,'Y',arcppts);
tprednn = cellfun(@(id) propList(id),n2c(ids),'UniformOutput',false);
% tprednn = propList(dsearchn(ppts,arcppts));

% errorbar(tpred,tsd)
x = linspace(0,rad2deg(d2),n);
hold on
ax = cell(3,1);
sz2 = zeros(n,K);

if nnQ2
    for i = 1:K
        sz = -rescale(nnd{i},-50,-2);
        %     sz2 = 2*rad2deg(-nnd{i}+2*min(nnd{i}));
        sz2(:,i) = 2*rad2deg(nnd{i});
        scatter(x,tprednn{i},sz,sz2(:,i),'o')
    end
end
i = 1;
[lgdax,lgdlbl]=deal({});
if nnQ
    ax{i} = plot(x,tprednn{1},'k--','LineWidth',1); %#ok<*UNRCH>
    lgdlbl = [lgdlbl,'$\overline{AB}$ 1st-NN'];
    lgdax = [lgdax,ax{i}];
    i = i+1;
end
axtmp=shadedErrorBar(x,tpred,tsd);
ax{i}=axtmp.mainLine;
ax{i}.LineWidth = 2.5;

pA = normr(arcpts(:,1:4));
pB = normr(arcpts(:,5:8));

lgdax = [lgdax,ax{i}];
lgdlbl = [lgdlbl,'$\overline{AB}$ GPR'];
i = i+1;
if brkQ
    brk = GB5DOF_setup(pA,pB,[0 0 1],'Ni',1);
    ax{i} = plot(x,brk,'r-','LineWidth',1);
    lgdax = [lgdax,ax{i}];
    lgdlbl = [lgdlbl,'$\overline{AB}$ BRK'];
end

if nnQ2
    cb = colorbar;
    cb.Label.String = 'distance to $\overline{AB}$ $(^\circ)$';
    cb.Label.Interpreter = 'latex';
    caxis([min(sz2,[],'all'),max(sz2,[],'all')])
    colormap jet
end
% kstr = num2cell(strcat(string(num2cell(1:K)),'-NN'));

legend(lgdax,lgdlbl,'Location','northoutside','Interpreter','latex')

xlabel('$\overline{AB}(t)$ $(^\circ)$','Interpreter','latex')
ylabel('$GBE (J m^{-2})$','Interpreter','latex')
axis square tight
end

%% CODE GRAVEYARD
%{
%% load data
% Sstr = 'gpr10000_gitID-9a0720b_puuID-51432301_rohrer-Ni-regularization.mat';
% S2str = 'nn10000_gitID-9a0720b_puuID-d2c68f10_rohrer-Ni-regularization.mat';

% Sstr = 'gpr1000_gitID-9a0720b_puuID-a92bb17f_rohrer-Ni-lamb1000.mat';
% S2str = 'nn1000_gitID-9a0720b_puuID-441cf59f_rohrer-Ni-lamb1000.mat';

% Sstr = 'gpr10000_gitID-9a0720b_puuID-2002d1ea_rohrer-Ni-lamb500m.mat';
% Sstr = 'gpr10000_gitID-9a0720b_puuID-2002d1ea_rohrer-Ni-lamb300m.mat';

% S = load(Sstr);
% S2 = load(S2str);

% mdl = S.mdl;
% gprMdl = mdl.gprMdl;
%}