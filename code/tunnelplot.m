function [tpredlist,tsdlist,propList,methodlist,A,B] = tunnelplot(mdls,A,B,n,nv)
arguments
    mdls = []
    A = []
    B = []
    n(1,1) double = 300
    nv.nnQ(1,1) logical = true
    nv.nnQ2(1,1) logical = true
    nv.brkQ(1,1) logical = true
    nv.lgdloc char = 'best' %'northoutside'
    nv.tpredlist = []
    nv.tsdlist = []
    nv.cilist = []
    nv.propList = []
    nv.methodlist = []
    nv.covmat = []
    nv.kfntmp = []
    nv.kfntmp2 = []
    nv.gprMdl2 = []
    nv.nsamp = []
    nv.extend = 0
end
% TUNNELPLOT  plot points along a 1D arc and NNs if specified (like going through a tunnel)
nnQ = nv.nnQ;
nnQ2 = nv.nnQ2;
brkQ = nv.brkQ;
extend = nv.extend;

if ~isempty(mdls) && ~iscell(mdls) && isscalar(mdls)
    mdls = {mdls};
end

if ~isempty(mdls)
    ppts = mdls{1}.mesh.ppts;
    pts = mdls{1}.mesh.pts;
end

% points to plot between
if isempty(A) && isempty(B)
    [A,B] = getAB(pts); %get points that are far apart
end

A = sqrt(2)*normr(A);
B = sqrt(2)*normr(B);

%% distance calculations/coordinate interpolation
d1 = get_omega(A,B);

slerptype = 'interparc'; %'oslerp', 'interparc'
switch slerptype
    case 'oslerp'
        arcpts = normr(get_octpairs(OSLERP(A,B,d1,n,extend),'dispQ',false));
        d2 = d1;
    case 'interparc'
        dx = B-A;
        Apre = sqrt2norm(-extend*dx+A,'quat','warnQ',false);
        Bpost = sqrt2norm(extend*dx+B,'quat','warnQ',false);
        t = arrayfun(@(i) linspace(Apre(i),Bpost(i),n).',1:size(A,2),'UniformOutput',false); %1D interpolation across each dimension
        tpts = [t{:}];
        tpts = 1/sqrt(2)*[normr(tpts(:,1:4)),normr(tpts(:,5:8))];
        t = n2c(tpts);
        arcpts = interparc(n,t{:},'spline');
        arcpts = normr(get_octpairs(sqrt2norm(arcpts,'quat'),[],1,'oref',mdls{1}.oref,'dispQ',false));
        %         arcpts = 1/sqrt(2)*sqrt2norm(arcpts);
        [d2,seg] = arclength(t{:},'spline');
%         d2 = 2*d2; %convert to omega
end

%% property interpolation
if isempty(mdls)
    nmdl = length(nv.tpredlist);
else
    nmdl = length(mdls);
end
if (isempty(nv.tpredlist) && isempty(nv.tsdlist) && isempty(nv.propList) && ...
        isempty(nv.methodlist)) || any(strcmp(mdls{1}.method,{'gprmix','egprm','egpr','gprm'}))
    [tpredlist,tsdlist,cilist,propList,methodlist] = deal(cell(nmdl,1));
    for i = 1:nmdl
        mdl = mdls{i};
        methodlist{i} = mdl.method;
        if mdl.projQ
            arcppts = proj_down(arcpts,mdl.projtol,mdl.usv,'zeroQ',mdl.zeroQ);
            %         ptstmp = proj_down([arcpts;mdl.mesh.pts],mdl.projtol,'zeroQ',mdl.zeroQ);
            %         arcppts = ptstmp(1:size(arcpts,1),:);
            %         ppts = ptstmp((size(arcpts,1)+1):end,:);
        else
            arcppts = arcpts;
        end
        switch mdl.method
            case 'pbary'
                mdl.data.pts = arcpts;
                mdl.data.ppts = arcppts;
                mdl.data.npts = size(arcpts,1);
                mdl.data.props = GB5DOF_setup(arcpts(:,1:4),arcpts(:,5:8),[0 0 1],'Ni',1);
                [intfacetIDs,databary,klist] = intersect_facet(mdl.mesh.ppts,mdl.mesh.sphK,arcppts,mdl.inttol,'nnMax',mdl.nnMax);
                tpredlist{i} = get_interp(mdl.mesh,mdl.data,intfacetIDs,mdl.barytype,mdl.barytol);
            case 'gpr'
                if isfield(mdl,'gprMdl')
                    gprMdl = mdl.gprMdl;
                else
                    gprMdl = mdl.cgprMdl;
                end
                [tpredlist{i},tsdlist{i},cilist{i}] = predict(gprMdl,arcppts);
                
                %     out = exec_argfn(mdl.mdlcmd,mdl,{'tpred','tsd'});
                %     tpred = out.tpred;
                %     tsd = out.tsd;
            case 'gprmix'
                alpha = 0.05;
                %                 [tpredlist{i},tsdlist{i}] = predict(mdl.gprMdl,arcppts);
                %                 if isempty(nv.kfntmp)
                %                     kfn = mdl.gprMdl.Impl.Kernel.makeKernelAsFunctionOfXNXM(mdl.gprMdl.Impl.ThetaHat);
                %                     if isempty(nv.covmat)
                %                         covmat = kfn(arcppts,arcppts);
                %                     else
                %                         covmat = nv.covmat;
                %                     end
                %                 else
                [tpredlist{1},tsdlist{1},cilist{1}] = gprmix(mdl,[],[],arcppts,'gprMdl2',nv.gprMdl2);
                kfn = kmix(nv.kfntmp,nv.kfntmp2);
                a = sigfn(tpredlist{1});
                covmat = kfn(arcppts,arcppts,a);
                covmat = nearestSPD(covmat);
                %                 end
                %                 [tpredlist{i},covmat,ci] = predictExactWithCov(mdl.gprMdl.Impl,arcppts,alpha);
                %                 tsdlist{i} = sqrt(diag(covmat));
                %                 T = cholcov(covmat);
                %                 tpredlist{i} = tpredlist{i} + T'*randn(n,nsamp);
                if ~isempty(nv.nsamp)
                    tpredlist{i} = mvnrnd(tpredlist{i}.',covmat,nv.nsamp);
                end
            case 'idw'
                tpredlist{i} = idw(mdl.mesh.ppts,arcppts,mdl.mesh.props,mdl.r,mdl.L);
                %     out = exec_argfn(mdl.mdlcmd,mdl,{'tpred'});
                %     tpred = out.tpred;
            case 'nn'
                tpredlist{i} = mdl.mesh.props(dsearchn(mdl.mesh.ppts,arcppts));
                
            case {'egprm','egpr','gprm'}
                K = mdl.K;
                mixQ = mdl.mixQ;
                [~,tpredlist{i},tsdlist{i},~,cilist{i}] = ...
                    egprm([],[],[],[],[],K,'mixQ',mixQ,...
                    'o2',sqrt2norm(arcpts,'quat'),'egprmMdl',mdl,...
                    'dispQ',false,'KdispQ',false,'egprmDispQ',false);
        end
        scl = 1;
        tpredlist{i} = tpredlist{i}*scl;
        tsdlist{i} = tsdlist{i}*scl;
        cilist{i} = cilist{i}*scl;
        
        % mdl2 = S2.mdl;
        % mdlcmd2 = mdl2.mdlcmd;
        propList{i} = mdl.mesh.props*scl;
    end
else
    tpredlist = nv.tpredlist;
    tsdlist = nv.tsdlist;
    cilist = nv.cilist;
    propList = nv.propList;
    methodlist = nv.methodlist;
end

%%
% tprednn = propList(dsearchn(ppts,arcppts));

% errorbar(tpred,tsd)
preExtend = rad2deg(d2);
start = -extend*rad2deg(d2);
finish = (1+extend)*rad2deg(d2);
x = linspace(start,finish,n);
% x = linspace(0,rad2deg(d2),n);
hold on
ax = cell(3,1);

if nnQ2
    K = 6;
    [~,nnd,~,~,ids] = get_knn(ppts,'norm',K,'Y',arcppts);
    tprednn = cellfun(@(id) propList{i}(id),n2c(ids),'UniformOutput',false);
    sz2 = zeros(n,K);
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

if brkQ
    pA = normr(arcpts(:,1:4));
    pB = normr(arcpts(:,5:8));
    brk = GB5DOF_setup(pA,pB,[0 0 1],'Ni',1);
    ax{i} = plot(x,brk,'k-','LineWidth',3);
    lgdax = [lgdax,ax{i}];
    lgdlbltmp = 'BRK';
    lgdlbl = [lgdlbl,lgdlbltmp];
    i = i+1;
end

errbartype = 'ci'; %'ci' or 'sd'
for j = 1:nmdl
    method = methodlist{j};
    switch method
        case 'gpr'
            switch errbartype
                case 'sd'
                    axtmp=shadedErrorBar(x,tpredlist{j},tsdlist{j},'lineProps','bo');
                case 'ci'
                    axtmp=shadedErrorBar(x,tpredlist{j},abs(cilist{j}-tpredlist{j}),'lineProps','bo');
            end
            ax{j}=axtmp.mainLine;
            ax{i}.LineWidth = 0.25;
            ax{j}.MarkerSize = 3;
        case {'gprmix','egprm','egpr','gprm'}
            if isempty(nv.nsamp)
                switch errbartype
                    case 'sd'
                        axtmp=shadedErrorBar(x,tpredlist{j},tsdlist{j},'lineProps','bo');
                    case 'ci'
                        axtmp=shadedErrorBar(x,tpredlist{j},abs(cilist{j}-tpredlist{j}),'lineProps','bo');
                end
                ax{j}=axtmp.mainLine;
                ax{i}.LineWidth = 0.25;
                ax{j}.MarkerSize = 3;
            else
                mkr = '-o';
                ax{j} = plot(x,tpredlist{j},mkr,'MarkerSize',3,'LineWidth',0.25);
            end
        otherwise
            switch method
                case 'pbary'
                    mkr = 'ro';
                case 'nn'
                    mkr = 'mo';
                case 'idw'
                    mkr = 'go';
            end
            ax{j} = plot(x,tpredlist{j},mkr,'MarkerSize',3,'LineWidth',0.25);
    end
    lgdax = [lgdax,ax{j}];
    switch method
        case 'pbary'
            methodtxt = 'Barycentric';
        case {'gpr','egprm','gprmix','egpr','gprm'}
            methodtxt = upper(char(method));
            switch errbartype
                case 'sd'
                    methodtxt = [methodtxt,'($\pm\sigma$)']; %#ok<AGROW>
                case 'ci'
                    methodtxt = [methodtxt,'($\pm$CI)']; %#ok<AGROW>
            end
        otherwise
            methodtxt = upper(char(method));
    end
    lgdlbl = [lgdlbl,methodtxt]; %#ok<AGROW>
    i = i+1;
end

if nnQ2
    cb = colorbar;
    cb.Label.String = 'distance to $\overline{AB}$ $(^\circ)$';
    cb.Label.Interpreter = 'latex';
    caxis([min(sz2,[],'all'),max(sz2,[],'all')])
    colormap jet
end

if extend ~= 0
    xline(0,'--')
    xline(preExtend,'--')
    %     lgdax = [lgdax,xline(0,'--')];
    %     lgdax = [lgdax,xline(preExtend,'--')];
    %     lgdlbl = [lgdlbl,'A','B'];
end

% kstr = num2cell(strcat(string(num2cell(1:K)),'-NN'));
if ~strcmp(mdl.method,'gprmix')
    legend(lgdax,lgdlbl,'Location',nv.lgdloc,'Interpreter','latex')
end

xlabel('$\overline{AB}(t)$ $(^\circ)$','Interpreter','latex')
ylabel('GBE $(J m^{-2})$','Interpreter','latex')
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


% a = ppts2(row,:);
% b = ppts2(col,:);
% a = [-0.204792525808597,-0.0645669209001381,-0.107955224738106,0.0381947287488522,0.0993759868054296,-0.0225722216770392,0.0189130134087235];
% b = [0.287494800345075,-0.00106803759647395,0.0779976779290583,-0.0258833561524949,-0.127050434093774,-0.00810139755838871,0.0258975652026086];

% pts2 = mdl.data.pts;
% A = pts2(row,:);
% B = pts2(col,:);
% A = proj_up(a,mdl.usv);
% B = proj_up(b,mdl.usv);


% % unpack
% if isfield(mdl,'gprMdl')
%     if ~isempty(mdl.gprMdl)
%         gprMdl = mdl.gprMdl;
%     end
% else
%     gprMdl = mdl.cgprMdl;
% end


% ids = knnsearch(pts2,arcpts);
% nnpts = uniquetol(pts2(ids,:),'ByRows',true);
% t = n2c(nnpts);
% arcpts = interparc(n,t{:},'spline');
% arcpts = 1/sqrt(2)*sqrt2norm(arcpts);


    %     lgdlbltmp = ['$\overline{AB}$ ' methodtxt];
    %     lgdlbltmp = methodtxt;

    %     lgdlbltmp = '$\overline{AB}$ BRK';

% ppts2 = mdl.data.ppts;

%}