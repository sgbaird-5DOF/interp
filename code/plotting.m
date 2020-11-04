% PLOTTING  plotting script for interp5DOF paper
%split apply & find groups
% fname = 'gitID-76fca8c_uuID-f51500cd_paper-data.mat';
% fname = 'gitID-396aaa2_uuID-6816f860_paper-data2.mat';
% fname = 'gitID-f585733_uuID-edf2fcc7_paper-data2.mat';
fname = 'gitID-c67a123_uuID-18d21f26_set4.mat';
files = dir(fullfile('**',fname));
fpath = fullfile(files(1).folder,files(1).name);
load(fpath);
disp('file loaded')

% slurmQ = 0;

%saving
files = dir(fullfile('**','interp5DOF-paper','figures'));
folder = files(1).folder;

addpathdir({'paperfigure.m','dist-parity.mat','olm_octonion_list.txt',...
    'oct50000.mat','gmat2q.m','PGnames.mat','var_names.m'})

set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'defaultAxesFontSize',12)

%% parity plot
methodlist = {'pbary','gpr','idw','nn'};
datatypelist = {'brk','kim'};
sgtitleQ = false;
for datatype = datatypelist
    for nmeshpts = [388 10000 50000]
        %extract parity and IDs
        tbl3 = mdlparstbl(...
            ismember(mdlparstbl.method,methodlist) & ...
            ismember(mdlparstbl.datatype,datatype) & ...
            mdlparstbl.nmeshpts==nmeshpts,:);
        [G3,ID3] = findgroups(tbl3.method);
        parity3 = splitapply(@(x){x(1)},tbl3.parity,G3);
        
        %plotting
        multiparity(parity3,ID3,'hex','titleQ',true)
        
        %extra
        if sgtitleQ
            sgtitle(['nmeshpts = ' int2str(nmeshpts)])
        end
        
        savefigpng(folder,[char(datatype) 'parity' int2str(nmeshpts)]);
    end
end

%% errors
methodlist = {'pbary','gpr','idw','nn','avg'};

for datatype = datatypelist
    tbltmp = mdlparstbl(mdlparstbl.datatype==datatype,:);
    multixyplots(tbltmp,methodlist,'nmeshpts',{'rmse','mae'},2,1,'ymin',0,'lgdloc','southwest')
    %saving
    savefigpng(folder,[char(datatype) 'error'])
end

%% timing
methodlist = {{'pbary','gpr'},{'idw','nn'}};
tbl3 = mdlparstbl(mdlparstbl.ncores==12,:);
%{ 
whether to multiple pbary by # cores since uses parfor (although rest
might use multiple cores via vectorization..)?
%}
% ids = ismember(tbl3.method,'pbary');
% tbl3(ids,'runtime') = tbl3(ids,'runtime').*tbl3(ids,'cores');
[G3,ID3] = findgroups(tbl3.method);
multixyplots(mdlparstbl,methodlist,'nmeshpts',{'runtime'},2,1,'YScale','linear','yunits','s','lgdloc','north')

% stackedplot
%saving
savefigpng(folder,'runtime');

%% arclength vs. euclidean distance parity
% addpathdir({'get_five.m','cu2qu.m','GBfive2oct.m','gmat2q.m','distance-parity.mat'})
addpathdir({'dist-parity.mat','paperfigure.m'})
% seed = 10; %#ok<*UNRCH>
% rng(seed);
% npts = 10000;
% five = get_five(npts);
% o = GBfive2oct(five);
% pts = get_octpairs(o);
% pts = normr(pts);
% pd1 = pdist(pts).';
% pd2 = real(pdist(pts,@get_omega).');
% pd3 = real(pdist(pts,@get_alen).');

fname = 'dist-parity';
load([fname,'.mat'],'pd1','pd2','pd3')

paperfigure()
parityplot(pd1,pd3,'scatter','cscale','linear','xname','Euclidean Distance',...
    'yname','Arc Length','xunits','','yunits','rad','sz',20,'mkr','.','c',[0.6350, 0.0780, 0.1840])

%saving
% savefigpng(folder,fname);
fpath = fullfile(folder,fname);
print(fpath,'-dpng','-r300')

    
%     load('olm_pairwise_oct.mat','oct_new')
%     t = num2cell(oct_new,3);
%     t = cellfun(@(x) squeeze(x).',t.','UniformOutput',false);
%     oAB = vertcat(t{:});
%     oA = get_octpairs(oAB(:,1:8));
%     oB = get_octpairs(oAB(:,9:16));

%% octonion ensemble
addpathdir({'olm_octonion_list.txt','olm_pairwise_distances_cubic.mat'})
A = importdata('olm_octonion_list.txt');
olmoct = A.data;

rng(5)

nptslist = [1 2 5 10];
npts = max(nptslist);
numnpts = length(nptslist);

load('olm_pairwise_distances_cubic.mat','olmpairwisedistancescubic')
pd_olmchesser = table2array(olmpairwisedistancescubic);
pd7 = rad2deg(pd_olmchesser(:));

[pd4c, pd5c, pd6c, errmetrics] = deal(cell(numnpts,1));
oref = zeros(numnpts,8);
paperfigure(2,2,14.5);
j = 0;
for i = 1:npts
    if i == 1
        oref(i,:) = get_ocubo(1,'random',[],10);
    else
        oref(i,:) = get_ocubo(); %random reference octonion
    end
    olmoctsym = get_octpairs(olmoct,'oref',oref(i,:)); %symmetrize w.r.t. to oref
    pd4c{i} = squareform(pdist(normr(olmoctsym)));
    pd5c{i} = squareform(pdist(olmoctsym,@get_omega)); %pairwise distances
    pd6c{i} = squareform(pdist(olmoctsym,@get_alen));
    pd4c{i} = pd4c{i}(:);
    pd5c{i} = rad2deg(pd5c{i}(:)); % flatten
    pd6c{i} = pd6c{i}(:);
    
    pd4tmp = [pd4c{:}];
    pd5tmp = [pd5c{:}]; %catenate vectors of pairwise distances
    pd6tmp = [pd6c{:}];
    pd4 = min(pd4tmp,[],2);
    pd5 = min(pd5tmp,[],2);
    pd6 = min(pd6tmp,[],2);
    
    pd4tow = 2*180/pi*pd4;
    
    errmetrics{i} = get_errmetrics(pd4tow,pd7,'dispQ',true);
    errmetrics{i}.ksize = i;
    
    if any(i == nptslist)
        j = j+1;
        nexttile
        yrange = ['[1,' int2str(i),']'];
        yname = ['$\mathop{}_{\scriptstyle{\forall k \in' yrange '}}^{\rm{min}}2\frac{180}{\pi}|\hat{o}_{i,k}^{sym}-\hat{o}_{j,k}^{sym}|$'];
        parityplot(pd7,pd4tow,'hex','cscale','log','xname','$\omega$ (traditional)',...
            'yname',yname,'xunits','$^{\circ}$','yunits','$^{\circ}$')
        papertext(j,'xypos',[-0.2,1.0])
    end
end
% Euclidean: '$2\frac{180}{\pi}|o_1^{sym}-o_2^{sym}|$ (this work)'

klbls = cellfun(@char,num2cell(string(nptslist)),'UniformOutput',false);
lbls = strcat(klbls,'-');
lbls{end}(end) = []; %remove last '-'
savefigpng(folder,['dist-ensemble-k' lbls{:}])

paperfigure()
enstbl = struct2table(structvertcat(errmetrics{:})); %ensemble table
% G = findgroups(enstbl.ksize);
plot(enstbl.ksize,enstbl.rmse,'-*')
xlabel('$k_{max}$','Interpreter','latex')
ylabel('Error ($^{\circ}$)','Interpreter','latex')
hold on
plot(enstbl.ksize,enstbl.mae,'-*')
legend('RMSE','MAE','Interpreter','latex')
savefigpng(folder,['dist-ensemble-rmse-mae'])


%% distance parity (3D)
% pts = normr(rand(388,3)-0.5);
% pts2 = normr(rand(388,3)-0.5);
% pd1 = pdist(pts);
% pd3 = pdist(pts,@get_alen);
% parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','','title','388 3D points')
% figure
% t=n2c(pts);
% scatter3(t{:})
% axis equal

%% distance histogram and knn vs. meshpoints
addpathdir({'gmat2q.m','oct50000.mat'})
setpaperdefaults()
% distance histogram
% seed = 10;
% rng(seed);
% npts = 50000;
% five = get_five(npts);
% o = GBfive2oct(five);
% pts = get_octpairs(o);

load('oct50000.mat','pts');

paperfigure(1,2);
nexttile

nnhist(pts,'omega');
papertext(1)
% savefigpng(folder,['nnhist',int2str(npts)])

% knn distances
nexttile
[nnpts,D,mu,sigma] = get_knn(pts,'omega',100);
% fig = figure;
% fig.Position = [460.2 245 498.4 472.8];
% paperfigure()
errorbar(mu,sigma)
xlabel('k-th NN','Interpreter','latex','FontSize',12)
ylbl = 'average $\omega$ (deg)';
ylabel(ylbl,'Interpreter','latex','FontSize',12)
papertext(2)

% savefigpng(folder,['knnhist',int2str(npts)])
savefigpng(folder,['nnhist-knn-',int2str(npts)])

%% barycentric methods
proj_down_test(1)
savefigpng(folder,'bary-remove-deg')
proj_down_test(2)
savefigpng(folder,'bary-delaunay')

%% voronoi example
addpathdir('PGnames.mat')
toBPFZ_test(1)
% nexttile
% toBPFZ_test(1)
% nexttile
% toBPFZ_test(4)
savefigpng(folder,'voronoi')
% savefigpng(folder,'voronoi-4NN')

%% bary interp
sphbary_test()
savefigpng(folder,'bary-interp')


%% Prior Work Error Summary
S = load('prior-work-error-summary.mat');
ypred = S.lobpcg.ypred;
ytrue = S.lobpcg.ytrue;
rmse = 0.0277;
rmse2 = 0.0076;
f1 = nnz(ytrue < 0.9)/90000;
f2 = 1-f1;
f1*rmse+f2*rmse2

%% Arbitrary path through GB space




%% CODE GRAVEYARD
%{
%split apply & find groups
fname = 'gitID-6a89f2b_uuID-585a81b1_paper-data.mat';
files = dir(fullfile('**',fname));
fpath = fullfile(files(1).folder,files(1).name);
load(fpath);
[G1,ID1] = findgroups(mdlparstbl.(:,{'method','ndatapts'}));
for i = 1:length(G1)
    G = G1(i);
    GIDs = find(G1==G);
[G2,ID2] = findgroups(mdlparstbl.nmeshpts);
parity = splitapply(@(x){x},mdlparstbl.parity,G);
medRMSE = splitapply(@(x)median(x),mdlparstbl.rmse,G);
medMAE = splitapply(@(x)median(x),mdlparstbl.mae,G);

    plottype = 'scatter';
%     switch plottype
%         case 'hex'
%             h = hexscatter(ptmp.ytrue,ptmp.ypred,[0 1.5],[0 1.5],'charlbl',charlbl{i},'title',t{i},'cscale','log','cbnds',[1 500]);
%         case 'scatter'
            parityplot(ptmp.ytrue,ptmp.ypred,plottype,'charlbl',charlbl{i},'title',t{i},'scatterOpts',struct('MarkerEdgeAlpha',0.1))
%     end




fig=figure;
fig.Position = [418.6 194.6 644.8 555.2];

tiledlayout('flow','TileSpacing','compact','Padding','compact')
nIDs = length(ID3);
for i = 1:nIDs
    nexttile
    ptmp = parity3{i};
    charlblQ = T;
    titleQ = F;
    if charlblQ
        alphabet = ('a':'z').';
        chars = num2cell(alphabet(1:nIDs));
        chars = chars.';
        charlbl = strcat('(',chars,')');
    else
        charlbl = repelem({''},length(ID3));
    end
    if titleQ
        t = ID3;
    else
        t = repelem({''},length(ID3));
    end
    plottype = 'hex';
    parityplot(ptmp.ytrue,ptmp.ypred,plottype,'charlbl',charlbl{i},'title',t{i},'scatterOpts',struct('MarkerEdgeAlpha',0.1))
end


parity = splitapply(@(x){x(1)},mdlparstbl.parity,G1);


    % tbl1 = tbl(:,{'method','nmeshpts','ndatapts'});


    tbl2 = mdlparstbl(:,{'method','ndatapts'});
    [G2,ID2] = findgroups(ID1.method);
    nmeshpts1 = splitapply(@(x){x},mdlparstbl.nmeshpts,G2);


fig1 = figure;
fig1.Position = [489 90.6 330.4 672.4];
t = tiledlayout(2,1);
t1 = nexttile;
hold(t1,'on')
t1.XScale = 'log';
t2 = nexttile;
hold(t2,'on')
t2.XScale = 'log';

fig2 = figure;
fig2.Position = [489 90.6 330.4 672.4];
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
hold on
for method = methodlist
    tbl = mdlparstbl(ismember(mdlparstbl.method,method),:);
    [G1,ID1] = findgroups(tbl.nmeshpts);
    
    [npts,rmse,mae,runtime] = deal(tbl.nmeshpts,tbl.rmse,tbl.mae,tbl.runtime);
    
    nmeshpts = splitapply(@(x)x(1),npts,G1);
    medRMSE = splitapply(@median,rmse,G1);
    stdRMSE = splitapply(@std,rmse,G1);
    medMAE = splitapply(@median,mae,G1);
    stdMAE = splitapply(@std,mae,G1);
    medruntime = splitapply(@median,runtime,G1);
    stdruntime = splitapply(@std,runtime,G1);
    
    semilogx(t1,nmeshpts,medRMSE,'-o')
    xlabel('nmeshpts')
    ylabel('RMSE (J/m^2)')
    axis square
    hold off
    
    semilogx(t2,nmeshpts,medMAE,'-o')
    xlabel('nmeshpts')
    ylabel('MAE (J/m^2)')
    axis square
    hold off
%     semilogx(nmeshpts,medruntime,'-o')

    figure(fig2)
    hold on
    semilogx(nmeshpts,medruntime,'-o')
    xlabel('nmeshpts')
    ylabel('runtime (s)')
    axis square
    hold off
end


%% result parity plot
methodlist = {'sphgpr','gpr','sphbary','pbary','nn','avg'};
cellfun(@(type) median(mdlparstbl(strcmp(mdlparstbl.method,type),:).rmse),methodlist)
cellfun(@(type) median(mdlparstbl(strcmp(mdlparstbl.method,type),:).mae),methodlist)
uuid = 'ad8908c4';
tbltmp = mdltbl(strcmp(mdltbl.uuid,uuid),:);
barypars = tbltmp.barypars{1};
yactual = padcat(tbltmp.data.props(barypars.ids),tbltmp.data.props(barypars.ilist));
ypred = padcat(tbltmp.propOut{1}(barypars.ids),tbltmp.propOut{1}(barypars.ilist));
parityplot(yactual,ypred,'title','sphbary','legend',{'interp','nn'})


pd1 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts).';
pd2 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_omega).';
pd3 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_alen).';
parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','')

    oA = vertcat(oct_new(:,:,1:8));
    oB = vertcat(oct_new(:,:,9:16));


A = importdata('olm_octonion_list.txt');
    olmoct = A.data;
    olmoctsym = get_octpairs(olmoct);
    pd4 = squareform(pdist(olmoctsym));
    pd5 = squareform(pdist(olmoctsym,@get_omega));
    pd6 = squareform(pdist(olmoctsym,@get_alen));
    pd4 = pd4(:);
    pd5 = rad2deg(pd5(:));
    pd6 = pd6(:);
    load('olm_pairwise_distances_cubic.mat','olmpairwisedistancescubic')
    pd_olmchesser = table2array(olmpairwisedistancescubic);
    pd7 = rad2deg(pd_olmchesser(:));
    figure
    parityplot(pd5,pd7,'hex','cscale','linear','xname','\omega (Johnson)','yname','\omega (Chesser)','xunits','deg','yunits','deg')


%}