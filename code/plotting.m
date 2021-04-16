% PLOTTING  plotting script for interp5DOF paper
%split apply & find groups
% fname = 'gitID-76fca8c_uuID-f51500cd_paper-data.mat';
% fname = 'gitID-396aaa2_uuID-6816f860_paper-data2.mat';
% fname = 'gitID-f585733_uuID-edf2fcc7_paper-data2.mat';
% fname = 'gitID-c67a123_uuID-18d21f26_set4.mat';
% fname = 'gitID-014bf70_uuID-3ed9cba0_paper-data3.mat';
% fname = 'gitID-6ede824_uuID-1cf78415_paper-data5.mat';
fname = 'gitID-0055bee_uuID-475a2dfd_paper-data6.mat';
files = dir(fullfile('**',fname));
fpath = fullfile(files(1).folder,files(1).name);
load(fpath);
disp('file loaded')

% slurmQ = 0;

%%
%saving
files = dir(fullfile('**','interp5DOF-paper','figures'));
figfolder = files(1).folder;
files = dir(fullfile('**','interp5DOF-paper','tables'));
tblfolder = files(1).folder;

addpath(genpath('.'))
setlatex()

%% parity plot
methodlist = {'pbary','gpr','idw','nn'};
% datatypelist = {'brk','kim'};
datatypelist = {'brk'};
sgtitleQ = false;
for datatype = datatypelist
    for ninputpts = [388 10000 50000]
        %extract parity and IDs
        tbl3 = mdlparstbl(...
            ismember(mdlparstbl.method,methodlist) & ...
            ismember(mdlparstbl.datatype,datatype) & ...
            mdlparstbl.ninputpts==ninputpts,:);
        [G3,ID3] = findgroups(tbl3.method);
        tlisttmp = cell(size(ID3));
        for i = 1:size(ID3,1)
            IDchar = char(ID3(i));
            tlisttmp{i} = regexprep(IDchar,{'pbary','gpr','idw','nn'},{'Barycentric','GPR','IDW','NN'});
        end
        titlelist = categorical(tlisttmp);
        parity3 = splitapply(@(x){x(1)},tbl3.parity,G3);

        if strcmp(datatype,'brk') && (ninputpts == 50000)
            cbnds = [1 500];
        else
            cbnds = [];
        end
        %plotting
        multiparity(parity3,ID3,'hex','cbnds',cbnds,'titlelist',titlelist)

        %extra
        if sgtitleQ
            sgtitle(['ninputpts = ' int2str(ninputpts)])
        end
        
        savefigpng(figfolder,[char(datatype) 'parity' int2str(ninputpts)]);
    end
end

%% errors
methodlist = {'pbary','gpr','idw','nn','avg'};
xytypelbls = ['Barycentric',upper(methodlist(2:4)),'Average'];
% datatypelist = {'brk','kim'};
datatypelist = {'brk'};
ytypes = {'rmse','mae'};
ytypelbls = upper(ytypes);
for datatype = datatypelist
    tbltmp = mdlparstbl(mdlparstbl.datatype==datatype,:);
    multixyplots(tbltmp,methodlist,'ninputpts',ytypes,1,2,'ymin',0,...
        'lgdloc','southwest','ytypelbls',ytypelbls,'xytypelbls',xytypelbls)
    fig = gcf;
    t1 = nexttile(1);
    t1.Legend.Position = [0.299964268608805 0.516990889931764 0.163876007802576 0.229411769754747];
    t2 = nexttile(2);
    t2.Legend.Position = [0.769986057377123 0.592250149191023 0.163876007802576 0.229411769754747];
%     fig.Children.Children(3).Position = [0.780569390710456 0.495824223265097 0.163876007802576 0.229411769754747];
    %saving
    savefigpng(figfolder,[char(datatype) 'error'])
end

%% timing
methodlist = {{'pbary','gpr','idw','nn'},{'pbary','gpr','idw','nn'}};
xytypelbls = {{'Barycentric','GPR','IDW','NN'},{'Barycentric','GPR','IDW','NN'}};
% methodlist = {{'pbary','gpr'},{'gpr','idw','nn'}};
% xytypelbls = {{'Barycentric','GPR'},{'GPR','IDW','NN'}};
tbl3 = mdlparstbl(mdlparstbl.ncores==12,:);
%{ 
whether to multiple pbary by # cores since uses parfor (although rest
might use multiple cores via vectorization..)?
%}
% ids = ismember(tbl3.method,'pbary');
% tbl3(ids,'runtime') = tbl3(ids,'runtime').*tbl3(ids,'cores');
[G3,ID3] = findgroups(tbl3.method);
multixyplots(mdlparstbl,methodlist,'ninputpts',{'runtime'},1,2,'XScale','linear',...
    'YScale','linear','yunits','s','lgdloc','southeast','xytypelbls',xytypelbls)
ax = gca;
% ax.YLim = [0 30];
ax.XScale = 'log';
ax.YScale = 'log';
% ax.YLim(2) = 1000;
ax.Legend.Location = 'southeast';
% stackedplot
%saving
savefigpng(figfolder,'runtime');

%% Runtime Table
% the data is too varied and a log-log plot is too misleading. Going with table instead.
methodlist = {'pbary','gpr','idw','nn'};
titlelist = {'Barycentric','GPR','IDW','NN'};
% tbltmp = tbl3(any(tbl3.method == methodlist,2),:);
[tmu,tsigma] = deal(zeros(8,length(methodlist)));
tstr = cell(8,length(methodlist));

files2 = dir(fullfile('**','interp5DOF-paper'));
folder2 = files2(1).folder;

for i = 1:length(methodlist)
    %unpack
    method = methodlist{i};
    %extract
    tbltmp = tbl3(tbl3.method == method,:);
    [Gtmp,IDtmp] = findgroups(tbltmp.ninputpts);
    %extract
    tmu(:,i) = round(splitapply(@mean,tbltmp.runtime,Gtmp),4);
    tsigma(:,i) = round(splitapply(@std,tbltmp.runtime,Gtmp),4);
    %convert
    tmustr = arrayfun(@(x)num2str(x,'%.4g'),tmu(:,i),'UniformOutput',false);
    tsigmastr = arrayfun(@(x)num2str(x,'%.4g'),tsigma(:,i),'UniformOutput',false);
    %prep
    tstr(:,i) = strcat('$',tmustr,{' '},'\pm',{' '},tsigmastr,'$');
    timeinfo = [titlelist; tstr];
    setsizeinfo = ['\gls{vfzo} Set Size';arrayfun(@(x)['\num{' int2str(x) '}'],IDtmp,'UniformOutput',false)];
    %concatenate
    catinfo = [setsizeinfo,timeinfo];
    %package
    T = table(catinfo);
    disp(T)
    %write
    writetable(T,fullfile(folder2,'runtime.xlsx'),'WriteVariableNames',false)
end

%% NN Mean and Standard Deviation vs. VFZO set size
tbltmp = mdlparstbl(mdlparstbl.datatype == 'brk',:);
[G,ID] = findgroups(tbltmp.ninputpts);
nnmu = splitapply(@mean,tbltmp.nnmu,G);
nnsigma = splitapply(@std,tbltmp.nnmu,G);
paperfigure();
errorbar(ID,nnmu,nnsigma);
set(gca,'XScale','log');
xlabel('VFZO Set Size','Interpreter','latex')
ylabel('VFZO $\omega_{\mathrm{NN}}$ ($^{\circ}$)','Interpreter','latex')
savefigpng(figfolder,'nndist-vs-setsize');

%% dist-parity
seed = 10; %#ok<*UNRCH>
rng(seed);
npts = 10000;
[five,q,nA] = get_five(npts);
o = five2oct(q,nA);
pts = get_octpairs(o);
pts = normr(pts);
pd1 = pdist(pts).';
pd2 = real(pdist(pts,@get_omega).');
pd3 = real(pdist(pts,@get_alen).');
save(fullfile(figfolder,'dist-parity.mat'),'pd1','pd2','pd3')

%% arclength vs. euclidean distance parity
% addpathdir({'get_five.m','cu2qu.m','GBfive2oct.m','gmat2q.m','distance-parity.mat'})
addpathdir({'dist-parity.mat','paperfigure.m'})

fname = 'dist-parity';
load([fname,'.mat'],'pd1','pd2','pd3')

paperfigure()
parityplot(pd1,pd3,'scatter','cscale','linear','xname','Euclidean Distance',...
    'yname','Arc Length','xunits','','yunits','rad','sz',20,'mkr','.','c',[0.6350, 0.0780, 0.1840])

%saving
% savefigpng(folder,fname);
print(fpath,'-dpng','-r300') %saving .fig takes too long, file is pretty big
clear pd1 pd2 pd3

    
%     load('olm_pairwise_oct.mat','oct_new')
%     t = num2cell(oct_new,3);
%     t = cellfun(@(x) squeeze(x).',t.','UniformOutput',false);
%     oAB = vertcat(t{:});
%     oA = get_octpairs(oAB(:,1:8));
%     oB = get_octpairs(oAB(:,9:16));

%% load data for octonion ensemble and fixed test
% addpathdir({'olm_octonion_list.txt','olm_pairwise_distances_cubic.mat'})
A = importdata('olm_octonion_list.txt');
olmoct = A.data;
olmoct = [qinv(olmoct(:,1:4)),qinv(olmoct(:,5:8))];

rng(5)

nptslist = [1 2 10 20];
npts = max(nptslist);
numnpts = length(nptslist);

load('olm_pairwise_distances_cubic.mat','olmpairwisedistancescubic')
pd_olmchesser = table2array(olmpairwisedistancescubic);
pd7 = rad2deg(pd_olmchesser(:));

%% octonion ensemble
[pd4c, pd5c, pd6c, errmetrics] = deal(cell(numnpts,1));
oref = zeros(numnpts,8);
% oref = get_orefs(8);
paperfigure(2,2,14.5);
j = 0;
for i = 1:npts
    if i == 1
        oref(i,:) = get_ocubo(1,'random',[],10);
    else
        oref(i,:) = get_ocubo(); %random reference octonion
%         oref(i,:) = orefs(i,:);
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
        yname = ['$\mathop{}_{\scriptstyle{\forall k \in' yrange '}}^{\rm{min}}2(\frac{180}{\pi})|\hat{o}_{i,k}^{sym}-\hat{o}_{j,k}^{sym}|$'];
        parityplot(pd7,pd4tow,'hex','cscale','log','xname','$\omega$ (traditional)',...
            'yname',yname,'xunits','$^{\circ}$','yunits','$^{\circ}$')
        papertext(j,'xypos',[-0.2,1.0])
    end
end
% Euclidean: '$2\frac{180}{\pi}|o_1^{sym}-o_2^{sym}|$ (this work)'

% klbls = cellfun(@char,num2cell(string(nptslist)),'UniformOutput',false);
% lbls = strcat(klbls,'-');
% lbls{end}(end) = []; %remove last '-'
lbls = lblcat(nptslist);
savefigpng(figfolder,['dist-ensemble-k' lbls])

paperfigure();
enstbl = struct2table(structvertcat(errmetrics{:})); %ensemble table
% G = findgroups(enstbl.ksize);
plot(enstbl.ksize,enstbl.rmse,'-*')
xlabel('$k_{max}$','Interpreter','latex')
ylabel('Error ($^{\circ}$)','Interpreter','latex')
hold on
plot(enstbl.ksize,enstbl.mae,'-*')
legend('RMSE','MAE','Interpreter','latex')
savefigpng(figfolder,['dist-ensemble-rmse-mae'])

%% fixed test - does holding one octonion fixed give same as both varied
tic
pd_fix = get_pd_fix(olmoct);
pdfixtime = toc;
tmp = pd_fix;
pd_fix = rad2deg(pd_fix);
pd_fix = pd_fix(:);

paperfigure();
parityplot(pd7,pd_fix,'scatter','xname','GBO distance','yname','Fixed GBO distance','xunits','deg','yunits','deg')
fixed_errmetrics = get_errmetrics(pd7,pd_fix);

savefigpng(figfolder,'pd-fix')

%% Ensemble Interpolation (Load)
fname = 'ensemble-interp';

seed = 10;
rng(seed);
ninputpts = 50000;
npredpts = 10000;
[~,qm,nA] = get_five(ninputpts);
o = five2oct(qm,nA);
y = GB5DOF_setup(o(:,1:4),o(:,5:8));
[~,qm2,nA2] = get_five(npredpts);
o2 = five2oct(qm2,nA2);
ytrue = GB5DOF_setup(o2(:,1:4),o2(:,5:8));
K = 10;
method = 'gpr';
[ypred,ypredlist,interpfnlist,mdllist,mdlparslist] = ...
    ensembleVFZO([],[],y,[],[],K,method,'ytrue',ytrue,'o',o,'o2',o2);
fpath = fullfile(figfolder,fname);
save(fpath,'ypred','ypredlist','ytrue','qm','nA','qm2','nA2','y','o','o2')
% save(fpath,'ypred','ypredlist','ytrue','qm','nA','qm2','nA2','y','o','o2',...
%     'interpfnlist','mdllist','mdlparslist') %ran into issues with saving, maybe out of memory

%% load Ensemble Interpolation
fname = 'ensemble-interp';
fpath = fullfile(figfolder,fname);
load(fname,'ypred','ypredlist','ytrue','qm','nA','qm2','nA2','y','o','o2')

%% Ensemble Interpolation
paperfigure(2,2,14.509833333333333)
nexttile
parityplot(ytrue,mean([ypredlist{:}],2),'hex')
title('mean')
nexttile
parityplot(ytrue,median([ypredlist{:}],2),'hex')
title('median')
nexttile
parityplot(ytrue,min([ypredlist{:}],[],2),'hex')
title('min')
nexttile
parityplot(ytrue,max([ypredlist{:}],[],2),'hex')
title('max')
savefigpng(figfolder,'ensemble-interp')

%% Ensemble Interpolation
K = length(ypredlist);
% [errmetmean(K),errmetmed(K),errmetmin(K),errmetmax(K)] = deal([]);
clear errmetmean errmetmed errmetmin errmetmax
for i = 1:K
    ypredmean = mean([ypredlist{1:i}],2);
    ypredmed = median([ypredlist{1:i}],2);
    ypredmin = min([ypredlist{1:i}],[],2);
    ypredmax = max([ypredlist{1:i}],[],2);
    errmetmean(i) = get_errmetrics(ypredmean,ytrue); %#ok<SAGROW>
    errmetmed(i) = get_errmetrics(ypredmed,ytrue); %#ok<SAGROW>
    errmetmin(i) = get_errmetrics(ypredmin,ytrue); %#ok<SAGROW>
    errmetmax(i) = get_errmetrics(ypredmax,ytrue); %#ok<SAGROW>
end
rmse.mean = [errmetmean.rmse];
rmse.med = [errmetmed.rmse];
rmse.min = [errmetmin.rmse];
rmse.max = [errmetmax.rmse];

mae.mean = [errmetmean.mae];
mae.med = [errmetmed.mae];
mae.min = [errmetmin.mae];
mae.max = [errmetmax.mae];

rmsecat = [rmse.mean.',rmse.med.',rmse.min.',rmse.max.'];
maecat = [mae.mean.',mae.med.',mae.min.',mae.max.'];
lgdlbl = {'mean','median','min','max'};
paperfigure(1,2);
nexttile
mkr = '*-';
plot(1:K,rmsecat,mkr)
% ax = gca;
% ax.YLim(1) = 0;
xlim([1,K])
legend(lgdlbl,'Location','best')
xlabel('Ensemble Size')
ylabel('RMSE ($J/m^2$)')
nexttile
plot(1:K,maecat,mkr)
% ax = gca;
% ax.YLim(1) = 0;
xlim([1,K])
legend(lgdlbl,'Location','best')
ylabel('MAE ($J/m^2$)')
xlabel('Ensemble Size')

fname = 'ensemble-interp-rmse-mae';
savefigpng(figfolder,fname)

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
% % setpaperdefaults()
% % distance histogram
% seed = 10;
% rng(seed);
% npts = 50000;
% five = get_five(npts);
% o = GBfive2oct(five);
% pts = get_octpairs(o);
% save('oct50000.mat','pts')

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
shadedErrorBar(1:length(mu),mu,sigma)
% errorbar(mu,sigma)
xlabel('k-th NN','Interpreter','latex','FontSize',12)
ylbl = 'average $\omega$ (deg)';
ylabel(ylbl,'Interpreter','latex','FontSize',12)
papertext(2)

% savefigpng(folder,['knnhist',int2str(npts)])
savefigpng(figfolder,['nnhist-knn-',int2str(npts)])

%% barycentric methods
proj_down_test(1)
savefigpng(figfolder,'bary-remove-deg',[130.5 19.5 1886 572])
proj_down_test(2)
savefigpng(figfolder,'bary-delaunay',[131.5 7.5 1872 432])

%% voronoi example
addpathdir('PGnames.mat')
toBPFZ_test(1)
% nexttile
% toBPFZ_test(1)
% nexttile
% toBPFZ_test(4)
fname = 'voronoi';
savefigpng(figfolder,fname)
% savefigpng(folder,'voronoi-4NN')

%% bary interp
sphbary_test()
fname = 'bary-interp';
savefigpng(figfolder,fname,[242.5 351.5 609 393])
% fpath = fullfile(folder,[fname,'.png']);
% I = imread(fpath);
% imshow(I)
% Icropped = imcrop(I,[242.5 351.5 609 393]);
% imwrite(Icropped,fpath)

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
%tunnelplot()

%% Kim GPR Mixture Interpolation, Teaching
rng(10)
[ypred,ysd,ytrue,ci,covmat,kfntmp,kfntmp2,mdl,gprMdl2] = gprmix_test();
savefigpng(figfolder,'kim-interp-teach')

%% Kim GPR Mixture Final Result
rng(10);
% paperfigure(2,2,13.6843);
paperfigure(1,2,7.46983333333333);
nexttile
parityplot(ytrue,ypred);
papertext(1);

n = 300;

nexttile
mdl.method = 'gprmix';
[~,~,~,~,A,B] = tunnelplot(mdl,[],[],n,'brkQ',false,'nnQ',false,'nnQ2',false,'tpredlist',ypred,...
    'kfntmp',kfntmp,'kfntmp2',kfntmp2,'gprMdl2',gprMdl2);
legend('GPR mixture model','Location','best','Interpreter','latex')
papertext(2);
disp(A)
disp(B)

% nsamp = 10;
% gprsamples = mvnrnd(ypred,covmat,nsamp);


savefigpng(figfolder,'kim-interp')

%% Kim Tunnel Plots and Posterior Sampling
paperfigure(2,2,13.6843);
rng(10)
i = 0;
nexttile
i = i+1;
npts = size(ytrue,1);
rids = randperm(npts,npts);
parityplot(ytrue(rids),ypred(rids),'scatter','c',ysd(rids),'cblbl','$\sigma_{\mathrm{pred}} (J m^{-2})$')
papertext(i);

nexttile
i = i+1;
mdl.method = 'gprmix';
tunnelplot(mdl,A,B,n,'brkQ',false,'nnQ',true,'nnQ2',true,'tpredlist',ypred,...
    'kfntmp',kfntmp,'kfntmp2',kfntmp2,'gprMdl2',gprMdl2);
ax = nexttile(i);
ch = ax.Children([1,2]);
legend(ch,'GPR mixture','1st (input) NN','Location','southwest','Interpreter','latex')
papertext(i);

nexttile
i = i+1;
mdl.method = 'gprmix';
tunnelplot(mdl,A,B,n,'brkQ',false,'nnQ',false,'nnQ2',false,'tpredlist',ypred,...
    'kfntmp',kfntmp,'kfntmp2',kfntmp2,'gprMdl2',gprMdl2,'nsamp',5);
papertext(i);

nexttile
i = i+1;
mdl.method = 'gprmix';
tunnelplot(mdl,[],[],n,'brkQ',true,'nnQ',false,'nnQ2',false,'tpredlist',ypred,...
    'kfntmp',kfntmp,'kfntmp2',kfntmp2,'gprMdl2',gprMdl2);
legend('Ni BRK','Fe GPR mixture','Location','southwest','Interpreter','latex')
papertext(i);
savefigpng(figfolder,'kim-interp-posterior')

%% GPR Mixture Sigmoid Function
paperfigure();
thr = 1.1;
scl = 30;
sigfn = @(x,scl,xshift) 1./(1+exp(-scl*(x-xshift)));
x = 0.5:0.01:1.5;
plot(x,sigfn(x,scl,thr),'LineWidth',1.5);
ylim([0 1])
xlabel('GBE ($J m^{-2}$)','Interpreter','latex')
ylabel('Mixing fraction (f)','Interpreter','latex')
axis square
savefigpng(figfolder,'gprmix-sigmoid')

%% tunnel plot (full)
ninputpts = 50000; %1000, 50000
n = 150;
rng(10)
paperfigure();
[tpredlist,tsdlist,propList,methodlist,A,B] = tunnelplot_test(2,ninputpts,n);
fname = ['tunnel-',int2str(ninputpts)];
fpath = fullfile(figfolder,fname);
save(fpath,'tpredlist','tsdlist','propList','methodlist','A','B')
disp(A)
disp(B)
savefigpng(figfolder,fname)

%% tunnel plot (shortcut)
ninputpts = 1000;
n = 150;
fname = ['tunnel-' int2str(ninputpts)];
fpath = fullfile(figfolder,fname);
load(fpath,'tpredlist','tsdlist','propList','methodlist','A','B')
paperfigure();
tunnelplot_test(2,ninputpts,n,tpredlist,tsdlist,propList,methodlist,A,B);
disp(A)
disp(B)
savefigpng(figfolder,fname)


%% Kim Interp Degeneracy Info
fname = 'kim-interp-degeneracy'; % produced in Kim2oct.m
fpath = fullfile(figfolder,fname);
load(fpath,'errmetrics')
paperfigure(2,2);
i = 0;
nexttile
i = i+1;
e = errmetrics.e;
ecat = vertcat(e{:});
nprops = errmetrics.nprops;

ids = nprops <= 4;
e1 = e(ids);
ecat1 = vertcat(e1{:});

nprops = errmetrics.nprops;
histogram(ecat1);
axis square
% set(gca,'YScale','log');
xlabel({'Error relative to (1-4 degeneracy)','set average ($J/m^2$)'},'Interpreter','latex')
ylabel('counts','Interpreter','latex')

min(ecat)
max(ecat)
papertext(i);

nexttile
i=i+1;
propavg = vertcat(errmetrics.propavg{ids});
props = vertcat(errmetrics.propsets{ids});
parityplot(propavg,props)
xlabel('Average set-wise GBE ($J/m^2$)','Interpreter','latex')
ylabel('Input GBE ($J/m^2$)','Interpreter','latex')

ids2 = nprops > 4;
e2 = e(ids2);
ecat2 = vertcat(e2{:});
papertext(i);

nexttile
i = i+1;
histogram(ecat2)
axis square
xlabel({'Error relative to (5+ degeneracy)','set average ($J/m^2$)'},'Interpreter','latex')
ylabel('counts','Interpreter','latex')
papertext(i);

nexttile
i = i+1;
propavg = vertcat(errmetrics.propavg{ids2});
props = vertcat(errmetrics.propsets{ids2});
parityplot(propavg,props)
xlabel('Average set-wise GBE ($J/m^2$)','Interpreter','latex')
ylabel('Input GBE ($J/m^2$)','Interpreter','latex')
papertext(i);

fname = 'kim-interp-degeneracy-results';
savefigpng(figfolder,fname)

paperfigure()
histogram(errmetrics.nprops)
axis square
xlabel('Number of degenerate GBs per set','Interpreter','latex')
ylabel('Number of sets','Interpreter','latex')
fname = 'kim-interp-degeneracy-sets';
savefigpng(figfolder,fname)
% 
% nexttile
% histogram(nprops(ids))

disp(['wtrmse: ' num2str(errmetrics.wtrmse)])
disp(['wtmae: ' num2str(errmetrics.wtmae)])

%% Olmsted Ni Dataset Comparison
%see also interp5DOF_setup.m
A = importdata('olm_octonion_list.txt');

o = A.data;
o = [qinv(o(:,1:4)),qinv(o(:,5:8))];
o = get_octpairs(o);

B = importdata('olm_properties.txt');
y = B.data(:,1);

[qm,nA] = oct2five(o);
qm2 = qm;
nA2 = nA;
ytrue = y;

[ypred,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,y,qm2,nA2,'ytrue',ytrue);
[ypred_lowSig,interpfn_lowSig,mdl_lowSig,mdlpars_lowSig] = interp5DOF(qm,nA,y,qm2,nA2,'ytrue',ytrue,...
    'mygpropts',{'PredictMethod','exact','Sigma',1e-4,'ConstantSigma',true,'SigmaLowerBound',5e-5});

y_brk = GB5DOF_setup(o(:,1:4),o(:,5:8));
gpr_errmetrics = get_errmetrics(ytrue,ypred);
gpr_lowSig_errmetrics = get_errmetrics(ytrue,ypred_lowSig);
brk_errmetrics = get_errmetrics(ytrue,y_brk);

parity{1} = struct('ytrue',ytrue,'ypred',y_brk);
parity{2} = struct('ytrue',ytrue,'ypred',ypred);
parity{3} = struct('ytrue',ytrue,'ypred',ypred_lowSig);

strs={'BRK','GPR','low-noise-GPR'};
for i = 1:3
    paperfigure();
    str = strs{i};
    parityplot(parity{i}.ypred,parity{i}.ytrue,'yname',['predicted ' str]);
    str = lower(str);
    savefigpng(figfolder,['resubloss-ni-' str])
end

S = structvertcat(mdlpars,mdlpars_lowSig);
KernelParameters = vertcat(S.KernelParameters);
SigmaL = 2*rad2deg(KernelParameters(:,1));
SigmaF = KernelParameters(:,2);
KernelParameterNames = S.KernelParameterNames;
Beta = vertcat(S.Beta);
Sigma = vertcat(S.Sigma);
ConstantSigma = {'no';'yes'};
T = table(ConstantSigma,SigmaL,SigmaF,Beta,Sigma,'VariableNames',...
    {'Constant $\sigma$','$\sigma_L$ ($^\circ{}$)',...
    '$\sigma_F$ ($J m^{-2}$)',...
    '$\beta$ ($J m^{-2}$)',...
    '$\sigma_\mathrm{in}$ ($J m^{-2}$)'});
caption = ['Fitted parameters for two \gls{gpr} models fitted to the 388 '...
    'simulated Ni \glspl{gbe} by \citet{olmstedSurveyComputedGrain2009a}. '...
    'The first model allows $\sigma$ to vary, whereas the second constrains '...
    '$\sigma$ to be fixed. $\sigma_L$, $\sigma_F$, $\beta$, and $\sigma$ are '...
    'the kernel length scale, signal standard deviation, constant basis function, '...
    'and input property standard deviation, respectively. See '...
    '\url{https://www.mathworks.com/help/stats/fitrgp.html} for additional details.'];
savetblstr(T,'resubloss-ni-pars',tblfolder,caption)

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
[G2,ID2] = findgroups(mdlparstbl.ninputpts);
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


    % tbl1 = tbl(:,{'method','ninputpts','ndatapts'});


    tbl2 = mdlparstbl(:,{'method','ndatapts'});
    [G2,ID2] = findgroups(ID1.method);
    ninputpts1 = splitapply(@(x){x},mdlparstbl.ninputpts,G2);


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
    [G1,ID1] = findgroups(tbl.ninputpts);
    
    [npts,rmse,mae,runtime] = deal(tbl.ninputpts,tbl.rmse,tbl.mae,tbl.runtime);
    
    ninputpts = splitapply(@(x)x(1),npts,G1);
    medRMSE = splitapply(@median,rmse,G1);
    stdRMSE = splitapply(@std,rmse,G1);
    medMAE = splitapply(@median,mae,G1);
    stdMAE = splitapply(@std,mae,G1);
    medruntime = splitapply(@median,runtime,G1);
    stdruntime = splitapply(@std,runtime,G1);
    
    semilogx(t1,ninputpts,medRMSE,'-o')
    xlabel('ninputpts')
    ylabel('RMSE (J/m^2)')
    axis square
    hold off
    
    semilogx(t2,ninputpts,medMAE,'-o')
    xlabel('ninputpts')
    ylabel('MAE (J/m^2)')
    axis square
    hold off
%     semilogx(ninputpts,medruntime,'-o')

    figure(fig2)
    hold on
    semilogx(ninputpts,medruntime,'-o')
    xlabel('ninputpts')
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


        tlisttmp = cellfun(@(ID) regexprep(ID,{'pbary','gpr','idw','nn'},{'Barycentric','GPR','IDW','NN'}),cellstr(ID3),'UniformOutput',false);

% addpathdir({'paperfigure.m','dist-parity.mat','olm_octonion_list.txt',...
%     'oct50000.mat','gmat2q.m','PGnames.mat','var_names.m'})

% multiparity(parity,{'BRK','VFZ','VFZ-lowSig'})
    %     ax=nexttile(i);
%     str=strs{i};
%     ax.XLabel.String = ['actual simulated GBE ($J m^{-2}$)'];
%     ax.YLabel.String = ['predicted ' str ' GBE ($J m^{-2}$)'];

%}