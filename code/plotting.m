%split apply & find groups
% fname = 'gitID-76fca8c_uuID-f51500cd_paper-data.mat';
% fname = 'gitID-396aaa2_uuID-6816f860_paper-data2.mat';
fname = 'gitID-f585733_uuID-edf2fcc7_paper-data2.mat';
files = dir(fullfile('**',fname));
fpath = fullfile(files(1).folder,files(1).name);
load(fpath);

slurmQ = 0;

%saving
files = dir(fullfile('**','interp5DOF-paper','figures'));
folder = files(1).folder;

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
        
        savefigpng(folder,[datatype 'parity' int2str(nmeshpts)]);
    end
end

%% errors
methodlist = {'pbary','gpr','idw','nn','avg'};

multixyplots(mdlparstbl,methodlist,'nmeshpts',{'rmse','mae'},2,1,'ymin',0,'lgdloc','southwest')
%saving
savefigpng(folder,[datatype 'error'])

%% timing
methodlist = {{'pbary','gpr'},{'idw','nn','avg'}};
tbl3 = mdlparstbl(mdlparstbl.ncores==24,:);
%{ 
whether to multiple pbary by # cores since uses parfor (although rest
might use multiple cores via vectorization..)?
%}
% ids = ismember(tbl3.method,'pbary');
% tbl3(ids,'runtime') = tbl3(ids,'runtime').*tbl3(ids,'cores');
[G3,ID3] = findgroups(tbl3.method);
multixyplots(mdlparstbl,methodlist,'nmeshpts',{'runtime'},2,1,'YScale','linear','yunits','s','lgdloc','north')
%saving
savefigpng(folder,'runtime');

%% arclength vs. euclidean distance parity
if slurmQ
    seed = 10; %#ok<*UNRCH>
    rng(seed);
    npts = 10000;
    five = get_five(npts);
    o = GBfive2oct(five);
    pts = get_octpairs(o);
    pts = normr(pts);
    pd1 = pdist(pts).';
    pd2 = real(pdist(pts,@get_omega).');
    pd3 = real(pdist(pts,@get_alen).');
    figure
    parityplot(pd1,pd3,'scatter','cscale','linear','xname','Euclidean Distance','yname','Arc Length','xunits','','yunits','rad')
    %saving
    fname = 'distance-parity';
    savefigpng(folder,fname);
end
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

%% distance histogram
addpathdir('gmat2q.m')
seed = 10;
rng(seed);
npts = 50000;
five = get_five(npts);
o = GBfive2oct(five);
pts = get_octpairs(o);
nnhist(pts,'omega');
savefigpng(folder,['disthist',int2str(npts)])

%% knn distances
[nnpts,D,mu,sigma] = get_knn(pts,'omega',100);
fig = figure;
fig.Position = [460.2 245 498.4 472.8];
errorbar(mu,sigma)
xlabel('k-th NN')
ylbl = 'average \omega (deg)';
ylabel(ylbl)
savefigpng(folder,['knnhist',int2str(npts)])

%% barycentric methods
proj_down_test(1)
savefigpng(folder,'baryRemoveDeg')
proj_down_test(2)
savefigpng(folder,'baryDelaunay')

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


%}