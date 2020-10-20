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

%split apply & find groups
fname = 'gitID-76fca8c_uuID-f51500cd_paper-data.mat';
files = dir(fullfile('**',fname));
fpath = fullfile(files(1).folder,files(1).name);
load(fpath);
tbl1 = mdlparstbl(:,{'method','nmeshpts','ndatapts'});
[G1,ID1] = findgroups(tbl1);
parity = splitapply(@(x){x(1)},mdlparstbl.parity,G1);
nmeshpts = splitapply(@(x)x(1),mdlparstbl.nmeshpts,G1);
medRMSE = splitapply(@(x)median(x),mdlparstbl.rmse,G1);
medMAE = splitapply(@(x)median(x),mdlparstbl.mae,G1);
medruntime = splitapply(@(x)median(x),mdlparstbl.runtime,G1);

tbl2 = mdlparstbl(:,{'method','ndatapts'});
[G2,ID2] = findgroups(ID1.method);
nmeshpts1 = splitapply(@(x){x},mdlparstbl.nmeshpts,G2);


%% parity plot
methodlist = {'pbary','gpr','idw','nn'};
sgtitleQ = false;
for nmeshpts = [388 10000 50000]
    %extract parity and IDs
    tbl3 = mdlparstbl(ismember(mdlparstbl.method,methodlist) & mdlparstbl.nmeshpts==nmeshpts,:);
    [G3,ID3] = findgroups(tbl3.method);
    parity3 = splitapply(@(x){x(1)},tbl3.parity,G3);
    
    %plotting
    multiparity(parity3,ID3)
    
    %extra
    if sgtitleQ
        sgtitle(['nmeshpts = ' int2str(nmeshpts)])
    end
    
    %saving
    files = dir(fullfile('**','interp5DOF-paper','figures'));
    folder = files(1).folder;
    fname = ['brkparity',int2str(nmeshpts)];
    fpath = fullfile(folder,fname);
    savefig(fpath)
    print(fpath,'-dpng')
end
    

%% distance parity
pd1 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts).';
pd2 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_omega).';
pd3 = pdist(mdltbl(strcmp(mdltbl.uuid,uuid),:).mesh{1}.pts,@get_alen).';
parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','')

%% distance parity
pts = normr(rand(388,3)-0.5);
pts2 = normr(rand(388,3)-0.5);
pd1 = pdist(pts);
pd3 = pdist(pts,@get_alen);
parityplot(pd1,pd3,'xname','euclidean','yname','arclength','units','','title','388 3D points')
figure
t=n2c(pts);
scatter3(t{:})
axis equal


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

%}