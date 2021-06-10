%% setup
%ensure that datafolder path contains interp\data
files = dir(fullfile('**','data'));
folders = {files.folder};
ids = find(cellfun(@(x) contains(x,fullfile('interp','data')),folders));
if isempty(ids)
    warning('interp/data directory not found')
end
id = ids(1);
datafolder = files(id).folder;

%ensure that figfolder path contains interp\results
files = dir(fullfile('**','results'));
folders = {files.folder};
ids = find(cellfun(@(x) contains(x,fullfile('interp','results')),folders));
if isempty(ids)
    warning('interp/results directory not found')
end
id = ids(1);
figfolder = files(id).folder;

fnames = {'gpr58604_gitID-3b7183c_puuID-ce759533_kim-Fe-train-all-data-fic.mat',...
    'gpr388_gitID-63ce950_puuID-7e693646_olmsted-Ni-rng11.mat'}; %'gpr388_gitID-57857cd_puuID-e9c787cd_olmsted-Ni'

urls = {'https://ndownloader.figshare.com/files/28256244?private_link=0ac714343cac8b2a5935',...
    'https://ndownloader.figshare.com/files/28256280?private_link=0ac714343cac8b2a5935'};

setlatex()

%% load (or download & load) data
for i = 1:length(fnames)
    fname = fnames{i};
    url = urls{i};
    datapath = fullfile(datafolder,fname);
    if exist(datapath,'file') ~= 2
        warning(['Could not find file at ' datapath])
        % Prompt to continue
        m = input('Download and continue? y/n:','s');
        if ~strcmp(m,'y') && ~strcmp(m,'Y')
            return
        else
            websave(datapath,url)
        end
    end
    load(datapath)
end

%% NBO pairs
nnbo = 2500;
nbo = repmat(get_cubo(nnbo),1,2);
nbo = normr(get_octpairs(nbo));
[idx,d] = knnsearch(nbo,nbo,'K',2,'IncludeTies',true);
d = vertcat(d{:});
d = 2*rad2deg(d(:,2)); %convert to "Toby" distance (i.e. traditional octonion)
paperfigure();
histogram(d) %histogram of pairwise distances
xlabel('$\omega_{\mathrm{NN}} (^\circ{})$','Interpreter','latex')
ylabel('Number of NBOs','Interpreter','latex')
savefigpng(figfolder,['nbo' int2str(nnbo) '-nnhist'])

%% posterior data setup
savenames = {'gpr58604-kim','gpr388-olmsted'};
nfnames = length(fnames);

load('oct50000.mat','pts')
npts = 15000; %can produce large matrices, at least (npts x npts)
pts = pts(1:npts,:);

A = importdata('olm_octonion_list.txt');
oolm = A.data;
oolm = oflip(oolm);
oolm = get_octpairs(oolm);

B = importdata('olm_properties.txt');
yolm = B.data(:,1);

load('Kim2011_Fe_oct_GBE','meshList','propList','specIDs','mechIDs');

okim = meshList(specIDs,:);
ykim = propList(specIDs);

olm_pts = normr([oolm;nbo;pts]);
kim_pts = normr([okim;nbo;pts]);

olm_ids = [3 169 32 21 33];
kim_ids = [7 162 259 315 406];

idlist = {kim_ids,olm_ids};
ptslist = {kim_pts,olm_pts};
ylist = {ykim,yolm};
Sigma = [3 5 7 9 11];

%% posterior data
%initialize
[covmats,ds,ypredlist,ysdlist,fivelist] = deal(cell(1,nfnames));
for i = 1:nfnames
    %unpack
    fname = fnames{i};
    ids = idlist{i};
    pts = ptslist{i};
    y = ylist{i};
    S = load(fname,'mdl');
    mdl = S.mdl;
    gprMdl = mdl.cgprMdl;
    
    %conversion to 5DOF
    [qm,nA] = oct2five(pts);
    five.qm = qm;
    five.nA = nA;
    
    %svd transformation
    ppts = proj_down(pts,mdl.projtol,mdl.usv,'zero',mdl.zeroQ);
    
    %covariance kernel/matrices
    kfcn = gprMdl.Impl.Kernel.makeKernelAsFunctionOfXNXM(gprMdl.Impl.ThetaHat);
    covmat = kfcn(ppts,ppts);
    
    % pairwise distances
    d = pdist(pts);
    
    %predictions
    [ypred,ysd] = predict(gprMdl,ppts);
    
    %saving
    savename = savenames{i};
    savepath = fullfile(datafolder,savename);
    
    savepath2 = [savepath '-cov'];
    warning(['saving ' savepath2 '.mat'])
    save(savepath2,'covmat','-v7.3')
    
    savepath3 = [savepath '-dist'];    
    warning(['saving ' savepath3 '.mat'])
    save(savepath3,'d')
    
    savepath4 = [savepath '-other'];
    warning(['saving ' savepath4 '.mat'])
    save(savepath4,'ids','Sigma','pts','five','y','ypred','ysd')
    
    %package
    %--comment these lines if you run out of memory
    fivelist{i} = five;
    covmats{i} = covmat;
    ds{i} = d;
    ypredlist{i} = ypred;
    ysdlist{i} = ysd;
    clear five covmat d ypred ysd
end

%% CODE GRAVEYARD
%{
%     writematrix(covmats{i},[savepath '-cov.csv']);

% savename = ['eucl' int2str(npts+nnbo)];
% % writematrix(d,fullfile(datafolder,[savename '-dist.csv']));
% % writematrix(pts,fullfile(datafolder,[savename '-dist.csv']));
% save(fullfile(datafolder,savename),'d','pts')

%     files = dir(fullfile('**',fname));
%     fpath = fullfile(files(1).folder,files(1).name);
%     load(fpath)

% qinv(oolm(:,1:4)),qinv(oolm(:,5:8))];
%}