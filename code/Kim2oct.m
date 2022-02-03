% KIM2OCT  load Kim2011 dataset, convert to norm-symmetrized octonions
%  and write to file with corresponding GB energies
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Dependencies:
%  addpathdir.m
%	'Kim2011_FeGBEnergy.txt'
%  GBfive2oct.m
%  get_octpairs.m
%
% Notes:
%  *
%--------------------------------------------------------------------------
clear; close all
%% setup

epsijk = 1;
avgQ = true;
removezeroQ = true;

addpathdir({'eu2qu.m','q2rod.m','get_octpairs.m','GBfive2oct.m','Kim'})

files = dir(fullfile('**','interp','code','Kim'));
kimfolder = files(1).folder;
% if exist(folder,'dir') ~= 7
% 	mkdir(folder)
% end

%load mechanically selected Kim data
% fname = 'Kim2011_FeGBEnergy.txt';
fname = 'Fe_BCC_Nor_DB.txt';
txt = fileread(fname);

%convert 'D' (meaning double precision, base 10) to 'e', as in 1e2 == 100
txt = strrep(txt,'D','e');
fname2 = [fname(1:end-4) '_matlab.txt'];
fpath = fullfile(kimfolder,fname2);
fid=fopen(fpath,'w');
fprintf(fid,txt);
fclose(fid);

%read in data from file
meshTable = readtable(fname2,'HeaderLines',15,'ReadVariableNames',true);
datatmp = table2array(meshTable);

%load intentionally selected Kim data
% fnamesym = 'Kim2011_FeGBEnergy_SymSubset.txt';
fnamesym = 'Fe_BCC_Spe_DB.txt';
txtsym = fileread(fnamesym);
txtsym = strrep(txtsym,'D','e');
fnamesym2 = [fname(1:end-4) '_matlab.txt'];
fpathsym = fullfile(kimfolder,fnamesym2);
fid=fopen(fpathsym,'w');
fprintf(fid,txtsym);
fclose(fid);

%read in data from file
meshTableSym = readtable(fnamesym2,'HeaderLines',15,'ReadVariableNames',true);
datatmpsym = table2array(meshTableSym);

%number of mechanically seleected points
npts = size(datatmp,1);
disp(['# mechanically selected pts: ' int2str(npts)])

nptssym = size(datatmpsym,1);
disp(['# intentionally selected points: ' int2str(nptssym)])

nptstot = npts+nptssym;
disp(['# total points: ' int2str(nptstot)])

%concatenate mechanical and intentional points
data5dof = [datatmp(:,1:end-1);datatmpsym(:,2:end-1)]; %ignore column that contains "Sigma" in datatmpsym
gbe = [datatmp(:,end);datatmpsym(:,end)];

%% pick out low GBE/low Sigma GBs
cl = readcell(fnamesym2);
Sigma = cl(17:end,1);

numIDs = cellfun(@isnumeric,Sigma);
for i = 1:length(Sigma)
    if isnumeric(Sigma{i})
        Sigma{i} = int2str(Sigma{i});
    end
end

% Sigma = datatmpsym(:,1);
id = unique(Sigma);
SigmaList = {'3' '5' '7' '9' '11'};
firstIDs = cellfun(@(S) find(strcmp(id,S)),SigmaList);
id = [id(firstIDs);id(setdiff(1:length(id),firstIDs))];
nid = length(id);
minIDs = cell(1,nid);
ykim = datatmpsym(:,end);
for i = 1:nid
    S = id(i);
    gbetmp = ykim; %reset
    gbetmp(~strcmp(Sigma,S)) = inf; %(Sigma ~= S)
    [~,minIDs{i}] = min(gbetmp);
end
Slist = id(1:5);
kim_ids = [minIDs{1:5}];

%SigmaNum = cellfun(@str2num, id, 'UniformOutput', false);
nSigma = length(Sigma);
SigmaNum = zeros(1,nSigma);
for i = 1:nSigma
    str = Sigma{i};
    SigmaNum(i) = str2double(str);
    if isnan(SigmaNum(i))
        str = str(1:end-1);
        SigmaNum(i) = str2double(str);
    end
end
paperfigure();
hexscatter(SigmaNum.', ykim, minmax(SigmaNum), minmax(ykim.')+1e-3*[-1,1], 'reflineQ', false, 'axis', 'square', 'cscale', 'linear')
xlabel('$\Sigma$','Interpreter','latex')
ylabel('GBE $(J m^{-2})$','Interpreter','latex')
hold on
ax = plot(SigmaNum([minIDs{:}]), ykim([minIDs{:}]),'k');
legend(ax,'Lower Bound','Interpreter','latex','Location','south')
savefigpng('C:\Users\sterg\Pictures\qual-exam', 'sigma-kim')

%% conversion to octonions
%extract 5DOF parameters
data5dof = deg2rad(data5dof);
t=n2c(data5dof);
[phi1,Phi,phi2,po,az] = t{:};
eulist = [phi1 Phi phi2]; %catenate euler angles

%convert to quaternions & cartesian normal pairs
% +1 or -1?, based on paper, misorientation seems to be defined in the active sense, but epsijk==-1 gives better results
initialepsijk = 1; %epsijk==1 also seems to give more expected results with GBdist4 comparisons, see "Extra" at bottom
qlist = eu2qu(eulist,initialepsijk);
el = po2el(po); %convert polar angle to elevation angle
[x,y,z] = sph2cart(az,el,ones(nptstot,1));
nAlist = [x y z];

%get octonion mesh
meshListTmp = five2oct(qlist,nAlist,epsijk);
meshListFull = get_octpairs(meshListTmp,epsijk);

mechIDs = [true(1,npts),false(1,nptssym)];
specIDs = ~mechIDs;

%% get property list
propListFull = gbe/1000; % convert from mJ/m^2 to J/m^2

%% save low-GBE/low-Sigma GBs
kimABCDE = num2cell(meshListTmp(kim_ids+npts,:),2);
ykimABCDE = gbe(kim_ids+npts);
files = dir(fullfile('**','interp','data'));
datafolder = files(1).folder;
save(fullfile(datafolder,'kim-oct-Sigma-3-5-7-9-11'),'kimABCDE','ykimABCDE','kim_ids','Slist')
lowSig_IDs = zeros(1,nptstot);
lowSig_IDs(kim_ids) = 1;

if removezeroQ
    %% remove GBs with GBE near 0
    ids = find(propListFull > 0.01);
    meshList = meshListFull(ids,:);
    propList = propListFull(ids);
    
    mechIDs = mechIDs(ids);
    specIDs = specIDs(ids);
    kimIDs = kimIDs(ids);
    
    qlist = qlist(ids,:);
    nAlist = nAlist(ids,:);
end

%% average properties for repeat octonions and remove repeats (except one)
if avgQ
    [meshList,propList,rmIDlist,keepIDs,rmIDcell,errmetrics] = avgrepeats(meshList,propList,'mean'); %#ok<*UNRCH>
    files = dir(fullfile('**','interp5DOF-paper','figures'));
    if ~isempty(files)
        figfolder = files(1).folder;
        save(fullfile(figfolder,'kim-interp-degeneracy.mat'),'errmetrics')
    end
    mechIDs(rmIDlist) = [];
    specIDs(rmIDlist) = [];
    kimIDs(rmIDlist) = [];
    
    qlist(rmIDlist,:) = [];
    nAlist(rmIDlist,:) = [];
    
else
    propList = propListFull;
    meshList = meshListFull;
end

%number of points after averaging
npts2 = size(meshList,1);
% disp(['# pts (after repeat and GBE=0 removal): ' int2str(npts2)])
disp(['# pts (removezeroQ==' int2str(removezeroQ) ', avgQ==' int2str(avgQ) '): ' int2str(npts2)])

nptsmech = nnz(mechIDs);
nptsspec = nnz(specIDs);

disp(['# mechanically selected pts (removezeroQ==' int2str(removezeroQ) ', avgQ==' int2str(avgQ) '): ' int2str(nptsmech)])
disp(['# special pts (removezeroQ==' int2str(removezeroQ) ', avgQ==' int2str(avgQ) '): ' int2str(nptsspec)])

%package q & nA pairs
five = struct('q',qlist,'nA',nAlist);

%% write files
%write octonions and GB Energy to .txt file
fname = 'Kim2011_Fe_oct_GBE.txt';
fpath = fullfile(kimfolder,fname);
fid = fopen(fpath,'w');
fprintf(fid,[...
	'#------------------------------------------------- \n' ...
	'#Calculated grain boundary energies of bcc Fe for \n' ...
    '#mechanically selected 66,339 grain boundaries and \n' ...
    '#intentionally selected 2,366 special boundaries \n' ...
	'#Columns 1:8 : Euclidean grain boundary octonion coordinates \n' ...
	'#Column 9 : Grain boundary energy (J/m^2) \n' ...
	'#------------------------------------------------- \n']); %6 header lines
% writematrix([meshList propList],fname,'WriteMode','append','Delimiter','tab') %only works in 2020a
ftmp = 'temp.txt';
ftmppath = fullfile(kimfolder,ftmp);
writematrix([meshList propList],ftmppath,'Delimiter','tab');
txtsym = fileread(ftmppath);
fprintf(fid,txtsym);
fclose(fid);
%save to .mat file
save(fpath(1:end-4),'meshList','propList','five','mechIDs','specIDs','kimIDs','meshTable')


%% Extra Commentary
%{
% initialepsijk = 1;
% 
% 38000	38001	0.1864
% 38000	65000	0.2298
% 38000	38003	0.2567
% 
% initialepsijk = -1;
% 
% 38000	38001	0.1864
% 38000	65000	0.3366
% 38000	38003	0.3378
% 
% Since norm(nAlist(ids(38000),:)-nAlist(ids(38003),:)) == 0.5002 = ~0.2567*2, this makes me think that epsijk == 1 is the correct interpretation. This is also supported by the description in the paper of the misorientation convention being qinv(qA)**qB
%}

%----------------------------------CODE GRAVEYARD--------------------------
%{

% varNames = meshTable.Properties.VariableNames; %Euler angles (misorientation, deg), the polar & azimuth (inclination, deg), GBE (mJ/m^2)

%[meshList,ia] = uniquetol(meshListFull,'ByRows',true);
%disp(['# unique pts: ' int2str(npts)])

%setGlobal_epsijk(-1);  

%','setGlobal_epsijk.m'


%get octonion mesh
meshListTmp = GBfive2oct(qlist,nAlist);
meshListFull = get_octpairs(meshListTmp);

%Kim references Bunge notation, so using that convention for rotation conversion


datatmp2 = deg2rad(datatmp(:,1:end-1));

% %convert from active to passive convention?
% qlist = qinv(qlist);


%split back into mechanically and intentionally selected GBs
% qmech = qlist(1:npts,:);
% nAmech = nAlist(1:npts,:);
% 
% qspec = qlist(npts+1:end,:);
% nAspec = nAlist(npts+1:end,:);

% meshmech = meshListFull(1:npts,:);
% meshspec = meshListFull(npts+1:end,:);


% propmech = propListFull(1:npts,:);
% propspec = propListFull(npts+1:end,:);


%     [meshmech,propmech] = avgrepeats(meshmech,propmech); %#ok<*UNRCH>
%     [meshspec,propspec] = avgrepeats(meshspec,propspec); %#ok<*UNRCH>



%     idsmech = find(propmech);
%     meshmech = meshmech(ids,:);
%     propmech = propmech(ids,:);
%     
%     idsspec = find(propspec);
%     meshspec = meshspec(ids,:);
%     propspec = propspec(ids,:);
    %     qmech = qmech(idsmech,:);
%     nAmech = nAmech(idsmech,:);
%     
%     qspec = qspec(idsspec,:);
%     nAspec = nAspec(idsspec,:);



        % mechIDs = 1:npts;
% specIDs = npts+1:nptstot;
%     mechIDs = setdiff(mechIDs,rmIDlist);
    
%     specIDs = setdiff(specIDs,rmIDlist);
%     nmechIDs = numel(mechIDs);
%     nspecIDs = numel(specIDs);
    
%     mechIDs = 1:nmechIDs;
%     specIDs = (nmechIDs+1):(nmechIDs+nspecIDs);

%     mechIDs = intersect(mechIDs,ids); %note that sorting occurs, but shouldn't matter because mechIDs is already sorted
%     specIDs = intersect(specIDs,ids);
    
%     nmechIDs = numel(mechIDs);
%     nspecIDs = numel(specIDs);
    
%     mechIDs = 1:nmechIDs;
%     specIDs = (nmechIDs+1):(nmechIDs+nspecIDs);
%}