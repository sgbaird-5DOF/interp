%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-07-27
%
% Description: load Kim2011 dataset, convert to norm-symmetrized octonions,
% and write to file with corresponding GB energies
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

addpathdir({'eu2qu.m','q2rod.m','get_octpairs.m','GBfive2oct.m'})

folder = 'Kim';
if exist(folder,'dir') ~= 7
	mkdir(folder)
end

%load Kim data
fname = 'Kim2011_FeGBEnergy.txt';
txt = fileread(fname);

%convert 'D' (meaning double precision, base 10) to 'e', as in 1e2 == 100
txt = strrep(txt,'D','e');
fname2 = [fname(1:end-4) '_matlab.txt'];
fpath2 = fullfile(folder,fname2);
fid=fopen(fpath2,'w');
fprintf(fid,txt);
fclose(fid);

%read in data from file
meshTable = readtable(fname2,'HeaderLines',9,'ReadVariableNames',true);
datatemp = table2array(meshTable);

%number of points
npts = size(datatemp,1);
disp(['# pts: ' int2str(npts)])

%% conversion to octonions
%extract 5DOF parameters
datatemp2 = datatemp(:,1:end-1);
t=n2c(datatemp2);
[phi1,Phi,phi2,pol,az] = t{:};
eulist = [phi1 Phi phi2]; %catenate euler angles

%convert to quaternions & cartesian normal pairs
%Kim references Bunge notation, so using that convention for rotation conversion
qlist = eu2qu(eulist,-1);
el = 90-pol; %convert polar angle to elevation angle
[x,y,z] = sph2cart(az,el,ones(npts,1));
nAlist = [x y z];

%package q & nA pairs
five = struct('q',qlist,'nA',nAlist);

%get octonion mesh
meshListFull = GBfive2oct(qlist,nAlist);

%% get property list
propListFull = datatemp(:,end)/1000; % convert from mJ/m^2 to J/m^2

%% average properties for repeat octonions, remove repeats (except one)
[meshList,propList] = avgrepeats(meshListFull,propListFull);

%% write files
%write octonions and GB Energy to .txt file
fname = 'Kim2011_Fe_oct_GBE.txt';
fpath = fullfile(folder,fname);
fid = fopen(fpath,'w');
fprintf(fid,[...
	'#------------------------------------------------- \n' ...
	'#Calculated grain boundary energies of bcc Fe for \n' ...
   '#mechanically selected 66,339 grain boundaries \n' ...
	'#Columns 1:8 : Euclidean grain boundary octonion coordinates \n' ...
	'#Column 9 : Grain boundary energy (J/m^2) \n' ...
	'#------------------------------------------------- \n']); %6 header lines
% writematrix([meshList propList],fname,'WriteMode','append','Delimiter','tab') %only works in 2020a
ftmp = 'temp.txt';
ftmppath = fullfile(folder,ftmp);
writematrix([meshList propList],ftmppath,'Delimiter','tab');
txt2 = fileread(ftmppath);
fprintf(fid,txt2);
fclose(fid);
%save to .mat file
save(fpath(1:end-4),'meshList','propList','five','meshTable')


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

%}
