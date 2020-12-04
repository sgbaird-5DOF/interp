
%% 1. LOAD OCTONION UTILITY FUNCTIONS 

addpathdir({'olmsted_xtal_info_numeric.csv','crystal_symmetry_ops',...
    'rotation_conversions','octonion_functions'})

% addpath([fileparts(pwd),'/Data']) %add Data directory to path
% addpath('**/crystal_symmetry_ops')
% addpath('**/octonion_functions/')
% addpath('**/rotation_conversions/')

%% PROBLEM 0: converting GB's from conventional representations to octonion representation 

% IMPORTANT: By convention, the GB plane normal needs to fall along the z direction 
% GBmat2oct(O1,O2): orientation matrices in respective crystal frames to octonion in boundary plane ref frame
% GBfive2oct(misorientation,n): (misorientation, boundary plane inclination) pair to octonion in boundary plane ref frame 

% To read more about the theory, view supplementary info of the first octonion paper

% We will start by considering two examples from the first GB octonion paper. 

%% Example 1. Define two symmetric tilt grain boundaries o1 = (qA,qB) and o2 = (qC,qD) %%%%%%
% this example highlights the meaning of the "boundary plane reference frame" 

% tilt boundaries in this case are defined as rotations relative to fixed boundary plane with normal <001> 

aaAbp = [0 1 0 atan(1/5)]; % grain orientations in GB plane reference frame
aaBbp = [0 1 0 -atan(1/5)];
qAbp = ax2qu(aaAbp); qBbp = ax2qu(aaBbp); % axis angle pair --> quaternion

o1 = 1/sqrt(2)*[qAbp qBbp];

aaCbp = [0 1 0 atan(1/2)]; 
aaDbp = [0 1 0 -atan(1/2)];
qCbp = ax2qu(aaCbp); qDbp = ax2qu(aaDbp);

o2 = 1/sqrt(2)*[qCbp qDbp];

Omega = rad2deg(2*acos(dot(o1,o2))); %in this case, the GBO misorientation angle is simply the difference in tilt angles!


%% Example 2: traditional BP / misorientation to octonion

format long 

qmis = [sqrt(2/3) 1/sqrt(6)*[1 1 0]]; %misorientation quaternion
nA = normr([3 1 2]);
nB = normr([1 3 2]);
nC = normr([7 -1 2]);
nD = normr([3 3 6]); %pair of BP normals in crystal frame of each grain

% It suffices only to consider the boundary planes of the upper grains (nA,nC): 

o1 = GBfive2oct(qmis,nA);
o2 = GBfive2oct(qmis,nC);

% Now let's look at the geodesic distance with cubic crystal symmetry
pgnum = 30;
[Omega_ex2, ~, zeta_min] = GBdist([o1 o2],30,false);

disp('minimizing U(1) angle (degrees):')
disp(rad2deg(zeta_min))

disp('minimizing GBOM angle (degrees):')
disp(rad2deg(Omega_ex2))

%% Example 3: conversion of orientation matrices to octonion
% here we convert 388 Olmsted GB's to octonions
% In Olmsted survey, BP is along x by convention
% Need to rotate BP from x --> z

% import Olmsted dataset
olmimp = importdata('olmsted_xtal_info_numeric.csv');
olmx = olmimp.data; ngb = length(olmx);
feature_names = olmimp.textdata;

% extract indices of orientation matrices as N_GB x 9 matrices
O1mat = olmx(:,[5:7 11:13 17:19]);
O2mat = olmx(:,[8:10 14:16 20:22]);

octlist = zeros(ngb,8);
for i = 1:ngb
    % O1 and O2 should have hkl directions along columns (x,y,z)
    O1 = reshape(O1mat(i,:),[3,3]);
    O2 = reshape(O2mat(i,:),[3,3]);
    
    aa_z = [0 1 0 pi/2];
    om_z = ax2om(aa_z); %rotation matrix, BP x --> z 
    
    OA = (om_z*O1')'; %rotate row-wise, transpose to column form
    OB = (om_z*O2')';
    
    oct = GBmat2oct(OA,OB);
    octlist(i,:) = oct;
end

disp([num2str(ngb),' GBs converted to octonions'])

%%%%%%%% To output octonions to a text file, uncomment below %%%%%%%

format long

fname = 'olmsted_octonions.txt';
fID = fopen(fname,'w');
fprintf(fID,'oct \n');
fprintf(fID,[num2str(ngb),' \n']);
fprintf(fID,'%6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f \n', octlist');





%% CODE GRAVEYARD
%{
filelist = ;
filepathgen = fullfile('**',filelist);
for i = 1:length(filepathgen)
	octfile = dir(filepathgen{i});
	if ~isempty(octfile)
		octfolder = octfile(1).folder;
		addpath(octfolder);
	end
end
%}