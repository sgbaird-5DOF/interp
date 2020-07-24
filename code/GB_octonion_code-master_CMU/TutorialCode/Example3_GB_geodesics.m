%% 1. LOAD OCTONION UTILITY FUNCTIONS 

addpath([fileparts(pwd),'/Data']) %add Data directory to path
addpath('crystal_symmetry_ops')
addpath('octonion_functions/')
addpath('rotation_conversions/')

%% PROBLEM: interpolate between two grain boundaries along geodesic

% IMPORTANT: By convention, the GB plane normal needs to fall along the z direction 
% IMPORTANT: octonions need to be symmetrized before computing OSLERP trajectory

test = importdata('../Data/olm_octonion_list.txt',' ',1); %list of GB octonions with number of octonions as first line in file
data = test.data;

id1 = 1;
id2 = 3;

o1i = data(id1,:); %unsymmetrized octonions
o2i = data(id2,:);

%%%%%%%%%%%%%%% Symmetrization routine for cubic point group symmetry %%%%%%%%%%%%%%%
pgnum = 30;

[omega,oct_symmetrized,~] = GBdist([o1i o2i],pgnum,false);
disp(['geodesic distance (degrees) between GB ids ',num2str(id1),' and ',num2str(id2)])
disp(rad2deg(omega))

o1 = oct_symmetrized(1:8); %symmetrized octonions
o2 = oct_symmetrized(9:16); 

%% perform OSLERP on symmetrized octonions

nt = 10; %number of interpolated points
oslerp_example = OSLERP(o1,o2,omega,nt);

%check that octonions have the same norm along trajectory: 
disp(['OSLERP trajectory connecting GBs ',num2str(id1),' and ',num2str(id2)])
disp(oslerp_example)
disp('norms (should be sqrt(2)):')
disp(vecnorm(oslerp_example,2,2))






