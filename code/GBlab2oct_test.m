%% Test from GB Octonion Tutorial script
% note: this requires functions from the GB_octonion_code
% repository, specifically the ax2qu.m function

addpathdir('ax2qu.m')

aaAbp = [0 1 0 atan(1/5)]; % grain orientations in GB plane reference frame
aaBbp = [0 1 0 -atan(1/5)];
qAbp = ax2qu(aaAbp); qBbp = ax2qu(aaBbp); % axis angle pair --> quaternion

o1 = 1/sqrt(2)*[qAbp qBbp];

aaCbp = [0 1 0 atan(1/2)]; 
aaDbp = [0 1 0 -atan(1/2)];
qCbp = ax2qu(aaCbp); qDbp = ax2qu(aaDbp);

o2 = 1/sqrt(2)*[qCbp qDbp];

Omega = rad2deg(2*acos(dot(o1,o2)))

%% Test using our functions with GBdist
qA_Lab = rot2q(atan(1/5),pi/2,pi/2);
qB_Lab = rot2q(-atan(1/5),pi/2,pi/2);
nA_Lab = [0 0 1];

oAB = GBlab2oct(qA_Lab,qB_Lab,nA_Lab,'francis');

qC_Lab = rot2q(atan(1/2),pi/2,pi/2);
qD_Lab = rot2q(-atan(1/2),pi/2,pi/2);
nC_Lab = [0 0 1];

oCD = GBlab2oct(qC_Lab,qD_Lab,nC_Lab,'francis');

OmegaTest = rad2deg(GBdist([oAB,oCD],32,false))

OmegaTest2 = rad2deg(GBdist4(oAB,oCD,32,'omega'))


%----------------------------CODE GRAVEYARD--------------------------------
%{

% note: this requires functions from Dr. Johnson's quaternion library, AND 
% you should remove the GB_octonion_code-master from the search path since 
% it has some files with the same name as those required here)

% addpath(genpath('GB_octonion_code-master_CMU'))
% rmpath(genpath('GB_octonion_code-master')) % remove from path so that GBlab2oct uses Dr. Johnson's qinv function
% addpath(genpath('GB_octonion_code-master')) % add back to path to use GBdist

nA_Lab = [0 0 1].';

oAB = GBlab2oct(qA_Lab,qB_Lab,nA_Lab(:),'francis');

nC_Lab = [0 0 1].';

oCD = GBlab2oct(qC_Lab,qD_Lab,nC_Lab(:),'francis');
%}