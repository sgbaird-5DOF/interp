%% Test from GB Octonion Tutorial script
% note: this requires functions from the GB_octonion_code
% repository, specifically the ax2qu.m function

addpathdir({'ax2qu.m','rot2q.m','GBdist.m','PGnames.mat'})

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

oAB = GBlab2oct(qA_Lab,qB_Lab,nA_Lab,'francis')

qC_Lab = rot2q(atan(1/2),pi/2,pi/2);
qD_Lab = rot2q(-atan(1/2),pi/2,pi/2);
nC_Lab = [0 0 1];

oCD = GBlab2oct(qC_Lab,qD_Lab,nC_Lab,'francis');

OmegaTest = rad2deg(GBdist([oAB,oCD],32,false))

OmegaTest2 = rad2deg(GBdist4(oAB,oCD,32,'omega'))

epsijk = 1;
euA=qu2eu(qAbp,epsijk);
euB=qu2eu(qBbp,epsijk);
euC=qu2eu(qCbp,epsijk);
euD=qu2eu(qDbp,epsijk);
qmconvention = 'johnson';
% [qm1,nA1] = eumA2five(euA,euB,nA_Lab,qmconvention)
% [qm2,nA2] = eumA2five(euC,euD,nC_Lab,qmconvention)

% o1 = GBfive2oct(qm1,nA1,qmconvention)
% o2 = GBfive2oct(qm2,nA2,qmconvention)

om1 = qu2om(qAbp);
om2 = qu2om(qBbp);
om3 = qu2om(qCbp);
om4 = qu2om(qDbp);

o1 = GBmat2oct(om1,om2);
o2 = GBmat2oct(om3,om4);

OmegaTest3 = rad2deg(GBdist4(o1,o2,32,'omega'))

om1a = qu2om(o1(1:4));
om2a = qu2om(o1(5:8));
om3a = qu2om(o2(1:4));
om4a = qu2om(o2(5:8));

[omcheck1,omcheck2] = constructGBMatrices(qA_Lab,qB_Lab,nA_Lab.','livermore');
GB5DOF(om1a,om2a,'Ni')
GB5DOF(omcheck1,omcheck2,'Ni')
o1a = GBmat2oct(om1a,om2a)
o2a = GBmat2oct(om3a,om4a)

OmegaTest4 = rad2deg(GBdist4(o1a,o2a,32,'omega'))

% [omcheck1a,omcheck2a] = constructGBMatrices(o1(1:4),o1(5:8),[0 0 1].','livermore');
% [omcheck1b,omcheck2b] = constructGBMatrices(o2(1:4),o2(5:8),[0 0 1].','livermore');
% 
% GB5DOF

qA = o1(1:4);
qB = o1(5:8);
qC = o2(1:4);
qD = o2(5:8);
%the two different statements for nA negate the second entry with respect to each
%other
conventionlist={'johnson','francis'};
for i=1:length(conventionlist)
    convention=conventionlist{i};
    switch convention
        case 'johnson'
            qm1 = qmult(qinv_francis(qA),qB)
            qm2 = qmult(qinv_francis(qC),qD)
        case 'francis'
            qm1 = qmult(qB,qinv_francis(qA))
            qm2 = qmult(qD,qinv_francis(qC))
    end
    switch convention
        case 'johnson'
            nA = qmult(qinv_francis(qA),qmult([0 0 0 1],qA));
            nC = qmult(qinv_francis(qC),qmult([0 0 0 1],qC));
        case 'francis'
            nA = qmult(qA,qmult([0 0 0 1],qinv_francis(qA)));
            nC = qmult(qC,qmult([0 0 0 1],qinv_francis(qC)));
    end
nA = nA(2:4)
nC = nC(2:4)
o1b = GBlab2oct([1 0 0 0],qm1,nA);
o2b = GBlab2oct([1 0 0 0],qm2,nC);

OmegaTest5 = rad2deg(GBdist4(o1b,o2b,32,'omega'))
end
%not sure why, but qm is the same between these two (just a fluke of the
%problem)


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