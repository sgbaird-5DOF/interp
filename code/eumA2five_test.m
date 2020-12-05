% EUMA2FIVE_TEST
clear; close all

addpathdir({'eu2qu.m','qinv.m'})
%% Conversion
npts = 1;
rng(12)
eA = normr(rand(npts,3));
eB = normr(rand(npts,3));
mA = normr(rand(npts,3));
epsijk = 1;

pA = eu2qu(eA,epsijk);
pB = eu2qu(eB,epsijk);

[qm,nA,nB] = eumA2five(eA,eB,mA,epsijk);

%% Plotting 5DOF conversion
paperfigure(1,2);
nexttile(1)
rotquivplot([1 0 0 0],qm(1,:),0.1,1,epsijk);
% rotquivplot(qA(1,:),qB(1,:),0.1,1,epsijk);
nexttile(2)
rotvecplot([1 0 0 0],pA(1,:),mA(1,:),0.1,1,epsijk);
w = [0 0 0];
wtmp=n2c(w);
t=n2c(0.5*nA(1,:));
quiver3(wtmp{:},t{:},0,'k','linewidth',1,'Autoscale','off');

%% Octonion Conversion
R = vecpair2rmat(nA(1,:),[0 0 1]);
% RB = vecpair2rmat(nB(1,:),[0 0 1]);
tmp=(R*(nA(1,:)).').';
assert(all(ismembertol(tmp,[0 0 1])),'did not map back to [0 0 1]')
% om1 = qu2om(pA,epsijk);
% om2 = qu2om(pB,epsijk);
om1 = qu2om([1 0 0 0],epsijk);
om2 = qu2om(qm,epsijk);
% om1 = {q2gmat(pA)};
% om2 = {q2gmat(pB)};

% omA = om1{1}*R;
% omB = om2{1}*R;
omA = R.'*om1(:,:,1);
omB = R.'*om2(:,:,1);
% omA = (R*om1{1}.').';
% omB = (R*om2{1}.').';

qR = om2qu(R,epsijk);

qA1 = om2qu(omA,epsijk);
qB1 = om2qu(omB,epsijk);

paperfigure(1,2);
nexttile(1)
rotquivplot(pA,qA1,0.1,1,epsijk);
nexttile(2)
rotquivplot(pB,qB1,0.1,1,epsijk);

qm1 = qlab2qm(qA1,qB1,epsijk)
nA1 = Lpr(qm1,[0 0 1],epsijk)

o = GBmat2oct(omA,omB,epsijk)

qA = o(1:4);
qB = o(5:8);
qm(1,:)
qm2 = qlab2qm(qA,qB,epsijk)

nA2 = Lpr(qm2,[0 0 1],epsijk)

%% Plotting Intermediate to Octonion
paperfigure(1,2);
nexttile(1)
rotquivplot(qA1,qB1,0.1,1,epsijk);
% rotquivplot(qA(1,:),qB(1,:),0.1,1,epsijk);
nexttile(2)
rotvecplot([1 0 0 0],qR,nA(1,:),0.1,1,epsijk);
w = [0 0 0];
wtmp=n2c(w);
t=n2c(0.5*[0 0 1]);
quiver3(wtmp{:},t{:},0,'k','linewidth',1,'Autoscale','off');

%% Symmetry Checks
% osym = get_octpairs(o);
% GBdist4(o,osym,32)

o2 = get_ocubo();
[~,osymorig] = GBdist([o o2],30);
osymorig1 = osymorig(1:8);
osymorig2 = osymorig(9:16);

qAsym = osym(1,1:4);
qBsym = osym(1,5:8);
qmsym = qlab2qm(qAsym,qBsym,epsijk)
% qmdis = disorientation(qm)
% qmsymdis = disorientation(qmsym)
nAsym = Lpr(qmsym,[0 0 1],epsijk)




