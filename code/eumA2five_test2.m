% EUMA2FIVE_TEST2
function tbl = eumA2five_test2()
%% Lab to Five
npts = 1;
epsijk = 1;
pgnum = 32;
rng(11)
plotQ = false;

type = 'rand'; %'rand', 'sym'

switch type
    case 'rand'
        %random data, more robust for error checking
        eA = normr(rand(npts,3));
        eB = normr(rand(npts,3));
        mA = normr(rand(npts,3));
        pA = eu2qu(eA,epsijk);
        pB = eu2qu(eB,epsijk);
        
    case 'sym'
        %higher symmetry, easier to visualize
        pA = ax2qu([0 0 1 deg2rad(20)],epsijk);
        pB = ax2qu([0 0 1 deg2rad(-10)],epsijk);
        eA = qu2eu(pA,epsijk);
        eB = qu2eu(pB,epsijk);
        mA = normr([cos(deg2rad(-20)),0.5,sin(deg2rad(-20))]);
end

gbe0 = GB5DOF_setup(pA,pB,mA,'Ni',epsijk);

[qm,nA] = eumA2five(eA,eB,mA,epsijk);
gbe = GB5DOF_setup([1 0 0 0],qm,nA,'Ni',epsijk);
qmdis = qdis(qm,pgnum,epsijk);

%% Plotting 5DOF conversion
if plotQ
    paperfigure(1,2);
    nexttile(1)
    title('Standard xyz to misorientation quaternion (qm)')
    rotquivplot([1 0 0 0],qm(1,:),0.1,1,epsijk);
    % rotquivplot(qA(1,:),qB(1,:),0.1,1,epsijk);
    nexttile(2)
    title('rotation of lab frame BP normal (mA) into grain A')
    rotvecplot([1 0 0 0],pA(1,:),mA(1,:),0.1,1,epsijk);
    
    % single quiver from [0 0 0] to nA(1,:) (i.e. the true direction it should map to)
    w = [0 0 0];
    wtmp=n2c(w);
    t=n2c(0.5*nA(1,:));
    quiver3(wtmp{:},t{:},0,'k','linewidth',1,'Autoscale','off');
end
%% Octonion Conversion
R = vecpair2rmat(mA(1,:),[0 0 1],1); %not sure why this has to stay in active interpretation for epsijk==1 *and* epsijk==-1

qR = om2qu(R,epsijk);

qA1 = qmult(qR,pA,epsijk);
qB1 = qmult(qR,pB,epsijk);

omA = qu2om(qA1,epsijk);
omB = qu2om(qB1,epsijk);

gbe1 = GB5DOF_setup(qA1,qB1,[0 0 1],'Ni',epsijk);

%plotting
if plotQ
    paperfigure(1,2);
    nexttile(1)
    title('lab frame [1 0 0 0] to octonion qA')
    rotquivplot([1 0 0 0],qA1,0.1,1,epsijk);
    nexttile(2)
    title('lab frame qm to octonion qB')
    rotquivplot(qm,qB1,0.1,1,epsijk);
end

[qm1,nA1] = qmA2five(qA1,qB1,[0 0 1],epsijk);
qmdis1 = qdis(qm1,pgnum,epsijk);

gbe1a = GB5DOF_setup([1 0 0 0],qm1,nA1,'Ni',epsijk);

% GBmat2oct
o = GBmat2oct(omA,omB,epsijk); %note that qA == qA1 and qB == qB1

qA = o(1:4);
qB = o(5:8);
[qm2,nA2] = qmA2five(qA,qB,[0 0 1],epsijk);
qmdis2 = qdis(qm2,pgnum,epsijk);

gbe2 = GB5DOF_setup(qA,qB,[0 0 1],'Ni',epsijk);

gbe2a = GB5DOF_setup([1 0 0 0],qm2,nA2,'Ni',epsijk);

%% Plotting Intermediate to Octonion
if plotQ
    paperfigure(1,2);
    nexttile(1)
    title('misorientation from octonion qA to qB')
    rotquivplot(qA1,qB1,0.1,1,epsijk);
    nexttile(2)
    title('grain A BP normal (lab) rotated to [0 0 1]')
    rotvecplot([1 0 0 0],qR,mA(1,:),0.1,1,epsijk);
    
    % single quiver from [0 0 0] to [0 0 1] (i.e. the true direction it should map to)
    w = [0 0 0];
    wtmp=n2c(w);
    t=n2c(0.5*[0 0 1]);
    quiver3(wtmp{:},t{:},0,'k','linewidth',1,'Autoscale','off');
end
%% Symmetry Checks
o2 = get_ocubo();
% [~,osym] = GBdist([o2 o],pgnum,false,epsijk);
% osym1 = osym(1:8);
% osym2 = osym(9:16);
waitbarQ = true;
[~, osym2] = GBdist4(o2,o,32,'norm',1e-12,waitbarQ,epsijk);
osym2 = osym2{1};

qAsym = osym2(1,1:4);
qBsym = osym2(1,5:8);
[qmsym,nAsym] = qmA2five(qAsym,qBsym,[0 0 1],epsijk);

gbesym = GB5DOF_setup(qAsym,qBsym,[0 0 1],'Ni',epsijk);
gbesyma = GB5DOF_setup([1 0 0 0],qmsym,nAsym,'Ni',epsijk);

%when active-active (qmult/qm) used, qmdis==qmdis2, qmdis2~=qmsymdis
%when active-passive (qmult/qm) used, qmdis~=qmdis2, qmdis2==qmsymdis
[qmsymdis,negQ,invQ,disSpair] = qdis(qmsym,pgnum,epsijk);

qmdata = [qmdis;qmdis1;qmdis2;qmsymdis];
nAdata = [nA;nA1;nA2;nAsym];
gbedata = [gbe0;gbe1;gbe2;gbesym];
gbedata2 = [gbe;gbe1a;gbe2a;gbesyma];
tbl = table(qmdata,nAdata,gbedata,gbedata2,'VariableNames',{'qm','nA','lab-gbe','five-gbe'},'RowNames',{'start','vecpair2rmat','GBmat2oct','symmetrized'});
disp(tbl)

%% Octonion distance check
[dmin,o2minsym] = GBdist4(o,osym2);
o2minsym = o2minsym{1};
disp(['dmin: ',num2str(dmin)])

%% checking BP symmetries
Spairs = get_sympairs(pgnum);
Spairs = [Spairs; -Spairs];
SA = Spairs(:,1:4);
SB = Spairs(:,5:8);
nsym = size(Spairs,1);
qmsymrep = repmat(qmsym,nsym,1);
qcheck = qmult(qinv(SA),qmult(qmsymrep,SB));
qIDs = find(ismembertol(qcheck,qmsym,1e-3,'ByRows',true));
end

%% CODE GRAVEYARD
%{
% if epsijk == 1
%     tmp = (R.'*(nA(1,:)).').';
% else
%     tmp = (R*(nA(1,:)).').';
% end
% assert(all(ismembertol(tmp,[0 0 1])),'did not map back to [0 0 1]')

% om1 = qu2om([1 0 0 0],epsijk);
% om2 = qu2om(qm,epsijk);
% 
% omA = R.'*om1(:,:,1); %should it be this or R*om1(:,:,1)?
% omB = R.'*om2(:,:,1); %should it be this or R*om2(:,:,1)?
% 
% qA1 = om2qu(omA,epsijk);
% qB1 = om2qu(omB,epsijk);

% nAsym = Lpr(qAsym,[0 0 1],epsijk);
% nBsym = Lpr(qBsym,[0 0 1],epsijk);

% qmsym = qlab2qm(qAsym,qBsym,epsijk);

% R = vecpair2rmat(nA(1,:),[0 0 1]);

% qA1 = qmult(qR,[1 0 0 0],epsijk);
% qB1 = qmult(qR,qm,epsijk);

% nA1 = Lpr(qA1,[0 0 1],epsijk);
% nA2 = Lpr(qA,[0 0 1],epsijk);
%}