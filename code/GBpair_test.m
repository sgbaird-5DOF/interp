%GBPAIR_TEST  test 
fname = {'misFZfeatures.mat'};
addpathdir(fname)

load('misFZfeatures.mat','qlist')

%% get two reference octonions (o1 and o2)
%interior point
name1 = 'interior';
%no boundary point
name2 = 'O';
%other high symmetry point
name3 = 'C';

%unpack quaternions
qA = qlist.(name1);
qB = normr(qlist.(name2));
qC = qlist.(name3);

%load normals (both are arbitrary set to [0 0 1])
[~,RA] = symaxis(qA,name1);
nA = (RA*[0 0 1].').';

[~,RB] = symaxis(qB,name2);
nB = normr((RB*[0 0 1].').');

[~,RC] = symaxis(qC,name3);
nC = normr((RC*[0 0 1].').');

o1 = GBfive2oct(qA,nA);
o2 = GBfive2oct(qB,nB);
o3 = GBfive2oct(qC,nC);
% o3 = sqrt2norm(round(get_ocubo(1),2),'quat');

t=n2c([o1;o2;o3]);
disp(table(t{:},'RowNames',{'o1','o2','o3'}))

[o3_out,omega3,o2_out] = GBpair(o1,o2,o3);