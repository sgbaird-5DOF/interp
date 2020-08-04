%optimize_zeta test
clear; close all
addpathdir({'misFZfeatures.mat','PGnames.mat','qu2ax.m','q2rod.m'})

seed = 10;
rng(seed);

npts = 2;
disp(['npts == ' int2str(npts)])
o = get_ocubo_r2018a(npts,'random',double.empty,seed);

% S = load('octvtx_pairmin.mat');
% o = S.octvtx;

savename = 'temp.mat';
NV.o2addQ = false;
NV.plotQ = false;
NV.method = 'standard';
NV.pgnum = 32;
NV.wtol = 1e-6;

o = get_octpairs_r2018a(o,savename,NV);
% npts = size(o,1);

dtype = 'omega';

[z,errmin,exitflag,output] = optimize_zeta_r2018a(o,dtype,[]);

disp(['errmin == ' num2str(rad2deg(errmin)) ' degrees'])

IDs = allcomb(1:npts,1:npts);
o1 = o(IDs(:,1),:);
o2 = o(IDs(:,2),:);
% "true" pairwise distances
pdtrue = GBdist4_r2018a(o1,o2,32,dtype,1e-6).';

err0 = pd_sse_r2018a(o,zeros(npts,1),dtype,pdtrue,function_handle.empty,'me');

disp(['err0 - errmin == ' num2str(rad2deg(err0-errmin)) ' degrees'])


%----------------------------CODE GRAVEYARD--------------------------------

