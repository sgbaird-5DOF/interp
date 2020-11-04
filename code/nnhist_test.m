%nnhist test
clear; close all
seed = 10;
rng(seed);
addpathdir({'sqrt2norm.m','normr.m','get_omega.m'})
pts = sqrt2norm(normr(rand(1000,8)));
dtype = 'omega';
[ax,D,nnpts,idx] = nnhist(pts,dtype);