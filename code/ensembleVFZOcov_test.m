% ENSEMBLEVFZOCOV_TEST

%rng
seed = 11;
rng(seed);

%generate octonions
o = get_ocubo(100);
% o2 = get_ocubo(100);
o2 = o;

%symmetrize octonions
o = get_octpairs(o);
o2 = get_octpairs(o2);

%reduce dimensionality
[Xcat,usv] = proj_down([o;o2]);
tol = 1e-6;
X = proj_down(o,tol,usv);
X2 = proj_down(o2,tol,usv);

%kernel parameters
theta = [1,1];

%get covariance matrix
cov = ensembleVFZOcov(X,X2,theta,usv);

%plot
paperfigure();
imagesc(cov);
axis equal