% mvrandn_test
Sig = [0.1 0.01; 0.01, 0.05];
m = [0.3 0.5];
l = -1*ones(1,2);
u = ones(1,2);
n = 10000;
rv = mvrandn(l,u,Sig,n);
paperfigure();
hexscatter(rv(1,:).',rv(2,:).','reflineQ',false)

% paperfigure();
% d=3;n=10^3;Sig=0.9*ones(d,d)+0.1*eye(d);l=(1:d)/d*4;u=l+2;
% X=mvrandn(l,u,Sig,n);boxplot(X','plotstyle','compact') % plot marginals
% 
% paperfigure();
