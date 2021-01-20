function testPosteriorCovarianceGP()
% 1. Some dummy data.
rng(0,'twister');
X = [-10,-8,-5,-1,1,3,7,10]';
N = numel(X);
noisescl = 1e-4;
y = 1 - X*5e-2 + sin(X)./X + noisescl*randn(N,1);

M = 1000;
Xnew = linspace(-10,10,M)';

% 2. Fit a GPR model.
gpr = fitrgp(X,y,'Sigma',0.01,'ConstantSigma',true,'SigmaLowerBound',1e-2);
% gpr = fitrgp(X,y);

% 3. Access the posterior covariance matrix by calling the undocumented
% predictExactWithCov method on the gpr.Impl object. Vectors ynew can be
% simulated from a Normal distribution with mean pred and covariance covmat
% using cholcov or svd and randn.
alpha = 0.05;
[pred,ysd,ci] = predict(gpr,Xnew);
kfcn = gpr.Impl.Kernel.makeKernelAsFunctionOfXNXM(gpr.Impl.ThetaHat);
covmat = kfcn(Xnew,Xnew);
% [pred,covmat,ci] = predictExactWithCov(gpr.Impl,Xnew,alpha);

paperfigure(1,2);
nexttile
% plot(Xnew,pred,'b');
% hold on;
% plot(Xnew,ci(:,1),'m');
% plot(Xnew,ci(:,2),'g');
shadedErrorBar(Xnew,pred,ysd)
hold on
plot(X,y,'ko','MarkerFaceColor','k');
xlabel('x');
ylabel('y');

nexttile
% T = cholcov(covmat);
numSamples = 50;
% ynew = pred + T'*randn(M,numSamples);
ynew = mvnrnd(pred,covmat,numSamples);
plot(Xnew,ynew);
hold on;
plot(X,y,'ko','MarkerFaceColor','k');
xlabel('x');
ylabel('y');
end