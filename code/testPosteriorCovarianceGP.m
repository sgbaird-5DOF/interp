function testPosteriorCovarianceGP()
% 1. Some dummy data.
rng(0,'twister');
X = [-15,-10,-8,-5,-1,1,3,7,10,15,20]';
N = numel(X);
noisescl = 1e-4;
y = 1 - X*5e-2 + sin(X)./X + noisescl*randn(N,1);

M = 1000;
Xnew = linspace(min(X),max(X),M)';

% 2. Fit a GPR model.
% gpr = fitrgp(X,y,'Sigma',0.01,'ConstantSigma',true,'SigmaLowerBound',1e-2);
gpr = fitrgp(X,y);

% 3. Access the posterior covariance matrix by calling the undocumented
% predictExactWithCov method on the gpr.Impl object. Vectors ynew can be
% simulated from a Normal distribution with mean pred and covariance covmat
% using cholcov or svd and randn.
alpha = 0.05;
[pred,ysd,ci] = predict(gpr,Xnew);
kfcn = gpr.Impl.Kernel.makeKernelAsFunctionOfXNXM(gpr.Impl.ThetaHat);
covmat = kfcn(Xnew,Xnew);
covmat = nearestSPD(covmat);
covmat = nearestSPD(covmat);
% [pred,covmat,ci] = predictExactWithCov(gpr.Impl,Xnew,alpha);

paperfigure(1,3);
i = 0;
nexttile
i = i+1;
% plot(Xnew,pred,'b');
% hold on;
% plot(Xnew,ci(:,1),'m');
% plot(Xnew,ci(:,2),'g');
shadedErrorBar(Xnew,pred,ysd)
hold on
plot(X,y,'ko','MarkerFaceColor','k');
xlabel('x');
ylabel('y');
axis tight square
papertext(i)
% T = cholcov(covmat);
n = 50;

methodlist = {'mvn','tmvn'};
for j = 1:length(methodlist)
    i = i+1;
    nexttile
    method = methodlist{j};
    % method = 'tmvn'; %'mvn', 'tmvn'
    switch method
        case 'mvn'
            % ynew = pred + T'*randn(M,n);
            % sig = mvnrnd(zeros(size(pred)).',covmat,n).';
            % % sig = exp(sig);
            % ynew = repmat(pred,1,n)+sig;
            ynew = mvnrnd(pred,covmat,n);
        case 'tmvn'
            % truncated normal distribution
            l = max([zeros(size(pred)),0.75*pred],[],2)-pred;
            u = 1.5*pred;
            zerofloorQ = true;
            ynew = tmvn(pred,covmat,l,u,n,zerofloorQ);
    end
    
    plot(Xnew,ynew);
    hold on;
    plot(X,y,'ko','MarkerFaceColor','k');
    axis tight square
    xlabel('x');
    ylabel('y');
    papertext(i)
end

end