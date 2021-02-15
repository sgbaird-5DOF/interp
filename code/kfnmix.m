function kfn = kfnmix(kfntmp,kfntmp2,ypred,scl,thr)
arguments
    kfntmp
    kfntmp2
    ypred %ypredtmp3 from gprmix
    scl(1,1) double = 30
    thr(1,1) double = 1.1
end

% sigmoid function values
A = sigfn(ypred,scl,thr);
B = 1-A;

kfn = @(X,X2) A.*kfntmp(X,X2)+B.*kfntmp2(X,X2); %make X == X2 for square covariance matrix

end