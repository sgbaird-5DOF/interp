function kfn = kmix(kfntmp,kfntmp2)
arguments
    kfntmp
    kfntmp2
end
% KMIX  sigmoid-based covariance matrix mixture function
kfn = @(X1,X2,A) A.*kfntmp(X1,X2)+(1-A).*kfntmp2(X1,X2);
end