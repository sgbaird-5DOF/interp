function covmat = covmix(X2,ypred,kfn,scl,thr)
arguments
   X2
   ypred
   kfn
   scl(1,1) double = 30
   thr(1,1) double = 1.1
end

A = sigfn(ypred,scl,thr);
covmat = kfn(X2,A);