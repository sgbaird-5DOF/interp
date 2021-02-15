% EGPRM_FAST_TEST
load('C:\Users\sterg\Documents\GitHub\posterior-sampling\figures\ensemble-interp-1000',...
    'mdllist','ypredlist');
n2 = 1000;
o2 = get_ocubo(1000);

mat = 'Ni';
epsijk = 1;
ytrue = GB5DOF_setup(o2(:,1:4),o2(:,5:8),[0 0 1],mat,epsijk);

thr = 1.1;
[ypost,ypred,ysd,ytrue,ci,covmat,kfntmp,kfntmp2,mdl,gprMdl2list,X2] = egprm(mdllist,o2,ytrue,thr);

npost = size(ypost,2);