function props = interp_gpr(gprMdl,qm2,nA2,projtol,usv,zeroQ)
% INTERP_GPR  interpolate using Gaussian Process Regression model (gprMdl) and predict()
pts2 = get_pts(qm2,nA2);
% ppts2 = get_ppts(pts2,projtol,usv,zeroQ);
props = predict(gprMdl,pts2);

end

%% CODE GRAVEYARD
%{
%command to do interpolation
mdlcmd = @(gprMdl,ppts2) predict(gprMdl,ppts2);
interpfn = @(qm2,nA2) mdlcmd(compact(gprMdl),get_ppts(qm2,nA2));
%}