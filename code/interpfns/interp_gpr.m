function props = interp_gpr(gprMdl,qm2,nA2,projtol,usv)

ppts2 = get_ppts(qm2,nA2,projtol,usv);
props = predict(gprMdl,ppts2);

end

%% CODE GRAVEYARD
%{
%command to do interpolation
mdlcmd = @(gprMdl,ppts2) predict(gprMdl,ppts2);
interpfn = @(qm2,nA2) mdlcmd(compact(gprMdl),get_ppts(qm2,nA2));
%}