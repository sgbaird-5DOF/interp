function props = interp_nn(ppts,qm2,nA2,projtol,usv,propList)

if isempty(usv)
    error('usv should not be empty')
end
props = propList(dsearchn(ppts,get_ppts(qm2,nA2,projtol,usv)));

end

%% CODE GRAVEYARD
%{
mdlcmd = @(ppts,ppts2,propList) propList(dsearchn(ppts,ppts2));
interpfn = @(qm2,nA2) mdlcmd(ppts,get_ppts(qm2,nA2),propList);
%}