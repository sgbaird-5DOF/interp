function ppts = get_ppts(qm,nA,projtol,usv)
ppts = proj_down(get_octpairs(GBfive2oct(qm,nA)),projtol,usv);
end