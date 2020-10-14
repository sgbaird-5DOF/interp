function ppts = get_ppts(pts,projtol,usv,zeroQ)
ppts = proj_down(pts,projtol,usv,'zeroQ',zeroQ);
end