function ppts = get_ppts(pts,projtol,usv,zeroQ)
% GET_PPTS  get projected points via proj_down.m
ppts = proj_down(pts,projtol,usv,'zeroQ',zeroQ);
end