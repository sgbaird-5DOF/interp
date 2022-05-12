% Sigma 5 disorientation
qm = rot2q(deg2rad(36.86),deg2rad(90),0);
qm = [qm; qm];

% BP normals
nA = [0.707106781186548  0.590122694337284  0.389557705132506;... % an orange point
    0.707106781186548 -0.0423055217160760 0.705840097212060]; % a blue point

% convert to octonions
o = five2oct(qm,nA,1)

% symmetrize
oref = get_orefs(1); % VFZ reference octonion
osymm = get_octpairs(o,'pgnum',32, 'oref',oref)

GBfive2oct()