function yq = interp_idw(X,qm2,nA2,y,r,L)
arguments
    X(:,8) double %input points
    qm2(:,4) double %query misorientations
    nA2(:,3) double %query BP normals
    y(:,1) double %property values
    r double = [] %radius
    L double = 2 % default is Euclidean norm
end
% INTERP_IDW interpolate using inverse-distance weighting and
% misorientation/BP normal query input pairs
Xq = get_pts(qm2,nA2);
yq = idw(X,Xq,y,r,L);

end