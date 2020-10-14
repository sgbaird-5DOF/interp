function alen = get_alen(pts,pts2)
arguments
    pts double
    pts2 double
end
% GET_ALEN  get arclength of sphere (acos of dot product)
npts = size(pts2,1);
if size(pts,1) == 1
    pts = repmat(pts,npts,1);
end
alen = acos(dot(pts,pts2,2));
end