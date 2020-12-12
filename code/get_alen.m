function alen = get_alen(pts,pts2)
arguments
    pts double
    pts2 double
end
% GET_ALEN  get arclength of sphere (acos of dot product, note that this is half the octonion distance)
npts = size(pts2,1);
if size(pts,1) == 1
    pts = repmat(pts,npts,1);
end
alen = acos(dot(pts,pts2,2));
thr = 1e-6;
maxim = max(abs(imag(alen)));
assert(maxim <= thr,['max imaginary component alen > ' num2str(thr) ', == ' num2str(maxim)])
alen = real(alen);
end