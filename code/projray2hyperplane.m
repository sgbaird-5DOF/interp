function a = projray2hyperplane(nvec,pt)
% PROJRAY2HYPERPLANE  Project ray (pt) from unit hypersphere to tangent hyperplane at another point (nvec)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
% 
% Inputs:
%	nvec	=== normalized normal to hyperplane
%
%	pt		=== point to project from hypersphere to hyperplane
%
% Intermediates:
%
%	d		===	dimensionality
%
%	dist	===	distance from origin to hyperplane (same as radius of
%					hypersphere)
%
%	t		===	parametric value that scales pt to the intersection of pt
%					and tangent hyperplane.
%
% Outputs:
%
%	a		===	intersection of pt and tangent hyperplane defined by nvec
%
% References:
%		https://math.stackexchange.com/q/1256236/798661
%
%		http://comprna.upf.edu/courses/Master_MAT/3_Optimization/U9_Hyperplanes.pdf
%--------------------------------------------------------------------------

%dimensionality
d = length(nvec);

%distance
dist = norm(nvec);  

%parametric value
t = dist/dot(nvec,pt);

%intersection
a = -pt*(-t);

end





	