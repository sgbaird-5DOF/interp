function rnew = Lpr(p,r,epsijk)
arguments
	p(:,4) double {mustBeReal,mustBeFinite}
	r(:,3) double {mustBeReal,mustBeFinite}
    epsijk(1,1) double = 1
end
% LPR quaternion rotation active (epsijk==1) and passive (epsijk==-1)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-07-27
% 
% Inputs:
%  qmis - misorientation quaternion
%  r - boundary plane normal in cartesian coordinates
%
% Outputs:
%  mA - rotated vector (i.e. r rotated by qmis)
%
% Usage:
%  mA = Lpr(qmis,r);
%  
% References:
%  https://dx.doi.org/10.1088/0965-0393/23/8/083501 (eqn 24)
%--------------------------------------------------------------------------
%definitions
p0 = p(:,1);
p = p(:,2:4);
pmag = vecnorm(p,2,2);

%equation
rnew = (p0.^2-pmag.^2).*r+2*dot(p,r,2).*p+2*epsijk*p0.*cross(p,r,2);

end