function [ximin,o2sym] = xi_min(o1,o2,epsijk)
arguments
	o1(:,8) double {mustBeFinite,mustBeReal,mustBeSqrt2Norm}
	o2(:,8) double {mustBeFinite,mustBeReal,mustBeSqrt2Norm}
    epsijk(1,1) double
end
% XI_MIN  Correction to CMU group function zeta_min(), written
% by Oliver Johnson & Sterling Baird
%--------------------------------------------------------------------------
% Date: 2024-01-25
% 
% Inputs:
%		(o1,o2)	- lists of octonions
%       epsijk - scalar indicating active (1) or passive (-1) rotation
%                convention
%
% Outputs:
%		ximin	- list of minimized xi angles
%
% Usage:
%		ximin = xi_min(o1,o2);
%
% Dependencies:
%		*
%
% Notes:
%		* Eqs. 28-29 in the Octonion paper [1] are incorrect, and are 
%         corrected herein.
%
% References:
%       [1] Francis, T., Chesser, I., Singh, S., Holm, E. A., & 
%           de Graef, M. (2019). A geodesic octonion metric for grain 
%           boundaries. Acta Materialia, 166, 135–147. 
%           https://doi.org/10.1016/j.actamat.2018.12.034
%--------------------------------------------------------------------------

%unpack quaternions
qA = o1(:,1:4);
qB = o1(:,5:8);
qC = o2(:,1:4);
qD = o2(:,5:8);

% compute numerator
switch epsijk
    case 1 % active
        v4 =  ( qA(:,4).*qC(:,1) - qA(:,1).*qC(:,4) ) + ( qA(:,3).*qC(:,2) - qA(:,2).*qC(:,3) ) +...
             ( qB(:,4).*qD(:,1) - qB(:,1).*qD(:,4) ) + ( qB(:,3).*qD(:,2) - qB(:,2).*qD(:,3) );
    case -1 % passive
        v4 = ( qA(:,4).*qC(:,1) - qA(:,1).*qC(:,4) ) - ( qA(:,3).*qC(:,2) - qA(:,2).*qC(:,3) ) +...
            ( qB(:,4).*qD(:,1) - qB(:,1).*qD(:,4) ) - ( qB(:,3).*qD(:,2) - qB(:,2).*qD(:,3) );
    otherwise
        error('epsijk must be either 1 (for active rotations) or -1 (for passive rotations).')
end

% compute denominator
v1 = dot(qA,qC,2) + dot(qB,qD,2);

% ensure unit quaternion condition is satisfied
% assert(all(v4.^2 + v1.^2 <= 4));
prec = 6;
assert(all( nlt(v4.^2 + v1.^2,4,prec) | neq(v4.^2 + v1.^2,4,prec) )); % test of v4.^2 + v1.^2 <= 4 to a precision of "prec"

% compute ximin
ximin = 2*atan2(v4,v1);

% put in [0,4*pi]
% NOTE: this is not strictly necessary, but is here to help people avoid incorrectly doing mod(ximin,2*pi) which would give incorrect results
ximin = mod(ximin,4*pi);

if nargout == 2
    % apply the rotation to get the o2s that give U(1) minimized distance
    qxizs = [cos(ximin/2),zeros(numel(ximin),1),zeros(numel(ximin),1),sin(ximin/2)];
    switch epsijk
        case 1 % active
            o2sym = [qmultiply(qxizs,qC), qmultiply(qxizs,qD)];
        case 2 % passive
            o2sym = [qmultiply(qC,qxizs), qmultiply(qD,qxizs)];
    end
end
