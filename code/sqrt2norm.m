function pts = sqrt2norm(pts,type)
% SQRT2NORM  take a set of octonions and give each row norm == sqrt(2) if (norm == 1) || (norm == sqrt(2))
arguments
	pts(:,8) double {mustBeNumeric,mustBeFinite}
	type char {mustBeMember(type,{'oct','quat'})} = 'oct'
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: a
% 
% Inputs:
%		a -	a
%
% Outputs:
%		b -	b
%
% Usage:
%		a = b(a);
%
% Dependencies:
%		*
%
% Notes:
%		*
%--------------------------------------------------------------------------
tol = 1e-1;
switch type
	case 'quat'
		% if norm(o) == 1 within tolerance, multiply by sqrt(2)
		if all(abs(vecnorm(pts(:,1:4),2,2)-1/sqrt(2)) < tol) ...
				&& ...
				all(abs(vecnorm(pts(:,5:8),2,2)-1/sqrt(2)) < tol)
			
			pts(:,1:4) = normr(pts(:,1:4));
			pts(:,5:8) = normr(pts(:,5:8));
			
		elseif all(abs(vecnorm(pts(:,1:4),2,2)-1) > tol) ...
				&& ...
				all(abs(vecnorm(pts(:,5:8),2,2)-1) > tol)
			
			error(['norm(qA) == ' num2str(norm(pts(1,1:4))) ', ' ...
				'norm(qB) == ' num2str(norm(pts(1,5:8))) ...
				'. norm of 1+ quaternions ~= 0.7071 || 1'])
			
		end
	case 'oct'
        ptnm = vecnorm(pts,2,2);
        if ~isempty(find(ptnm == 0, 1))
%            warning('identity octonion(s) present')
           ptnm(ptnm == 0) = [];
        end
		if all(abs(ptnm - 1) < tol)
			pts = normr(pts)*sqrt(2);
		elseif any(abs(ptnm - sqrt(2)) > tol)
			error('norm of octonions ~= 1 || sqrt(2)')
		end
end

end
