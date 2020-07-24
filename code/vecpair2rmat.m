function R = vecpair2rmat(v1,v2)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Description: Compute rotation matrix. If v1 == v2 or v1 == -v2 within the
% given numerical precision, then the identity matrix or -1*identity matrix
% is given, respectively.
%
% Inputs:
%		v1, v2	===	two vectors to find the rotation matrix to go from
%							v1->v2
%
% Outputs:
%		R			=== Rotation matrix such that R*v1 == v2, and R\v2 == v1
%
% References
%	https://math.stackexchange.com/a/897677
%	https://math.stackexchange.com/a/476311/76513

precision = 12;
r = @(a,b) round(a-b,precision);

isEqual = all(r(v1,v2)==0); %pointing in same direction
isOpposite = all(r(v1,-v2)==0); %pointing in opposite direction
isParallel = isEqual || isOpposite; %either of previous two cases

if ~isParallel
	ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
	R = eye(3) + ssc(cross(v1,v2)) + ssc(cross(v1,v2))^2*(1-dot(v1,v2))/(norm(cross(v1,v2))^2);
	
elseif isEqual
	% same vector
	R = eye(3);
	
elseif isOpposite
	% vectors pointing in opposite directions
	R = -eye(3);
end

end
%----------------------------vecpair2rmat----------------------------------
