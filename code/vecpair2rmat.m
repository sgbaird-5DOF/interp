function R = vecpair2rmat(v1,v2)
arguments
   v1(1,3) double
   v2(1,3) double
end
% VECPAIR2RMAT  Compute a (non-unique) rotation matrix to go from v1 to v2.
%  If v1 == v2 or v1 == -v2 within the given numerical precision, then the
%  identity matrix or -1*identity matrix is given, respectively. Passive
%  rotation matrix
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-01
%
% Inputs:
%  v1, v2 - two vectors to find the rotation matrix to go from v1->v2
%
% Outputs:
%  R - Rotation matrix such that R*v1 == v2, and R\v2 == v1
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

% if epsijk == -1 %-1 seems to produce consistent results with section 4.2 of 10.1088/0965-0393/23/8/083501
% R = R.'; %change from active to passive rotation matrix
% end

end
%----------------------------vecpair2rmat----------------------------------
