function TF = inmisFZ(dlist,A,b,tol)
arguments
	dlist(:,3) double
	A = []
	b = []
	tol = 1e-3
end
% INMISFZ  check which rodrigues vectors fall inside the misorientation fundamental zone via vert2con (FEX)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-08-13
%
% Inputs:
%  A,b - constraint output from vert2con.m (FEX)
%
% Outputs:
%  TF - logical column vector, 1 if in misFZ (within tol), 0 if outside
%  misFZ (within tol)
%
% Usage:
%	[A,b] = misFZcon(); %pre-compute constraints
%  TF = inmisFZ(A,b,dlist,1e-3); %faster when A and b are specified if
%											%inmisFZ is called multiple times
%
% Dependencies:
%  misFZcon.m
%		--vert2con.m
%		--q2rod.m
%
% Notes:
%  *
%--------------------------------------------------------------------------
%compute constraints if not supplied
if isempty(A) || isempty(b)
	[A,b] = misFZcon();
end

% check dlist against constraints
TF = A*dlist.' <= b + tol;

% find dlist entries that meet all constraints (and transpose into column vector)
TF = all(TF,1).';

end