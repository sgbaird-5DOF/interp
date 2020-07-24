function R = slightRot(d,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-17
%
% Description: rotate points slightly in every dimension
% 
% Inputs: 
%		d			=== dimensionality
%
%		alpha		=== rotation angle in degrees (optional)
%
% Outputs: R === compound rotation matrix
%
% Dependencies: RotMatrix.m (FEX)
%
%--------------------------------------------------------------------------
usual = 1;
if nargin - usual == 1
	alpha = deg2rad(varargin{1});
else
	alpha = deg2rad(0.01);
end

bases = 1-eye(d);
R = eye(d); %initial rotation matrix (i.e. identity matrix)
for i = 1:d-1
	basis1 = bases(i,:);
	basis2 = bases(i+1,:);
	R = R*RotMatrix(alpha,basis1,basis2); %compounded rot matrix
end

end %slightRot