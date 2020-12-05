function [A,b] = misFZcon(epsijk)
arguments
   epsijk(1,1) double = 1 
end
% MISFZCON  get constraints for misorientation fundamental zone (misFZ)
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-08-13
%
% Inputs:
%  N/A
%
% Outputs:
%  A,b - outputs from vert2con.m given vertices of misFZ
%
% Usage:
%  [A,b] = misFZcon();
%
% Dependencies:
%  q2rod.m
%	vert2con.m
%
% Notes:
%  (1) Patala, S.; Schuh, C. A. Symmetries in the Representation of Grain
%  Boundary-Plane Distributions. Philosophical Magazine 2013, 93 (5),
%  524–573. https://doi.org/10.1080/14786435.2012.722700.

%--------------------------------------------------------------------------
%define quaternions (corresponding to vertices of misFZ in rodrigues space), mostly from (1)
k = sqrt(2)-1;

qA = [cos(pi/8),sin(pi/8),0,0];

qB = [1/sqrt(1+2*k^2),k/sqrt(1+2*k^2),k/sqrt(1+2*k^2),0];

qC = 1/(2*sqrt(2))*[1/k,1,1,k];

qD = [(7+(-4).*2.^(1/2)).^(-1/2),...
	((1/17).*(5+(-2).*2.^(1/2))).^(1/2),...
	((1/34).*(5+(-2).*2.^(1/2))).^(1/2),...
	((1/34).*(5+(-2).*2.^(1/2))).^(1/2)];

qE = [sqrt(3)/2,1/(2*sqrt(3)),1/(2*sqrt(3)),1/(2*sqrt(3))];

qO = [1 0 0 0];

%catenate quaternions
qlist = [qA;qB;qC;qD;qE;qO];

%convert to rodrigues vectors
% ro = qu2ro(qlist,epsijk);
% dlist = ro(:,4).*ro(:,1:3);

dlist = q2rod(qlist);

%compute constraints
[A,b] = vert2con(dlist); %A*x <= b

end
