function misFZ_ids = InCubicFZ(R,tol)
arguments
	R(:,3) double
	tol(1,1) double = 1e-3
end
% INCUBICFZ  Test whether or not a set of rodrigues vectors is within the cubic misorientation fundamental zone.
%--------------------------------------------------------------------------
% Authors: Eric Homer, Sterling Baird
%
% Date: 2020-07-20
% 
% Inputs: R - rows of rodrigues vector
%
% Outputs: inmisFZ - logical array, whether or not InCubicFZ
%
% Dependencies:
%
% References:
%		Orientation mapping [1] F.C. Frank, Orientation mapping, Metall Trans
%		A.19 (1988) 403–408. doi:10.1007/BF02649253.
%
% Notes:
% 
% function received from Eric Homer, modified for matrix input
%--------------------------------------------------------------------------

assert(size(R,2) == 3,['input must be a rodrigues vector of 3 elements, not ' int2str(size(R,2))])

R1 = R(:,1);
R2 = R(:,2);
R3 = R(:,3);

negcheck = (R3 < -tol) | (R2 < R3-tol) | (R1 < R2-tol);

sqrtcheck = R1 > (sqrt(2)-1) + tol;

distcheck = (R1+R2+R3) > 1 + tol;

misFZ_ids = ~(negcheck | sqrtcheck | distcheck);

end


%---------------------------CODE GRAVEYARD---------------------------------
%{
Original code:
%function received from Eric Homer

TF=true;

if R(3) >= 0 && R(2) >= R(3) && R(1) >= R(2)
	%do nothing
else
	TF=false;
	return;
end

if R(1) > sqrt(2)-1
	TF=false;
	return;
end

if R(1) + R(2) + R(3) > 1
	TF=false;
	return;
end

end
%}





