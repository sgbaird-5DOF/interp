function out = qmult(pp,qq)

global epsijk
if isempty(epsijk)
	error('need to specify epsijk as global variable, e.g. via setGlobal_epsijk.m')
end
% multiply lists of quaternion pairs (input: rows of quaternions),
% vectorized by SGB 2020-07-27
p = pp(:,2:4); q = qq(:,2:4);
qr = pp(:,1).*qq(:,1) - dot(p,q,2);
qi = pp(:,1).*q + qq(:,1).*p + epsijk*cross(p,q,2);

out = [qr qi];

end

%----------------------------CODE GRAVEYARD--------------------------------
%{
CMU group original function (single pair multiplication)
%%% multiply two quaternions p*q
% p = pp(2:4); q = qq(2:4);
% 
% qr = pp(1)*qq(1)-dot(p,q);
% qi = pp(1)*q + qq(1)*p + cross(p,q);
% 
% out = [qr qi];
%}
