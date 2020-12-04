function out = qmult(pp,qq,epsijk)
arguments
   pp(:,4) double {mustBeReal,mustBeFinite}
   qq(:,4) double {mustBeReal,mustBeFinite}
   epsijk = 1
end
% epsijk = 1;
% global epsijk
% if isempty(epsijk)
%     setGlobal_epsijk(1)
% end

% multiply lists of quaternion pairs (input: rows of quaternions),
% vectorized by SGB 2020-07-27
p = pp(:,2:4); q = qq(:,2:4);
qr = pp(:,1).*qq(:,1) - dot(p,q,2);
qi = pp(:,1).*q + qq(:,1).*p + epsijk*cross(p,q,2);

out = [qr qi];

% %negate quaternions if scalar part is negative
% ids = qr < 0;
% out(ids,:) = -out(ids,:);

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
