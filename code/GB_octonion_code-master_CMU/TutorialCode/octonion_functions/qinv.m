function out = qinv(q)
%vectorized by SGB 2020-08-15
%%% Lp(r) = vec(prp*) active. For passive, implement Lp*(r) = vec(p*rp)
% also see https://dx.doi.org/10.1088/0965-0393/23/8/083501 eq. 24 in
% regards to Lp(r)

q0 = q(:,1);
q = q(:,2:4);
out = [q0 -q];
% out = normr(out1);
% out = out1./vecnorm(out1,2);

% qr = pp(1)*qq(1)-dot(p,q);
% qi = pp(1)*q + qq(1)*p + cross(p,q);


end
