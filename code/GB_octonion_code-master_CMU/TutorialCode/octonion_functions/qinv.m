function out = qinv(q)
%%% Lp(r) = vec(prp*) active. For passive, implement Lp*(r) = vec(p*rp)

q0 = q(1);
q = q(2:4);
out1 = [q0 -q];
out = out1/norm(out1);

% qr = pp(1)*qq(1)-dot(p,q);
% qi = pp(1)*q + qq(1)*p + cross(p,q);


end
