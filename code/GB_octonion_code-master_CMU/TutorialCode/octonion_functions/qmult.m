function out = qmult(pp,qq)
%%% multiply two quaternions p*q
p = pp(2:4); q = qq(2:4);

qr = pp(1)*qq(1)-dot(p,q);
qi = pp(1)*q + qq(1)*p + cross(p,q);

out = [qr qi];

end
