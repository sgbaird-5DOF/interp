% convert rodriguez vector to quaternion
function q = rod2q(r)

q(:,1) = sqrt(1./(1+sum(r.^2,2)));
q(:,2) = q(:,1).*r(:,1);
q(:,3) = q(:,1).*r(:,2);
q(:,4) = q(:,1).*r(:,3);