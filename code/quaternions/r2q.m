% r should be an npts-by-3 array

function q = r2q(r)

%---convert to rotation angle parameters---%
w = 2*atan(sqrt(sum(r.^2,2)));
[p,t,~] = cart2sphZ(r(:,1),r(:,2),r(:,3));

%---convert to quaternions---%
q = rot2q(w,t,p);

%---fix results for the identity---%
is_id = w == 0;
q(is_id,2:4) = zeros(sum(is_id),3);