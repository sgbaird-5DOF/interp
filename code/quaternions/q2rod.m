% convert quaternion to rodriguez vector
function r = q2rod(q)

r = q(:,2:4)./q(:,1);