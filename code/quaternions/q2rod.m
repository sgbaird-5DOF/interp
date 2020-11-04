function r = q2rod(q)
% Q2ROD  convert quaternion to rodriguez vector
r = q(:,2:4)./q(:,1);