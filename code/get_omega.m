function omega = get_omega(o1,o2)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
% 
% Inputs: o1, o2 form an octonion pair, each can be rows of octonions
%
% Outputs:
%
% Dependencies:
%
%--------------------------------------------------------------------------
qA = o1(:,1:4);
qB = o1(:,5:8);
qC = o2(:,1:4);
qD = o2(:,5:8);

dot1 = dot(qA,qC,2);
dot2 = dot(qB,qD,2);

omega = 2*acos(abs(dot1+dot2)/2);

end
%----------------------------END get_omega()-------------------------------
