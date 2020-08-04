function omega = get_omega_r2018a(o1,o2)
% arguments
% 	o1 (:,8) double {mustBeReal,mustBeFinite,normMustBeSqrt2(o1),quatNormsMustBeOne(o1)}
% 	o2 (:,8) double {mustBeReal,mustBeFinite,normMustBeSqrt2(o2),quatNormsMustBeOne(o2)}
% end
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

% adjust values that are close to 1 or -1
tol = 1e-6;
dot1(abs(dot1 - 1) < tol) = 1;
dot2(abs(dot2 - 1) < tol) = 1;

dot1(abs(dot1 + 1) < tol) = -1;
dot2(abs(dot2 + 1) < tol) = -1;

omega = 2*acos(abs(dot1+dot2)/2);

end
%----------------------------END get_omega()-------------------------------


%--------------------------VALIDATION FUNCTIONS----------------------------
% Custom validation function
function normMustBeSqrt2(val)

tol = 0.1;
normcheck = abs(norm(val(1,:))-sqrt(2)) < tol;
errmsg = ['octonion norm must be equal to sqrt(2) within tol = ' num2str(tol)];
assert(normcheck,errmsg)

end

function quatNormsMustBeOne(val)

tol = 0.1;
qAcheck = abs(norm(val(1,1:4))-1) < tol;
qBcheck = abs(norm(val(1,5:8))-1) < tol;
errmsg = ['quaternion norms must be equal to 1 within tol = ' num2str(tol)];
assert(qAcheck && qBcheck,errmsg)

end
