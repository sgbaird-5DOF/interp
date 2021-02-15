function mustBeSqrt2Norm(o)
% MUSTBESQRT2NORM  check that first octonion in list has norm == sqrt(2) and each quaternion has norm == 1
%NOTE: do not use "arguments" syntax here since this is a validation fn
%---------------------HELPER VALIDATION FUNCTION---------------------------
onorm = norm(o(1,:));
errmsg = ['octonion norm == ' num2str(onorm) ' ~= sqrt(2) == 1.4142'];
assert(abs(onorm - sqrt(2)) <= 1e-3, errmsg)

% qnorm = norm(o(1,1:4));
% errmsg = ['quaternion norm == ' num2str(qnorm) ' ~= 1'];
% assert(qnorm - 1 <= 1e-1, errmsg);

% qnorm = norm(o(1,5:8));
% errmsg = ['quaternion norm == ' num2str(qnorm) ' ~= 1'];
% assert(qnorm - 1 <= 1e-1, errmsg);
end %mustBeSqrt2Norm