function omega = get_omega(o1,o2)
arguments
	o1 (:,8) double {mustBeReal,mustBeFinite}
	o2 (:,8) double {mustBeReal,mustBeFinite}
end
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

% enables use as custom distance function to e.g. pdist2
if size(o1,1) == 1
    npts = size(o2,1);
    o1 = repmat(o1,npts,1);
end

qA = normr(o1(:,1:4));
qB = normr(o1(:,5:8));
qC = normr(o2(:,1:4));
qD = normr(o2(:,5:8));

dot1 = dot(qA,qC,2);
dot2 = dot(qB,qD,2);

omega = real(2*acos(abs(dot1+dot2)/2)); %added real() 2020-08-03

end
%----------------------------END get_omega()-------------------------------


% %--------------------------VALIDATION FUNCTIONS----------------------------
% % Custom validation function
% function normMustBeSqrt2(val)
% 
% tol = 0.1;
% normcheck = abs(norm(val(1,:))-sqrt(2)) < tol;
% errmsg = ['octonion norm must be equal to sqrt(2) within tol = ' num2str(tol)];
% assert(normcheck,errmsg)
% 
% end

% function quatNormsMustBeOne(val)
% 
% tol = 0.1;
% qAcheck = abs(norm(val(1,1:4))-1) < tol;
% qBcheck = abs(norm(val(1,5:8))-1) < tol;
% errmsg = ['quaternion norms must be equal to 1 within tol = ' num2str(tol)];
% assert(qAcheck && qBcheck,errmsg)
% 
% end


%----------------------------CODE GRAVEYARD--------------------------------
%{
% adjust values that are close to 1 or -1
tol = 1e-6;
dot1(abs(dot1 - 1) < tol) = 1;
dot2(abs(dot2 - 1) < tol) = 1;
dot1(dot1 > 1) = 1;
dot2(dot2 > 1) = 1;

dot1(dot1 < -1) = -1;
dot2(dot2 < -1) = -1;

dot1(abs(dot1 + 1) < tol) = -1;
dot2(abs(dot2 + 1) < tol) = -1;

,normMustBeSqrt2(o1)
,normMustBeSqrt2(o2)

%}
