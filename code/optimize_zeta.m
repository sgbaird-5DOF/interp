function [xopt, fopt, exitflag, output] = optimize_zeta(o,dtype,z0)
arguments
	o(:,8) double {mustBeReal,mustBeFinite}
	dtype char {mustBeMember(dtype,{'omega','norm'})} = 'omega'
	z0 double {mustBeReal,mustBeFinite} = double.empty
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: Minimize the pairwise distance matrix of a set of GBs w.r.t.
% zeta (twist angle via U(1) symmetry)
%
% Inputs:
%	o - list of input octonions (pass into get_octpairs.m first)
%
% Outputs:
%	b - b
%
% Usage:
%	c = b(a);
%
% Dependencies:
%	*
%
% Notes:
%	*
%--------------------------------------------------------------------------

npts = size(o,1);

% starting point
if isempty(z0)
	z0 = zeros(npts,1);
end

% bounds
ub = 2*pi*ones(npts,1);
lb = zeros(npts,1);

% linear constraints
A = [];
b = [];
Aeq = [];
beq = [];
nlc = [];

IDs = allcomb(1:npts,1:npts);

o1 = o(IDs(:,1),:);
o2 = o(IDs(:,2),:);

% "true" pairwise distance metrics
pdtrue = GBdist4(o1,o2,32,dtype,1e-6,true).';

assert(size(pdtrue,1) == size(o1,1),'pdtrue is wrong size. Maybe transpose?')

%objective function
obj = @(z) pd_sse(o,z,dtype,pdtrue);

%minimization
opttype = 'fmincon';
switch opttype
    case 'fmincon'
        opts = optimoptions('fmincon','Display','iter');
        [xopt, fopt, exitflag, output] = fmincon(obj,z0,A,b,Aeq,beq,lb,ub,nlc,opts);
    case 'ga'
        opts = optimoptions('ga','Display','iter');
        ga(obj,npts)
end

end

%------------------------CODE GRAVEYARD------------------------------------
%{

% ------------Starting point and bounds------------
% var =
npts = size(o,1);
x0 = zeros(1,npts);
ub = 2*pi*ones(npts,1);
lb = zeros(npts,1);

% ------------Linear constraints------------
A = [];
b = [];
Aeq = [];
beq = [];

% ------------Objective and Non-linear Constraints------------
function [f, c, ceq] = objcon(x)
		
		% design variables
		z = x;
		
		% other analysis variables
		qC = o(:,1:4);
		qD = o(:,5:8);
		
		%analysis functions
		qz = [cos(z/2) zeros(n,2) sin(z/2)];
		qCz = qmult(qlist,qz);
		qDz = qmult(qD,qz);
		o_out = [qCz qDz];
		
		pd = zeros(npts);
		
		tmp = num2cell(octvtx,2);
		opairs = allcomb(tmp,tmp);
		o1 = vertcat(opairs{:,1});
		o2 = vertcat(opairs{:,2}); %if I want to use the vectorized versions, I'll need to deal with preserving PD-formating
		
		% objective function
		F =
		
		%inequality constraints (c<=0)
		c = zeros(3,1);         % create column vector
		c(1) = maxwidth-0.75;
		
		%equality constraints (ceq=0)
		ceq = [];
	end

% ------------Call fmincon------------
options = optimoptions(@fmincon, 'display', 'iter-detailed');
[xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);


% ------------Separate obj/con (do not change)------------
function [f] = obj(x)
		[f, ~, ~] = objcon(x);
	end
function [c, ceq] = con(x)
		[~, c, ceq] = objcon(x);
	end
end


tmp = num2cell(o,2);
opairs = allcomb(tmp,tmp);
o1 = vertcat(opairs{:,1});
o2 = vertcat(opairs{:,2}); %if I want to use the vectorized versions, I'll need to deal with preserving PD-formating


tmp = num2cell(o_out,2);
opairs = allcomb(tmp,tmp);
o1 = vertcat(opairs{:,1});
o2 = vertcat(opairs{:,2}); %if I want to use the vectorized versions, I'll need to deal with preserving PD-formating


pd = zeros(npts^2,1);

if isempty(pdtrue)
	IDs = allcomb(1:npts,1:npts);
	o1 = o(IDs(:,1),:);
	o2 = o(IDs(:,2),:);
	% "true" pairwise distance metrics
	pdtrue = GBdist4(o1,o2,32,dtype);
end

%}