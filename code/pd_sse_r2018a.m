function err = pd_sse_r2018a(o,z,dtype,pdtrue,pdtruefn,errtype)
% arguments
% 	o double
% 	z(:,1) double
% 	dtype char {mustBeMember(dtype,{'omega','norm'})} = 'omega'
% 	pdtrue(:,1) double {mustBeFinite,mustBeReal} = double.empty
% 	pdtruefn function_handle = @(o1,o2,dtype) GBdist4(o1,o2,32,dtype)
% 	errtype char {mustBeMember(errtype,{'rmse','mse','me'})} = 'me'
% end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: get the error of a pairwise distance matrix relative to the
% true pairwise distance matrix
%
% Inputs:
%	a - a
%
% Outputs:
%	b - b
%
% Usage:
%	c = b(a);
%
% Dependencies:
%  qmult.m
%
%  allcomb.m
%
% Notes:
%	It's generally more efficient to supply pdtrue rather than calculate it
%	each time, especially if the function call is expensive, as in GBdist
%	and counterparts.
%--------------------------------------------------------------------------
npts = size(o,1);

assert(any(size(z) == npts),'z is incorrect size.')
if size(z,2) == npts
    z = z.';
end

%calculate pdtrue if not supplied
if isempty(pdtrue)
	IDs = allcomb(1:npts,1:npts);
	o1 = o(IDs(:,1),:);
	o2 = o(IDs(:,2),:);
	% "true" pairwise distances
	pdtrue = pdtruefn(o1,o2,dtype);
end

%unpack quaternions
qC = o(:,1:4);
qD = o(:,5:8);

%get and apply rotations
qz = [cos(z/2) zeros(npts,2) sin(z/2)];
qCz = qmult(qC,qz);
qDz = qmult(qD,qz);
onew = [qCz qDz];

%new combinations for pd dists
IDs = allcomb(1:npts,1:npts);
o1 = onew(IDs(:,1),:);
o2 = onew(IDs(:,2),:);

%pairwise distances
switch dtype
	case 'omega'
		pd = get_omega_r2018a(o1,o2);
	case 'norm'
		pd = vecnorm(o1-o2,2,2);
end

%mean error
switch errtype
	case 'rmse'
		err = sqrt(immse(pd,pdtrue));
	case 'mse'
        err = immse(pd,pdtrue);
	case 'me'
		err = mean(pd-pdtrue);
    case 'se'
        err = sum(pd-pdtrue);
end

end

%------------------------CODE GRAVEYARD------------------------------------
%{
% 		err = mean((pd-pdtrue).^2);
%}
