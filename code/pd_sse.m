function [err,onew,pd,pdtrue] = pd_sse(o,z,dtype,pdtrue,pdtruefn,errtype)
arguments
	o double
	z(:,1) double
	dtype char {mustBeMember(dtype,{'omega','norm'})} = 'omega'
	pdtrue(:,1) double {mustBeFinite,mustBeReal} = double.empty
	pdtruefn function_handle = @(o1,o2,dtype) GBdist4(o1,o2,32,dtype,1e-6,true)
	errtype char {mustBeMember(errtype,{'rmse','mse','me'})} = 'me'
end
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
%  get_errmetrics.m
%
% see also PDIST
%
% Notes:
%	It's generally more efficient to supply pdtrue rather than calculate it
%	each time, especially if the function call is expensive, as in GBdist
%	and counterparts.
%--------------------------------------------------------------------------
npts = size(o,1);

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
		pd = get_omega(o1,o2);
	case 'norm'
		pd = vecnorm(o1-o2,2,2);
end

%get error
err = get_errmetrics(pd,pdtrue,errtype);
end

%---------------------------------CODE GRAVEYARD---------------------------
%{
% switch errtype
% 	case 'rmse'
% 		err = sqrt(immse(pd,pdtrue));
% 	case 'mse'
% 		err = immse(pd,pdtrue);
% % 		err = mean((pd-pdtrue).^2)
% 	case 'me'
% 		err = mean(pd-pdtrue);
% end

%}
