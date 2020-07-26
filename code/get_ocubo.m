function o = get_ocubo(n,method,sidelength)
arguments
	n {mustBeNonNegIntegerOrEmpty} = 1
	method string {mustBeMember(method,{'random','uniform'})} = 'random'
	sidelength {mustBeNonNegIntegerOrEmpty} = double.empty
end

if strcmp(method,'uniform') && isempty(sidelength)
	sidelength = ceil(n^(1/3)); % auto-calculate sidelength
	n = [];
elseif isempty(n)
	n = sidelength^3;
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-25
%
% Description: get octonions formed by pairs of quaternions from randomly
% or uniformly sampled cubochoric points. In general, for random, no two
% quaternions will be the same.
%
% Inputs:
%		n - # of octonions to output (re-calculated if using 'uniform' method
%		and sidelength is specified
%
%		method - sampling method, 'random', 'uniform'
%
%		sidelength - # of points along edge of cube used in cubochoric
%		sampling (automatically calculated if n is given and sidelength is
%		not specified)
%
% Outputs:
%
%		o - list of octonions
%
% Usage:
%		o = get_ocubo(); %generate a single octonion formed by two
%		quaternions sampled randomly from cubochoric space
% 
%		o = get_ocubo(5,'random') %generate 5 random octonions
%
%		o = get_ocubo(100,'uniform') %generate ~100 uniform octonions with
%		automatically calculated sidelength (ceil(100^(1/3))
%
%		o = get_ocubo(~,'uniform',5) %generate all combinations of quaternion
%		pairs (i.e. octonions) using 5^3 uniformly sampled quaternions (7750
%		octonions)
%
%		o = get_ocubo(100,'random',5) %generate 100 random quaternion pairs
%		(i.e. octonions) sampled from 5^3 uniformly sampled quaternions
%
% Dependencies:
%
%--------------------------------------------------------------------------

switch method
	case 'random'
		% get 2*n cubochoric points
		q = get_cubo(2*n,method,sidelength);
		
		%unpack
		qA = q(1:n,:);
		qB = q(n+1:2*n,:);
		
	case 'uniform'
		%get n cubochoric points
		q = get_cubo(n,method,sidelength);
		
		%convert to cell array of quaternions
		q = num2cell(q,2);
		
		%form pairs
		nboQ = true; %whether to include no-boundary octonions
		if nboQ
			qpairs = allcomb(q,q);
		else
			qpairs = nchoosek(q,2);
		end
		
		% unpack
		qA = vertcat(qpairs{:,1});
		qB = vertcat(qpairs{:,2});		
end

%catenate
o = [qA qB];

end %get_ocube.m

%----------------------CUSTOM VALIDATION FUNCTIONS-------------------------
function mustBeNonNegIntegerOrEmpty(arg)
errmsg = 'must be non-neg integer or empty []';
if ~isempty(arg)
	if floor(arg) ~= arg
		error(errmsg)
	elseif arg < 0
		error(errmsg)
	end
end
end

