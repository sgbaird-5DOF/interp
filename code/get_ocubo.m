function o = get_ocubo(n,method,sidelength,seed)
arguments
	n {mustBeNonNegIntegerOrEmpty} = 1
	method char {mustBeMember(method,{'random','uniform'})} = 'random'
	sidelength {mustBeNonNegIntegerOrEmpty} = double.empty
	seed = 'shuffle'
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
%		o = get_ocubo(5) %generate 5 octonions from pairs of quaternions
%		randomly sampled from cubochoric space
%
%		o = get_ocubo(5,'random') %generate 5 octonions from pairs of
%		quaternions randomly sampled from cubochoric space
%
%		o = get_ocubo(100,'uniform') %generate 100 octonions randomly sampled
%		from list of pairs of quaternions generated via uniform cubochoric
%		sampling with automatically calculated sidelength (ceil(100^(1/3))
%
%		o = get_ocubo([],'uniform',5) %generate all combinations of
%		quaternion pairs (i.e. octonions) using 5^3 == 125 uniformly sampled
%		quaternions (15625 octonions)
%
% Dependencies:
%		allcomb.m (optional if nboQ == false)
%
%		ocubo.m
%			--cu2qu.m (De Graef group)
%
% Note: specifying 'random' and a non-empty sidelength will error, as these
% are two contradictory options.
%
%--------------------------------------------------------------------------
%set random number generator
rng(seed)

% argument validation (cont.)
if strcmp(method,'random') && ~isempty(sidelength)
	error('sidelength should not be specified for random sampling')
end
%--------------------------------------------------------------------------
%setup

if strcmp(method,'uniform') && isempty(sidelength)
	sidelength = ceil(n^(1/3)); % auto-calculate sidelength
	nq = [];
else
	if isempty(sidelength)
		nq = n;
	else
		nq = sidelength^3; %auto-calculate # of quaternions
	end
end
%--------------------------------------------------------------------------

switch method
	case 'random'
		% get 2*n cubochoric points
		q = get_cubo(2*nq,method,sidelength);
		
		%unpack
		qA = q(1:n,:);
		qB = q(n+1:2*n,:);
		
	case 'uniform'
		%get n cubochoric points
		q = get_cubo(nq,method,sidelength);
		
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
		
		if ~isempty(n)
			if (n < length(q)) && ~isempty(sidelength)
				% get a random list of the quaternions
				randlist = randi(size(qA,1),n,1);
				qA = qA(randlist,:);
				qB = qB(randlist,:);
			end
		end
		
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

