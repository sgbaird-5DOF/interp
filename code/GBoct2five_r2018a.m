function five = GBoct2five_r2018a(octlist,disQ)
% arguments
% 	octlist(:,8) double {mustBeNumeric,mustBeFinite}
% 	disQ(1,1) logical = true
% end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description: Do the inverse operation of GBfive2oct.m (CMU group,
% octonion code), and package into "five" structure
%
% Inputs:
%		octlist - list of octonions (rows)
%
%		opts - structure with fields parforQ and disQ
%
% Outputs:
%		five - struct with fields q, nA, d, and geometry (quaternion, BP
%		normal, rodrigues vector, and misFZ geometry type, respectively)
%
% Dependencies:
%		findgeometry.m
%
%--------------------------------------------------------------------------
npts = size(octlist,1);

%initialize
five(npts) = struct;
five(1).q = [];
five(1).nA = [];
five(1).d = [];
five(1).geometry = '';

%textwaitbar setup
waitbarQ = true;
slurmQ = true;
if waitbarQ
	D = parallel.pool.DataQueue;
	afterEach(D, @nUpdateProgress);
	N=npts;
	p=1;
	reverseStr = '';
	if slurmQ
		nintervals = 20;
	else
		nintervals = 100;
	end
	if npts > nintervals
		nreps = floor(npts/nintervals);
		nreps2 = floor(npts/nintervals);
	else
		nreps = 1;
		nreps2 = 1;
	end
else
	D = [];
	nreps2 = 0;
end
%convert subdivided points to 5DOF
disp(' ')
disp('GBoct2five ')
for i = 1:npts %parfor compatible
	%textwaitbar
	if waitbarQ
		if mod(i,nreps2) == 0
			send(D,i);
		end
	end
	
	%unpack octonion
	oct = octlist(i,:);
	
	%get quaternion and BP normal
	[q,nA] = GBoct2five_once(oct);
	
	%convert
	d = q2rod(q);
	
	%package
	five(i).q = q;
	five(i).nA = nA;
	five(i).d = d;	
end
if waitbarQ
	disp(' ')
end

%call to disorientation might be expensive
if disQ
	geometry = findgeometry(disorientation(vertcat(five.q),'cubic'));
else
	geometry = findgeometry(vertcat(five.q));
end

[five.geometry] = geometry{:};

	function nUpdateProgress(~)
		percentDone = 100*p/N;
		msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
		fprintf([reverseStr, msg]);
		reverseStr = repmat(sprintf('\b'), 1, length(msg));
		p = p + nreps;
	end

end

function [qm,nA] = GBoct2five_once(o)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
%
% Inputs:
%
% Outputs:
%
% Dependencies:
%
%--------------------------------------------------------------------------
qA = normr(o(1:4)); %has to be normalized, otherwise it's doesn't go back to same 5DOF parameters
qB = normr(o(5:8));

qm = qmult(qB,qinv(qA));

pA = normr(qmult(qm,qA)); %normr added 2020-07-16 to avoid NaN values in qu2ax
ax = qu2ax(pA);

axisA = ax(1:3);
phiA = ax(4);

mA0(1) = 0;
scl = sqrt(1-cos(phiA)^2);
mA0(2) = -axisA(2)*scl;
mA0(3) = axisA(1)*scl;
mA0(4) = cos(phiA);

nA = qmult(qinv(qm),qmult(mA0,qm));
nA = normr(nA(2:4));

end

%------------------------CODE GRAVEYARD------------------------------------
%{

%*norm([axisA(2),axisA(1)])
 %*norm([axisA(2),axisA(1)])

% qm*.mA0.qm = qm*.qm.[0 nA].qm*.qm
%
% nA = qm*.mA0.qm

% nA = [0 0 1];
% qA = o(1:4);
% z = [0 0 0 1];
% nA = qmult(-qm,qmult(z,qinv(qm)));
% nA2 = qmult(-qB,qmult(z,qinv(qB)));
%
% qmult(qinv(qA),qmult(z,qA))
%
% qmult((qmis),qmult([0 nA],qinv(qmis)));
%


% qs = disorientation([qA;qB],'cubic');
% qA = qs(1,:);
% qB =qs(2,:);

% qA = o(1:4);
% qB = o(5:8);

% if nargin == 2
% 	parforQ = varargin{1};
% else
% 	parforQ = true;
% end


% if npts == 1
% 	[five.q,five.nA] = GBoct2five_once(octlist);
% 	five.d = q2rod(five.q);
% 	geometry = findgeometry(disorientation(q,'cubic'));
% 	five(i).geometry = geometry;
% 	return
% end

		%call to disorientation might be expensive
		if disQ
			geometry = findgeometry(disorientation(q,'cubic'));
		else
			geometry = findgeometry(q); %if outside FZ, still registers as 'interior' with current implementation
		end


defaultparforQ = true;
defaultdisQ = true;

P = inputParser;
addRequired(P,'octlist',@isnumeric);
addParameter(P,'parforQ',defaultparforQ,@islogical);
addParameter(P,'disQ',defaultdisQ,@islogical);
parse(P,octlist,varargin{:});

parforQ = P.Results.parforQ;
disQ = P.Results.disQ;

if size(octlist,2) == 7
	octlist = [octlist zeros(npts,1)];
end


	if ~disQ
		five(i).geometry = findgeometry(q);
	end

%}
