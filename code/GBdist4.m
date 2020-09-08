function [dmin, o2minsyms] = GBdist4(o1,o2,pgnum,dtype,wtol,waitbarQ)
arguments
	o1(:,8) double {mustBeFinite,mustBeReal,mustBeSqrt2Norm}
	o2(:,8) double {mustBeFinite,mustBeReal,mustBeSqrt2Norm}
	pgnum(1,1) double {mustBeInteger} = 32 % default == cubic Oh point group
	dtype char {mustBeMember(dtype,{'omega','norm'})} = 'norm'
	wtol(1,1) double {mustBeFinite,mustBeReal} = 1e-6 %omega tolerance
	waitbarQ logical = false
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: modified version of GBdist function by CMU group. Keeps o1
% constant.
%
% Inputs:
%		o1, o2 -	octonions
%
%		pgnum - point group number
%
%		dtype - distance type ('omega' arc length or euclidean 'norm')
%
% Outputs:
%		dmin -	minimized distance metric
%
%		o2minsyms - minimized octonions
%
% Usage:
%
%		[dmin, o2minsyms] = GBdist4(o1,o2);
%
%		[dmin, o2minsyms] = GBdist4(o1,o2,32);
%
%		[dmin, o2minsyms] = GBdist4(o1,o2,32,'norm');
%
%		[dmin, o2minsyms] = GBdist4(o1,o2,32,'omega');
%
% Dependencies:
%		osymsets.m
%			--osymset.m
%				--qmult.m
%
%		get_omega.m
%
%		zeta_min2.m (naming distinct from 'zeta_min' to prevent conflicts in
%						GBdist.m)
%
% Notes:
%
%--------------------------------------------------------------------------
prec = 12; %precision
tol = 1e-6; %tolerance

%number of octonion pairs
npts = size(o1,1);

grainexchangeQ = true;
doublecoverQ = true;
%get symmetric octonions (SEOs)
% osets = osymsets(o2,pgnum,struct,grainexchangeQ,doublecoverQ); %out of memory 2020-08-03

%assign distance fn to handle
switch dtype
	case 'omega'
		distfn = @(o1,o2) get_omega(o1,o2);
	case 'norm'
		distfn = @(o1,o2) vecnorm(o1-o2,2,2);
end

dmin = zeros(1,npts);
o2minsyms = cell(1,npts);

%textwaitbar setup
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateProgress);
nsets = npts;
ninterval = 20;
N=nsets;
p=1;
reverseStr = '';
if nsets > ninterval
	nreps2 = floor(nsets/ninterval);
	nreps = nreps2;
else
	nreps2 = 1;
	nreps = nreps2;
end

function nUpdateProgress(~)
	percentDone = 100*p/N;
	msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
	fprintf([reverseStr, msg]);
	reverseStr = repmat(sprintf('\b'), 1, length(msg));
	p = p + nreps;
end

%loop through octonion pairs
parfor i = 1:npts %parfor compatible
	%text waitbar
	if mod(i,nreps2) == 0
		if waitbarQ
			send(D,i);
		end
	end
	
	%% setup	
	%unpack SEOs
	oset = osymsets(o2(i,:),pgnum,struct,grainexchangeQ,doublecoverQ);
	o2tmp = oset{1};
	% 	o2tmp = osets{i};
	
	%number of CSEOs
	nsets = size(o2tmp,1);
% 	nsets = size(osets{i},1);
	
	%unpack first octonion (held constant)
	o1tmp = o1(i,:);
	
	%copy octonion
	o1rep = repmat(o1tmp,nsets,1);
	
	%unpack quaternions
	qSC = o2tmp(:,1:4);
	qSD = o2tmp(:,5:8);
	
	%% apply U(1) symmetry
	% get minimum zeta & sigma values (zm)
	zm = zeta_min2(o1rep,o2tmp);
	qzm = [cos(zm/2) zeros(nsets,2) sin(zm/2)];
	
	% get minimized quaternions
	qCz = qmult(qSC,qzm);
	qDz = qmult(qSD,qzm);
	
	%package quaternions
	o2syms = [qCz qDz];
	
	%package quaternions
% 	o2syms = [...
% 		qCz	qDz
% 		-qCz	qDz
% 		qCz	-qDz
% 		-qCz	-qDz
% 		qDz	qCz
% 		-qDz	qCz
% 		qDz	-qCz
% 		-qDz	-qCz];
	
	%% compute distances
	%give the octonions a norm of sqrt(2)
	o1rep = sqrt2norm(o1rep,'oct');
	
% 	o1rep = repelem(o1rep,8,1);
	
	%compute all distances
	dlist = distfn(o1rep,o2syms); %#ok<PFBNS> %either omega or euclidean norm (see disttype arg)

	%% find minimum distances & octonions
	%get first instance of minimum omega
	dmin(i) = min(dlist);
	
	%find logical indices of all minimum omegas
	minIDs = ismembertol(dlist,dmin(i),wtol,'DataScale',1); %loosened tol for min omegas, 2020-07-28

	%find corresponding symmetrized octonions (with duplicates)
	o2minsymsTmp = o2syms(minIDs,:);
	
	%delete duplicate rows (low tol OK b.c. matching 8 numbers)
	[~,minIDs] = uniquetol(round(o2minsymsTmp,prec),tol,'ByRows',true,'DataScale',1); %generally repeats will be the same within ~12 sig figs
	o2minsyms{i} = o2minsymsTmp(minIDs,:);
	
end

end %GBdist4.m

%----------------------------CODE GRAVEYARD--------------------------------
%{
%from CMU group GBdist.m function
%now we implement U(1) and grain exchange symmetry
	
	%1. (A B C'(zeta) D'(zeta))
	zm1 = zeta_min(qA,qB,qSC,qSD);
	qzm1 = [cos(zm1/2) 0 0 sin(zm1/2)];
	qCz1 = qmult(qSC,qzm1);
	qDz1 = qmult(qSD,qzm1);
	
	w1 = norm([qA,qB]-[qCz1,qDz1]);
	w5 = norm([qA,qB]-[-qCz1,qDz1]);
	w9 = norm([qA,qB]-[qCz1,-qDz1]);
	w13 = norm([qA,qB]-[-qCz1,-qDz1]);
	

	
	sm1 = zeta_min(qB,qA,qSC,qSD);
	qsm1 = [cos(sm1/2) 0 0 sin(sm1/2)];
	qCs1 = qmult(qSC,qsm1);
	qDs1 = qmult(qSD,qsm1);
	
	w2 = norm([qA,qB]-[qCs1,qDs1]);
	w6 = norm([qA,qB]-[-qCs1,qDs1]);
	w10 = norm([qA,qB]-[qCs1,-qDs1]);
	w14 = norm([qA,qB]-[-qCs1,-qDs1]);
	
	%3. (A -B C'(zeta') D'(zeta'))
	
	zm2 = zeta_min(qA,-qB,qSC,qSD);
	qzm2 = [cos(zm2/2) 0 0 sin(zm2/2)];
	qCz2 = qmult(qSC,qzm2);
	qDz2 = qmult(qSD,qzm2);
	
	
	w3 = norm([qA,qB]-[qCz2,qDz2]);
	w7 = norm([qA,qB]-[-qCz2,qDz2]);
	w11 = norm([qA,qB]-[qCz2,-qDz2]);
	w15 = norm([qA,qB]-[-qCz2,-qDz2]);
	
	%4. (B -A C'(sigma') D'(sigma'))
	
	sm2 = zeta_min(qB,-qA,qSC,qSD);
	qsm2 = [cos(sm2/2) 0 0 sin(sm2/2)];
	qCs2 = qmult(qSC,qsm2);
	qDs2 = qmult(qSD,qsm2);
	
	w4 = norm([qA,qB]-[qCs2,qDs2]);
	w8 = norm([qA,qB]-[-qCs2,qDs2]);
	w12 = norm([qA,qB]-[qCs2,-qDs2]);
	w16 = norm([qA,qB]-[-qCs2,-qDs2]);
	
	%store candidate omega values
	wvec = [w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11 w12 w13 w14 w15 w16];
	wveclist{k}(ctrange) = wvec.';


	%unpack quaternions
% 	qCtmp = o2tmp(:,1:4);
% 	qDtmp = o2tmp(:,5:8);


	%unpack SEOs
	oset = osets{i};
	
	%unpack oset
	qSC = oset(:,1:4);
	qSD = oset(:,5:8);
	
	%% apply U(1) and grain exchange
	%1. (A B C'(zeta) D'(zeta))
	%2. (A B D'(sigma) C'(sigma))
	%3. (A B C'(zeta') -D'(zeta'))
	%4. (A B D'(sigma') -C'(sigma'))
	
	%get quaternion pairs for zeta_min()
	o2tmp = ...
		[qSC	qSD		%zeta1
		qSD	qSC		%sigma1
		qSC	-qSD		%zeta2
		qSD	-qSC];	%sigma2

% 	%copy quaternions
% 	qSC = repmat(qSC,4,1);
% 	qSD = repmat(qSD,4,1);


% 	%get quaternion pairs for zeta_min()
% 	o2tmp = ...
% 		[qSC	qSD		%zeta1
% 		qSD	qSC		%sigma1
% 		qSC	-qSD		%zeta2
% 		qSD	-qSC];	%sigma2
	
% 	o2tmp = [qSC,qSD];



	
	%% apply U(1) and grain exchange
	%1. (A B C'(zeta) D'(zeta))
	%2. (A B D'(sigma) C'(sigma))
	%3. (A B C'(zeta') -D'(zeta'))
	%4. (A B D'(sigma') -C'(sigma'))
	
	
	% 	%apply double cover
% 	o2syms = ...
% 		[qCz qDz
% 		-qCz qDz
% 		qCz -qDz
% 		-qCz -qDz]; %symmetrically equivalent candidate octonion pairs


% 	o1rep = repmat(o1rep,4,1); % expand to total rows == nsets*16

%}