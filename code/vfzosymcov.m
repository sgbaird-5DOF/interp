function [cov,o,oref] = vfzosymcov(o,pgnum,covtype,covpars,grainexchangeQ,doublecoverQ,epsijk)
arguments
   o(:,8) = get_ocubo(100)
   pgnum(1,1) double = 32
   covtype char = 'squaredexponential'
   covpars double = [deg2rad(5),0.01]
   grainexchangeQ(1,1) logical = false
   doublecoverQ(1,1) logical = false
   epsijk(1,1) double = 1
%    NV.oref(1,8) double = get_ocubo(1,'random',[],10)
end
% VFZOSYMCOV  get symmetrized covariance matrix for an octonion set
%--------------------------------------------------------------------------
% Inputs:
%  o - list of octonions (these could be VFZOs or not)
%
%  pgnum - point group number
%
%  covtype - covariance type, e.g. 'squaredexponential'
%
%  covpars - parameters for kernel covariance, e.g. for 'squaredexponential':
%   covpars = [L, sigma];, where L === smoothness length, sigma === marginal
%   uncertainty value. Follows documentation in 'KernelParameters' section
%   of documentation for fitrgp, but check implemented kernels in d2cov.m
%
%  grainexchangeQ, doublecoverQ - whether or not to consider these types of
%  symmetries
%
%  epsijk - rotation convention (1 or -1)
%
% Outputs:
%  C - symmetrized covariance matrix
%
%  o - VFZ-symmetrized octonions
%
%  oref - reference octonion used to symmetrize octonions into VFZ
%
% Usage:
%  C = vfzosymset(o); %use all defaults
%  C = vfzosymset(o,L,sigma,pgnum,covtype,grainexchangeQ,doublecoverQ,epsijk);
%
% Dependencies:
%  *
%
% Notes:
%  See "Explanation of Covariance Matrix Averaging" section at end of file,
%  before "CODE GRAVEYARD" section
%
%  The VFZ-symmetrized octonion output should be identical within numerical
%  error to the original octonion set if the input is a VFZO set and the
%  same 'oref' is used to symmetrize.
%
% Author(s): Sterling Baird
%
% Date: 2020-12-12
%--------------------------------------------------------------------------

[o,oref] = get_octpairs(o,[],epsijk);

%% get symmetrically equivalent octonions (SEOs)
npts = size(o,1);
osets = osymsets(o,pgnum,[],grainexchangeQ,doublecoverQ,epsijk);
nsyms = cellfun(@(oset) size(oset,1),osets);
nsym = nsyms(1);
assert(all(nsyms-nsym==0),'all entries of osets should have same # of pts')

%% compile symmetrically equivalent octonions (SEOs)
osets = vertcat(osets{:});
osets = normr(osets); %for use with get_alen

%% U(1) symmetry
% setup
nptstot = npts*nsym;
orefrep = repmat(oref,nptstot,1);
zm = zeta_min2(orefrep,osets,-epsijk);
mA = [0 0 1]; %octonion convention that BP normal is [0 0 1] in lab frame
mArep = repmat(mA,nptstot,1);
qzm = ax2qu([mArep zm],-epsijk);
% apply
osets(:,1:4) = qmult(qzm,osets(:,1:4),epsijk);
osets(:,5:8) = qmult(qzm,osets(:,5:8),epsijk);

%% reshape the osets array
% IDs for interleaving vectors
ids = cell(nsym,1);
for i = 1:nsym
    ids{i} = i:nsym:nptstot;
end
intIDs = [ids{:}];
osets(intIDs,:) = osets;

%% textwaitbar setup
D = parallel.pool.DataQueue;
afterEach(D, @nUpdateProgress);
nsets = npts; % # parfor loop iterations
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
waitbarQ = true;

ofull = osets;
osets = mat2cell(osets,repelem(nsym,npts));


%% main loop
PDvertavgs = cell(1,npts);
parfor i = 1:npts
    %text waitbar
    if mod(i,nreps2) == 0
        if waitbarQ
            send(D,i);
        end
    end
    %unpack
    oset2 = osets{i};
    %compute distances as set of "vertically" stacked covariance matrices
    PDvert = pdist2(ofull,oset2);
    
    %split into individual arrays and average
    PDvertsplit = mat2cell(PDvert,repelem(npts,nsym));
    PDvertsplit = cat(3,PDvertsplit{:});
    PDvertavgs{i} = mean(PDvertsplit,[2 3]);
end
%average along the "horizontal" stack of covariance matrices
pd = [PDvertavgs{:}];

cov = d2cov(pd,covpars,covtype);
end

%% Explanation of Covariance Matrix Averaging
%{
N == npts

Take a "master" covariance matrix given by:

| C11 C12 C13 ... C1N |
| C21 C22 C23 ... C2N |
| C31 C32 C33 ... C3N |
|  .   .   .       .  |
| CN1 CN2 CN3 ... CNN |

o1 is the first SEO of the first GB (just the identity rotation, i.e. o(1,:) )
o1' is the second SEO of the first GB
o1['...] represents the final SEO of the first GB
C11 expanded, expressed as pairwise distance matrix:
            o1 o1' o1'' ... o1['...]
         ---------------------------
o1       |                         |
o1'      |                         |
o1''     |                         |
 .       |                         |
o1['...] |_________________________|

C21 expanded, expressed as pairwise distance matrix:
            o1 o1' o1'' ... o1['...]
         ---------------------------
o2       |                         |
o2'      |                         |
o2''     |                         |
 .       |                         |
o2['...] |_________________________|

etc.

"Symmetrized" covariance matrix is given by:
C(1,1) = mean(C11,'all')
C(2,1) = mean(C21,'all')
C(2,2) = mean(C22,'all')
...
C(N,N) = mean(CNN,'all')

Unfortunately, the master covariance matrix is generally too large to
compute at once, and rather slow to compute with a nested for loop.

So instead we calculate averages within vertical chunks of the matrix, where
the first column chunk looks like:
| C11 |
| C21 |
| C31 |
|  .  |
| CN1 |

and C11, C21, etc. are averaged individually to produce a vector:
| c11 |
| c21 |
| c31 |
|  .  |
| cN1 |

That is, for the i-th column:
Cvertstack = [C1i; C2i; C3i; ... ; CNi];
Cvertavg{i}(1) = mean(C1i,'all');
Cvertavg{i}(2) = mean(C2i,'all');
etc.

This is applied for each column to obtain the symmetrized covariance matrix C

This two-step averaging approach provides a better trade-off between
runtime and memory.

In our implementation (2020-12-12), the pairwise distance matrices are
computed, followed by a post-processing step of converting to a covariance
matrix.

%}

%% CODE GRAVEYARD
%{
% npts = size(o,1);
% % qA = normr(o(:,1:4));
% % qB = normr(o(:,5:8));
% nsympairs = size(Spairs,1);
% for i = 1:npts
% %     qAtmp = qA(i,:);
% %     qBtmp = qB(i,:);
%     otmp = o(i,:);
%     oset = osymsets(otmp,pgnum,[],grainexchangeQ,doublecoverQ,epsijk);
% %     oset = osymsets(qAtmp,qBtmp,Spairs);
% end


% %         s1 = (i-1)*npts+1; %start ID1
% %         f1 = i*npts; %finish ID2
%         oset1 = osets(s1:f1,:);
%                 s2 = (j-1)*npts+1; %start ID2
%         f2 = j*npts; %finish ID2
%         oset2 = osets(s2:f2,:);


%% main loop
Cmaster = cell(nsym,nsym);
for i = 1:nsym
    %text waitbar
    if mod(i,nreps2) == 0
        if waitbarQ
            send(D,i);
        end
    end
    oset2 = osets{j};
    for j = i:nsym
        oset1 = osets{i};
        oset2 = osets{j};
        Cmaster{i,j} = pdist2(osetfull,oset2,@get_alen);
    end
end



        %repeat array
%     oset2rep = repmat(oset2,nsym,1);
    %reshape into stacked sub-arrays
%     Csubstack = reshape(Ctmp,npts,nsym);

%     PDvertsplit = shiftdim(PDvertsplit,[2 3 1]); % i.e. 100x576x576 -> 576x576x10


% PDhorz = cat(3,PDvertavg{:});
% pd = mean(PDhorz,3);
Cvertavg{i} = (C1i + C2i + C3i + ... + CNi)/N;

Chorz contains Cvertavg{1}, Cvertavg{2}, Cvertavg{3}, ..., Cvertavg{N}
arranged as an (N * N * nsym) array

Applying the horizontal average:
C = (Cvertavg{1} + Cvertavg{2} + Cvertavg{3} + ... + Cvertavg{N})/N


% osets = get_octpairs(osets,[],epsijk);

% nsympts = size(osets,1); %replaced by nptstot

%}