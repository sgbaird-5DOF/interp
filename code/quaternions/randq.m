% RANDQ generates uniformly distributed positive quaternions.
%-------------------------------------------------------------------------%
%Filename:  randq.m
%Author:    Oliver Johnson
%Date:      6/21/2011 (validated 7/6/2011)
%
% Inputs:
%   N - A scalar integer indicated the desired number of quaternions.
%   domain - (optional) A string indicating that the quaternions are to be
%            taken either from 'SO(3)' or from 'S3'.
%
% Outputs:
%   q - An N-by-4 array of random quaternions.
%
% [1] Weisstein, Eric W. "Hypersphere Point Picking." From MathWorld--A
%     Wolfram Web Resource.
%     http://mathworld.wolfram.com/HyperspherePointPicking.html
%-------------------------------------------------------------------------%

function q = randq(N,domain,cs)

if nargin > 0
    assert(isscalar(N) && mod(N,1) == 0,'N must be a scalar integer.')
else
    N = 1;
end
if nargin < 2
    domain = 'SO(3)';
end
if nargin > 2
    assert(any(strcmpi(cs,{'triclinic','c1','1','-1','ci','monoclinic','m','cs','c2h','c2','2/m','2','orthorhombic','mmm','mm2','d2h','d2','c2v','222','trigonal','rhombohedral','d3d','d3','c3v','c3i','c3','3m','32','3','-3m','-3',' -3','tetragonal','s4','d4h','d4','d2d','c4v','c4h','c4','4mm','422','4/mmm','4/m','4','-42m','-4','hexagonal','d6h','d6','d3h','c6v','c6h','c6','c3h','6mm','622','6/mmm','6/m','6','-6m2','-6','cubic','m-3m','m-3','th','td','t','oh','o','432','23','-43m'})), 'Unrecognized crystal system.');
    nsymm = nsymmdis(cs);
end

q = zeros(0,4); %initialize
n = 0;
while n < N
        
    %---limit to domain---%
    switch domain
        case 'SO(3)'
            
            %---sample points---%
            qtemp = randn(N,4);
            
            %---normalize---%
            r = qnorm(qtemp);
            qtemp = qtemp./(r(:,ones(4,1)));
            qtemp(qtemp(:,1) < 0,1) = -qtemp(qtemp(:,1) < 0,1);
            
        case 'S3'
            
            %---sample points---%
            qtemp = randn(N,4);
            
            %---normalize---%
            r = qnorm(qtemp);
            qtemp = qtemp./(r(:,ones(4,1)));
            
        case 'FZ'
            error('Not yet implemented.');
        case 'disorientation'

            %---sample points---%
            if nsymm*(N-n) < 1e7
                qtemp = randn(nsymm*(N-n),4);
            else
                qtemp = randn(1e7,4);
            end
            
            %---normalize---%
            r = qnorm(qtemp);
            qtemp = qtemp./(r(:,ones(4,1)));
            qtemp(qtemp(:,1) < 0,1) = -qtemp(qtemp(:,1) < 0,1); %all disorientations require a >= 0 so save yourself some time
            idx = isdisorientation(qtemp,cs);
            qtemp = qtemp(idx,:);
    end
    
    %---store---%
    q = [q; qtemp]; %#ok<AGROW>
    
    %---check size---%
    n = size(q,1);
end
q = q(1:N,:);

end

% IDENTIFY QUATERNIONS IN THE FUNDAMENTAL ZONE FOR A GIVEN SYMMETRY
function tf = isdisorientation(q,cs)

a = q(:,1); % may be faster if you don't separate, since 2d arrays are stored contiguously
b = q(:,2);
c = q(:,3);
d = q(:,4);

switch lower(cs)
    case {'triclinic','c1','1','-1','ci'}
        tf = (a >= 0) & (d >= 0);
    case {'monoclinic','m','cs','c2h','c2','2/m','2'}
        tf = (a >= b) & (b >= 0) & (d >= 0);
    case {'orthorhombic','mmm','mm2','d2h','d2','c2v','222'}
        tf = (a >= b) & (b >=0) & (a >= c) & (c >=0) & (a >= d) & (d >=0);
    case {'trigonal','rhombohedral','d3d','d3','c3v','c3i','c3','3m','32','3','-3m','-3',' -3'}
        tf = (a >= b) & (b >= sqrt(3)*abs(c)) & (a >= sqrt(3)*d) & (d >= 0);
    case {'tetragonal','s4','d4h','d4','d2d','c4v','c4h','c4','4mm','422','4/mmm','4/m','4','-42m','-4'}
        tf = (a >= b) & (b >= c) & (c >= 0) & (a >= (b+c)/sqrt(2)) & (a >= (sqrt(2)+1)*d) & (d >=0);
    case {'hexagonal','d6h','d6','d3h','c6v','c6h','c6','c3h','6mm','622','6/mmm','6/m','6','-6m2','-6'}
        tf = (a >= b) & (b >= sqrt(3)*c) & (c >= 0) & (a >= (sqrt(3)*b+c)/2) & (a >= (2+sqrt(3))*d) & (d >= 0);
    case {'cubic','m-3m','m-3','th','td','t','oh','o','432','23','-43m'}
        tf = (a >= (sqrt(2)+1)*b) & (a >= b+c+d) & (b >= c) & (c >= d) & (d >= 0);
end

end

function nsymm = nsymmdis(cs)

switch lower(cs)
    case {'triclinic','c1','1','-1','ci'}
        neq = 1; %number of equations used for computing symmetryically equivalent rotations
        nperms = 1;
        nsigns = 4;
    case {'monoclinic','m','cs','c2h','c2','2/m','2'}
        neq = 1;
        nperms = 2;
        nsigns = 8;
    case {'orthorhombic','mmm','mm2','d2h','d2','c2v','222'}
        neq = 1;
        nperms = 4;
        nsigns = 16;
    case {'trigonal','rhombohedral','d3d','d3','c3v','c3i','c3','3m','32','3','-3m','-3',' -3'}
        neq = 9;
        nperms = 2;
        nsigns = 8;
    case {'tetragonal','s4','d4h','d4','d2d','c4v','c4h','c4','4mm','422','4/mmm','4/m','4','-42m','-4'}
        neq = 2;
        nperms = 8;
        nsigns = 16;
    case {'hexagonal','d6h','d6','d3h','c6v','c6h','c6','c3h','6mm','622','6/mmm','6/m','6','-6m2','-6'}
        neq = 9;
        nperms = 4;
        nsigns = 16;
    case {'cubic','m-3m','m-3','th','td','t','oh','o','432','23','-43m'}
        neq = 6;
        nperms = 24;
        nsigns = 16;
end

nsymm = neq*nperms*nsigns;

end