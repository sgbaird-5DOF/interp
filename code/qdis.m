function qd = qdis(q,pgnum,epsijk)
arguments
    q(:,4)
    pgnum(1,1) double = 32 %only set up for cubic 2020-12-01
    epsijk(1,1) double = 1
end

%% load symmetry operators
%load operators
symops = load('PGsymops.mat');
%unpack point group
qpt = symops.Q{pgnum};
nsym = size(qpt,1);

%% setup
npts = size(q,1);
qd = zeros(npts,4);

%% map to disorientation
for i = 1:npts
    % get symmetrically equivalent quaternions
    Spairs = get_sympairs(pgnum);
    S1 = Spairs(:,1:4);
    S2 = Spairs(:,5:8);
%     S1 = [S1; -S1];
%     S2 = [S2; -S2];
    nsympairs = size(S1,1);
    qrep = repmat(q(i,:),nsympairs,1);
    qset = qmult(S1,qmult(qrep,qinv(S2),epsijk),epsijk); %passive convention?
    qset = [qset; qinv(qset)]; %#ok<AGROW>
    % disorientation IDs
%     r = qu2ro(qset,epsijk);
    r = q2rod(qset);
    isdis = inmisFZ(r,1e-3);
    
    % take first non-zero entry
    disID = find(isdis);
    disID = disID(1);
    
    % assign disorientation
    qd(i,:) = qset(disID,:);
    if qd(i,1) < 0
        qd(i,:) = -qd(i,:);
    end
end