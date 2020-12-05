function [qd,negQ,invQ,disSpair] = qdis(q,pgnum,epsijk)
arguments
    q(:,4)
    pgnum(1,1) double = 32 %only set up for cubic 2020-12-01
    epsijk(1,1) double = 1
end

%% setup
npts = size(q,1);
qd = zeros(npts,4);

%% map to disorientation
for i = 1:npts
    % get symmetrically equivalent quaternions
    Spairs = get_sympairs(pgnum);
    
    % separate symmetry pairs
    S1 = Spairs(:,1:4);
    S2 = Spairs(:,5:8);
    
    % apply symmetry operators
    nsympairs = size(S1,1);
    qrep = repmat(q(i,:),nsympairs,1);
    qset1 = qmult(S1,qmult(qrep,S2,epsijk),epsijk); %passive convention?
    qset2 = qmult(S1,qmult(qinv(qrep),S2,epsijk),epsijk);
    % add inverse as well (should check with Grimmer H. Acta crystallographica section a 1980;36:382-389)
%     qset = [qset; qinv(qset)]; %#ok<AGROW>
    qset = [qset1; qset2];
    
    %rodrigues vectors
    %     r = qu2ro(qset,epsijk);
    r = q2rod(qset);
    
    % disorientation IDs
    isdis = inmisFZ(r,1e-6);
    
    % take first non-zero entry
    disID = find(isdis);
    disID = disID(1);
    
    % assign disorientation
    qd(i,:) = qset(disID,:);
    if qd(i,1) < 0
        %convention that q0 is positive
        qd(i,:) = -qd(i,:);
        negQ = true;
    else
        negQ = false;
    end
    
    if disID > nsympairs
        invQ = true;
        disID2 = disID - nsympairs;
    else
        invQ = false;
        disID2 = disID;
    end
    disSpair = Spairs(disID2,:);
    
end

%% CODE GRAVEYARD
%{
%% load symmetry operators
%load operators
% symops = load('PGsymops.mat');
% for point group symmetry names, see PGsymnames.mat
%unpack point group
qpt = symops.Q{pgnum};

%     S1 = [S1; -S1];
%     S2 = [S2; -S2];
%}