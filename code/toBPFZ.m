function [nAout,nArot,nR] = toBPFZ(qlist,nAlist,NV)
arguments
    qlist(:,4)
    nAlist(:,3)
    NV.nint(1,1) double = 1
    NV.ctrcuspQ(1,1) logical = true
    NV.pgnum(1,1) double = 32 %cubic Oh
    NV.plotQ(1,1) logical = false
end
% TOBPFZ  Rotate an arbitrary boundary plane normal for quaternion into standard boundary plane fundamental zone
%  similar to what disorientation.m does for misorientations.
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-15
%
% Inputs:
%		qlist === quaternion
%
%		nAlist === GB normal
%
% Outputs:
%
%		nA_out === corrected grain boundary normal
%
% Dependencies:
%		meshBP.m
%			findgeometry.m
%			symaxis.m
%
%		qu2om.m
%
%		PGsymops.mat, PGnames.mat (from Grain Boundary Octonion code)
%
% Notes:
%		Implemented for Oh crystal symmetry only (2020-07-15). Computes the
%		24 rotations of the Oh point group and outputs the one that meets the
%		spherical coordinate inequalities defined by the BP FZ spherical
%		wedge.
%--------------------------------------------------------------------------
symnames = load('PGnames.mat'); %need to add crystal_symmetry_ops to path in order for this to work
symops = load('PGsymops.mat');
qpt = symops.Q{NV.pgnum}; %choose point group symmetry

t = num2cell(qpt,2);

%convert quaternions to rotation matrices
Rmats = cellfun(@(q) qu2om(q),t,'UniformOutput',false);

npts = size(qlist,1);
nAout = zeros(npts,3);

nR = length(Rmats);
nArot = repmat({zeros(2*nR,3)},npts,1);

for i = 1:npts
    %unpack
    q = qlist(i,:);
    nA = nAlist(i,:);
    
    %get mesh points
    [pts,A,R,TRI,af] = meshBP(q,NV.nint,NV.ctrcuspQ);
    t = n2c(vertcat(pts.sub));
    [az1,el1,~] = cart2sph(t{:});
    
    %define bounds
    azmin = min(az1);
    elmin = min(az1);
    
    azmax = max(el1);
    elmax = max(az1);
    
    %loop through rotation matrices
    %initialize
    for j = 1:nR
        nArot{i}(j,:) = (Rmats{j}*nA.').';
        nArot{i}(j+nR,:) = (Rmats{j}*-nA.').';
    end
    nArot{i} = uniquetol(nArot{i},'ByRows',true);
    
    %convert to spherical
    t = n2c(nArot{i});
    [az,el,~] = cart2sph(t{:});
    
    %test bounds
    azcheck = (azmin <= az) & (az <= azmax);
    elcheck = (elmin <= el) & (el <= elmax);
    keepIDs = find(azcheck & elcheck);
    switch size(keepIDs,1)
        case 0
            warning(['no triangle found for iteration ' int2str(i)])
        case 1
            nAout(i,:) = nArot{i}(keepIDs,:);
        case 2
            warning(['multiple triangles found for iteration ' int2str(i)])
            nAout(i,:) = nArot{i}(keepIDs(1),:);
    end
    if NV.plotQ
        figure
        %         sphplot()
        %         hold on
        t=n2c(nArot{i});
        %         plot3(t{:},'*')
        test_voronoisphere(nArot{i}.');
    end
end


%------------------------------CODE GRAVEYARD------------------------------
%{
% pgname = symnames.PG_names{pgnum};
% disp(['loading point group: ' pgname])

% Rcombs = allcomb(Rmats,Rmats);

% for i = 1:length(Rcombs)
% 	nA(i,:) = (Rcombs{i,1}*Rcombs{i,2}*az.').';
% end

z_axis = [0 0 1];

		nArot(:,i,j) = (Rmats{j}*z_axis.').';


% 	[az,el] = deal(zeros(nR,1));

%}

