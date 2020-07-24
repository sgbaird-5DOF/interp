%function nA_out = toBPFZ(qlist,nAlist)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-15
%
% Description: Rotate an arbitrary boundary plane normal for a given
% quaternion into a standard boundary plane fundamental zone, similar to
% what disorientation.m does for misorientations.
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

nint = 1;
ctrcuspQ = 1;

pgnum = 32; %cubic Oh
symnames = load('PGnames.mat'); %need to add crystal_symmetry_ops to path in order for this to work
symops = load('PGsymops.mat');
qpt = symops.Q{pgnum}; %choose point group symmetry

t = num2cell(qpt,2);

%convert quaternions to rotation matrices
Rmats = cellfun(@(q) qu2om(q),t,'UniformOutput',false);

npts = size(qlist,1);
for i = 1:npts
	%unpack
	nA = nAlist(i,:)
	Rmat = Rmats{i};
	
	%get mesh points
	[pts,A,R,TRI,af] = meshBP(q,nint,ctrcuspQ);
	t = num2cell(pts);
	[az1,el1,~] = cart2sph(t{:});
	
	%define bounds
	azmin = min(az1);
	elmin = min(az1);
	
	azmax = max(el1);
	elmax = max(az1);
	
	%loop through rotation matrices
	%initialize
	nR = length(Rmats);
	az = zeros(nR,1);
	el = az;
	nArot = zeros(nR,3);
	for j = 1:nR
		nArot(j,:) = (Rmat*nA.').';
	end
	
	%convert to spherical
	t = num2cell(nArot,1);
	[az(j),el(j),~] = cart2sph(vertcat(t{:}));
	
	%test bounds
	(azmin <= az) && (az <= azmax)
	nA_out(i,:) = nArot(:,i,:);
	
end

t=num2cell(nArot,1);
plot3(t{:},'*')


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
%}

