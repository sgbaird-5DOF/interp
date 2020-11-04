% OSLERP_SETUP  script to call OSLERP (deprecated)

%% setup
clear; close all;

fnames = {'PGnames.mat','misFZfeatures.mat','q2rod.m','ax2qu.m','GBfive2oct.m'};
for i = 1:length(fnames)
	fname = fnames{i};
	fpath = fullfile('**',fname);
	file = dir(fpath);
	if ~isempty(file)
		addpath(file.folder);
	end
end

load('misFZfeatures.mat','qlist','dlist')

% fList = {'A','B','C','D','E','O'}; %feature list
fList = {'A','B','C'};

%generate 5DOF mesh

% five = struct;
% five(1).q = {};
% five(1).nA = {};
% five(1).d = {};
% for i = 1:length(fList)
% 	feature = fList{i};
% 	q = qlist.(feature);
% 	[A,R] = symaxis(q,feature);
% 	[nlist,nK{i}] = orthoplex(3); %normals & triangulation
% 	for j = 1:size(nlist,1)
% 		five.q{i,j} = q;
% 		nA = nlist(j,:);
% 		five.nA{i,j} = (A*(nA.')).'; %see meshBP.m
% 		five.d{i,j} = q2rod(q);
% 		
% 		o{i,j} = GBfive2oct(q,nA);
% 	end
% end
fig = figure;
fig.Position = [100  50 1000 750];
% t = tiledlayout(5,4);

seed = 10;
rng(seed)

name1 = 'interior';
name2 = 'D';
qA = qlist.(name1);
qB = qlist.(name2);
% qA = normr(qlist.(name1)+0.1*rand(1,4));
% qB = normr(qlist.(name2)+0.1*rand(1,4));

five(1).q = qA;
five(2).q = qB;
five(1).d = dlist.(name1);
% five(2).d = dlist.(name2);
five(2).d = q2rod(qB);

[~,RA] = symaxis(qA,name1);
nA = (RA*[0 0 1].').';
[~,RB] = symaxis(qB,name2);
nB = (RB*[0 0 1].').';

% n0 = normr([-1 1 1]);
% nA = normr([n0+0.3*rand(1,3)]);
% nB = normr([n0+0.3*rand(1,3)]);
% nA = [0 0 1];
% nB = [0 0 1];

five(1).nA = nA;
five(2).nA = nB;

% plot5DOF(five,'original')

%convert to octonions
o1 = GBfive2oct(qA,nA);

% name1 = 'identity';
% o1 = [1 0 0 0 0 0 0 0];
o2 = GBfive2oct(qB,nB);

%calculate distance
% [omega, oct_sym, zeta] = GBdistdis([o1 o2],32,false);
[omega,oct_sym,zeta,wveclist,octonion_pair_sym_list] = GBdistdis([o1 o2],32,false);

get_omega = @(o1,o2) 2*acos(abs(sum(o1(1:4).*o2(1:4))-sum(o1(5:8).*o2(5:8)))/2);

o1_sym = oct_sym(1:8);
o2_sym = oct_sym(9:16);
% disp([o1_sym; o2_sym])
% disp(omega)
get_omega(oct_sym(1:8),oct_sym(9:16));

o1_IDs = find(ismembertol(octonion_pair_sym_list(:,1:8),o1,1e-12,'ByRows',true));
o2s = octonion_pair_sym_list(o1_IDs,9:16);

minIDs = find(ismembertol(wveclist,omega,1e-6));
octpairsymlist = octonion_pair_sym_list(minIDs,:);
eucldist = vecnorm(octpairsymlist(:,1:8)-octpairsymlist(:,9:16),2,2);

disp('unique o1s')
disp(uniquetol(octpairsymlist(:,1:8),1e-12,'ByRows',true))

nsyms = size(octpairsymlist,1);
for i = 1:nsyms
	oct_sym = octpairsymlist(i,:);
	%interpolate octonions
	o1_sym = oct_sym(1:8);
	o2_sym = oct_sym(9:16);
	oct_interp = OSLERP(o1_sym,o2_sym,omega,30);
	
	%convert back to 5DOF
	nocts = size(oct_interp,1);
	% five(nocts) = struct;
	five(1).q = [];
	five(1).nA = [];
	five(1).d = [];
	for j = 1:nocts
		oct = oct_interp(i,:);
		five(j) = GBoct2five(oct);
	end
	
	%plot
	plot5DOF(five,'OSLERP')
	
	five2(1).q = [];
	five2(1).nA = [];
	five2(1).d = [];
	
	qdis = disorientation(vertcat(five.q),'cubic');
	tmp = num2cell(qdis,2);
	[five2.q] = tmp{:};
	for j = 1:length(five2)
		qm = qmult(five(i).q,qinv_francis(qA));
		Rn = qu2om(qm);
		five2(j).d = q2rod(five(i).q);
		five2(j).nA = (Rn*five(i).nA.').';
	end
% 	plot5DOF(five2,'OSLERP')
	
	sgtitle([name1 ' to ' name2])
	
% 	pause
% 	drawnow
end

tmp = uniquetol(abs(octpairsymlist(:,1:8)),1e-15);
tmp2 = uniquetol(abs(octpairsymlist(:,9:16)),1e-15);
tmp3 = uniquetol(abs(octpairsymlist),1e-15);
disp(size(tmp))
disp(size(tmp2))
disp(size(tmp3))



%------------------------------CODE GRAVEYARD------------------------------
%{
nAlist = normr([1 1 1]);
nB = nA;


pA = qmult(qm,qA);
ax = qu2ax(pA);

axisA = ax(1:3);
phiA = ax(4);

%}
