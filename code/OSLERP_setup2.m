%OSLERP_SETUP2  another script to call OSLERP (deprecated)

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

seed = 10;
rng(seed)

% fList = {'A','B','C','D','E','O'}; %feature list
% name1 = 'interior';
% name2List = {'A'};
% name2List = {'A','B','C','D','E','O'};
name1 = 'O';
name2List = {'A','C'};

for i = 1:length(name2List)
	name2 = name2List{i};

	% t = tiledlayout(5,4);
	qA = qlist.(name1);
	qB = qlist.(name2);
	% qA = normr(qlist.(name1)+0.1*rand(1,4));
	% qB = normr(qlist.(name2)+0.1*rand(1,4));
	
	five(1).q = qA;
	five(2).q = qB;
	five(1).d = dlist.(name1);
	% five(2).d = dlist.(name2);
	five(2).d = q2rod(qB);
	
% 	[~,RA] = symaxis(qA,name1);
% 	nA = (RA*[0 0 1].').';
% 	[~,RB] = symaxis(qB,name2);
% 	nB = (RB*[0 0 1].').';
	
	% n0 = normr([-1 1 1]);
	% nA = normr([n0+0.3*rand(1,3)]);
	% nB = normr([n0+0.3*rand(1,3)]);
	nA = [0 0 1];
	nB = [0 0 1];
	
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
% 	[omega,oct_sym,zeta] = GBdist([o1 o2],32,false);

	
	get_omega = @(o1,o2) 2*acos(abs(sum(o1(1:4).*o2(1:4))-sum(o1(5:8).*o2(5:8)))/2);
	get_norm = @(o1,o2) norm(o1-o2);
	
	o1_sym = oct_sym(1:8);
	o2_sym = oct_sym(9:16);
	% disp([o1_sym; o2_sym])
	% disp(omega)
	get_omega(o1_sym,o2_sym);
	disp(omega)
	disp(get_norm(o1_sym,o2_sym))
	
	o1_IDs = find(ismembertol(octonion_pair_sym_list(:,1:8),o1,1e-3,'ByRows',true));
	o2_IDs = find(ismembertol(octonion_pair_sym_list(:,9:16),o2,1e-3,'ByRows',true));

	o1s = octonion_pair_sym_list(o1_IDs,1:8);
	o2s = octonion_pair_sym_list(o1_IDs,9:16);
	
	minIDs = find(ismembertol(wveclist,omega,1e-6));
	octpairsymlist = octonion_pair_sym_list(minIDs,:);
	eucldist{i} = vecnorm(octpairsymlist(:,1:8)-octpairsymlist(:,9:16),2,2);
	
	minsymlist{i} = uniquetol(round(octpairsymlist,12),1e-6,'ByRows',true);
	
	disp('unique o1s')
	disp(minsymlist{i}(:,1:8))
	disp('corresponding o2s')
	disp(minsymlist{i}(:,9:16))
end

tmp = vertcat(minsymlist{:});
disp(tmp(:,9:16))