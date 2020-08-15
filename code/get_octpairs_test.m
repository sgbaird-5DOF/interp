%get_octpairs test
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
% Notes:
%
%--------------------------------------------------------------------------
%% load files
clear; close all;

addpathdir({'q2rod.m','GBfive2oct.m','ax2qu.m'})

% fname = '5DOF_vtx_deleteOz_octsubdiv1.mat';
fname = '5DOF_vtx_octsubdiv1.mat';
load(fname,'pts')
savename = 'octvtx_pairmin.mat';
NV = {'o2addQ',false,'plotQ',true};
[octvtx,usv,five] = get_octpairs(pts,savename,NV{:});
% [octvtx,usv,five,omega23_pair,omega23_GBdist] = get_octpairs2(pts,five,savename,o2addQ,plotQ);

[vtxnames,ia] = sort({five.geometry});
octvtx = octvtx(ia,:);

%% pairwise distance matrices
npts2 = size(octvtx,1);
% f = waitbar(0);
pd1 = zeros(npts2);
pd2 = pd1;

% tmp = num2cell(octvtx,2);
% octvtx_pairs = allcomb(tmp,tmp);
% octvtx1 = vertcat(octvtx_pairs{:,1});
% octvtx2 = vertcat(octvtx_pairs{:,2}); %if I want to use the vectorized versions, I'll need to deal with preserving PD-formating

parfor i = 1:npts2
% 	waitbar(i/npts2,f);
	for j = 1:npts2
		pd1(i,j) = get_omega(octvtx(i,:),octvtx(j,:)); %pairwise distance matrix of omega23_pair
		pd2(i,j) = GBdist([octvtx(i,:) octvtx(j,:)],32,false); % 32 == Oh point group
	end
end
% close(f)

%get rid of imaginary part, if applicable
pd1 = rad2deg(abs(pd1));
pd2 = rad2deg(abs(pd2));

%% plotting
fig = figure;
fig.Position = [249.0000  134.5000  849.5000  668.0000];
% nexttile
% plot(omega23_pair)
% hold on
% plot(omega23_GBdist)
% % axis tight
% ax = gca;
% ax.YLim(1) = 0;
% legend('oct pair','GBdist','Location','southeast')
% xlabel('vertex #')
% ylabel('\omega (rad) rel. to point \it{O}\rm (+z)')

nexttile
histogram(reshape(pd1,1,[]))
hold on
histogram(reshape(pd2,1,[]))
legend('oct pair','GBdist','Location','best')
xlabel('\omega (deg)')
ylabel('counts')

nexttile
histogram(reshape(pd1-pd2,1,[]))
legend('excess','Location','best')
xlabel('\omega (deg)')
ylabel('counts')

nexttile
pcolor(pd1)
cb = colorbar;
cb.Label.String = 'pair-referenced \omega (deg)';
% clims = cb.Limits;
clims = [0 60];
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')
% caxis([0 pi/2])

nexttile
pcolor(pd1-pd2)
cb = colorbar;
cb.Label.String = 'excess \omega (deg)';
% caxis([0 pi/2])
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')


%%
%---------------------------CODE GRAVEYARD---------------------------------
%{
% ic2 = ic12(ia23(minID23,1));
% ic3 = ic13(ia23(minID23,2));

% ic2 = ia23(minID23,1);
% ic3 = ia23(minID23,2);

name2list = {'A','B','C','D','E','O'};

name3list = {'A','B','C','D','E','O'};

namepairlist = nchoosek({'A','B','C','D','E','O'},2);

for i = 1:length(namepairlist)
	
	namepair = namepairlist{i,:};
	name2 = namepair{1};
	name3 = namepair{2};

end

load('misFZfeatures.mat','qlist','dlist')

S = load('vtxmesh.mat','meshpts');
pts = S.meshpts;


	[~,RC] = symaxis(qC,name3);
	nC = (RC*[0 0 1].').';

	% n0 = normr([1 1 1]);
	% nA = normr([nA+0.1*rand(1,3)]);
	% nB = normr([n0+0.3*rand(1,3)]);
	% nC = normr([n0+0.3*rand(1,3)]);
	% nA = [0 0 1];
	% nB = [0 0 1];
	% nC = [0 0 1];

	% %load octonions
	% S = load('olist.mat','olist');
	% olist = S.olist;
	%
	% o1 = olist(1,:);
	% o2 = olist(2,:);
	% o3 = olist(3,:);

	% o1 = [-1 0 0 0 1 0 0 0];
% 	o2 = GBfive2oct(qB,nB);
% 	o3 = GBfive2oct(qC,nC);
	
	% get_omega(o1,o2)
	% get_omega(o1,o3)
	% get_omega(o2,o3)


	qB = qlist.(name2);
	qC = qlist.(name3);
	
	seed = 10;
	rng(seed);




% pairlist = nchoosek(1:npts,2); %all unique pairs

get_omega = @(o1,o2) 2*acos(abs( dot(o1(1:4),o2(1:4)) + dot(o1(5:8),o2(5:8)) )/2);
get_omega2 = @(o1,o2) 2*acos(abs(sum(o1(1:4).*o2(1:4))-sum(o1(5:8).*o2(5:8)))/2);

	

%remove duplicate rows (low tol OK b.c. matching against 16 numbers)
[min13,ia13,ic13] = uniquetol(round(octpairsymlist13,15),1e-3,'ByRows',true);
[min12,ia12,ic12] = uniquetol(round(octpairsymlist12,15),1e-3,'ByRows',true);


clabel('\omega (rad)')


%		GBdistdis.m

nexttile
histogram(omega23_pair)
hold on
histogram(omega23_GBdist)
legend('oct pair','GBdist','Location','best')
xlabel('\omega (rad)')
ylabel('counts')
%}
