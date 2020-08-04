%optimize_zeta test

%% setup
clear; close all

%add to path
addpathdir({'misFZfeatures.mat','PGnames.mat','qu2ax.m'})

%set rng
seed = 10;
rng(seed);

%random octonions
npts = 30;
o = get_ocubo([],'uniform',4,10);

method = 2;
disp(['method: ' int2str(method)])
switch method
	case 1
		%symmetrized octonions
		o1 = o(1,:);
		o1rep = repmat(o1,npts-1,1);
		[w,o2syms] = GBdist4(o1rep,o(2:end,:),32,'norm');
		
		% unpack symmetrized octonions
		o2out = zeros(npts-1,8);
		for i = 1:npts-1
			o2tmp = o2syms{i};
			if size(o2tmp,1) > 1
				nm = vecnorm(o2tmp-repmat(o2syms{end}(1,:),2,1),2,2);
				[~,id] = min(nm);
				o2out(i,:) = o2tmp(id,:);
				% 		o2out(i,:) = o2tmp(o2tmp(:,1) > o2tmp(:,5),:);
				% 		warning('more than one octonion output')
			else
				o2out(i,:) = o2tmp;
			end
		end
		
		% repackage
		onew = [o(1,:); o2out];
	case 2
		onew = get_octpairs(o);
end

onew = uniquetol(onew,'ByRows',true);

onewtmp1 = proj_down(onew,1e-5);
[~,~,onewtmp2] = hsphext_subdiv(onewtmp1,1,true);
ids = ismembertol(onewtmp1,onewtmp2,'ByRows',true);
onew = onew(ids,:);

npts = size(onew,1);


%% optimization
dtype = 'norm';
[z,errmin,exitflag,output] = optimize_zeta(onew,dtype,[])

% initial data
[err0,~,pdzero,pdtrue] = pd_sse(onew,zeros(npts,1),dtype);

% final error
[err,of,pd] = pd_sse(onew,z,dtype,pdtrue);

% convert to pairwise distance matrix
IDs = allcomb(1:npts,1:npts);
pd0 = zeros(npts,npts);
pd1 = zeros(npts,npts);
pd2 = zeros(npts,npts);
idx = sub2ind([npts npts],IDs(:,1),IDs(:,2));
for i = 1:npts^2
	pd0(i) = pdzero(idx(i));
	pd1(i) = pd(idx(i));
	pd2(i) = pdtrue(idx(i));
end
% 
% pd0 = flipud(abs(pd0));
% pd1 = flipud(abs(pd1));
% pd2 = flipud(abs(pd2));
if strcmp(dtype,'omega')
	pd0 = rad2deg(pd0);
	pd1 = rad2deg(pd1);
	pd2 = rad2deg(pd2);
end
	
%% plotting
% setup
fig=figure;
fig.Position=[269.0000  183.5000  694.5000  552.0000];
tiledlayout(2,2)
cmax = 40;
switch dtype
	case 'omega'
		str = '\omega';
		unit = ' (deg)';
		clim = [0 cmax];
	case 'norm'
		str = 'norm';
		unit = '';
		clim = [0 deg2rad(cmax)];
end

% initial distances
nexttile
imagesc(pd0)
ax=gca;
ax.XAxisLocation = 'top';
cb = colorbar;
cb.Label.String = [str '_0' unit];
clims = clim;
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')
axis equal tight

% initial excess
nexttile
imagesc(pd0-pd2)
ax=gca;
ax.XAxisLocation = 'top';
cb = colorbar;
cb.Label.String = ['excess ' str '_0' unit];
clims = clim;
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')
axis equal tight

% final distances
nexttile
imagesc(pd1)
ax=gca;
ax.XAxisLocation = 'top';
cb = colorbar;
cb.Label.String = [str '_f' unit];
clims = clim;
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')
axis equal tight

% final excess
nexttile
imagesc(pd1-pd2)
ax=gca;
ax.XAxisLocation = 'top';
cb = colorbar;
cb.Label.String = ['excess ' str '_0' unit];
clims = clim;
caxis(clims);
xlabel('vertex #')
ylabel('vertex #')
axis equal tight




%----------------------------CODE GRAVEYARD--------------------------------
%{
% "true" pairwise distances
IDs = allcomb(1:npts,1:npts);
o1 = o(IDs(:,1),:);
o2 = o(IDs(:,2),:);
pdtrue = GBdist4(o1,o2,32,dtype);
%}
