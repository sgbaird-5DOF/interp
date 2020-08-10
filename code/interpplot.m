function interpplot(fname)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-28
%
% Description: Create parity plot with 5DOF plots on side
% 
% Inputs: fname - filename that contains relevant variables for plotting
%
% Outputs: none
%
% Usage: interpplot(fname)
%
% Dependencies: plot5DOF.m
%
% Notes:
%
%--------------------------------------------------------------------------

%load info
addpathdir({fname})
vars = {'mesh','data','datainterp','nnID','ilist','nndistList','nonintDists',...
	'interpSE','nnSE','interpRMSE','nnRMSE'};
load(fname,vars{:})

%figure setup
fig = figure;
fig.Position = [125,125,1000,600];

t = tiledlayout(2,3);
nexttile(1)

xmin = min([mesh.props(nnID);data.props],[],'all');
xmax = max([mesh.props(nnID);data.props],[],'all');

xlims = [xmin,xmax];

alphaval = 0.5;
%parity plot
scatter(data.props,datainterp,2,'k','filled','markerfacealpha',alphaval);
hold on
% axis tight
scatter(data.props(ilist),mesh.props(nnID),2,'r','filled','markerfacealpha',alphaval);

xlabel('BRK Energy')
ylabel('interpolated BRK Energy')
title(['RMSE (J/m^2): interp == ' num2str(interpRMSE,'%3.4f') ', NN == ' num2str(nnRMSE,'%3.4f')])

plot(xlims,xlims,'c')
axis tight

legend({'interp data','NN data'},'Location','northoutside')


% distance histogram
nexttile(4)

% specify number of bins and edges of those bins; this example evenly spaces bins
nbins = 20;
edges = linspace(0,max([interpSE;nnSE]),nbins);

histogram(interpSE,edges);
hold on
histogram(nnSE,edges);

title(['# non-intersections: ' int2str(length(nnSE)) ...
	'/' int2str(length(nnSE)+length(interpSE))]);
legend('intersecting','non-intersecting','Location','northeast')
xlabel('multi-GB symmetrized SE (J/m^2)')
ylabel('counts')
% 5DOF plots

%convert to disorientation
meshtempq = disorientation(vertcat(mesh.five.q),'cubic');
datatempq = disorientation(vertcat(data.five.q),'cubic');

meshtempd = q2rod(meshtempq);
datatempd = q2rod(datatempq);

t = num2cell(meshtempq,2);
[mesh.five.q] = t{:};
t = num2cell(datatempq,2);
[data.five.q] = t{:};

t = num2cell(meshtempd,2);
[mesh.five.d] = t{:};
t = num2cell(datatempd,2);
[data.five.d] = t{:};

if contains(fname,'ocubo')
	plot5DOF(mesh.five,'mesh')
	plot5DOF(data.five,'data',[],ilist)
else
	plot5DOF(mesh.five,'mesh',mesh.opts)
	plot5DOF(data.five,'data',data.opts,ilist)
end

nptsmesh = size(mesh.pts,1);
nptsdata = size(data.pts,1);

% nexttile(5)
% legend('intersecting','non-intersecting','Location','northoutside')
% nexttile(6)
% legend('intersecting','non-intersecting','Location','northoutside')

sgtitle({...
	['mesh == ' mesh.sampleMethod '_octsubdiv' int2str(mesh.opts.octsubdiv) ', npts == ' int2str(nptsmesh)],...
	['data == ' data.sampleMethod '_octsubdiv' int2str(data.opts.octsubdiv) ', npts == ' int2str(nptsdata)]}, ...
	'Interpreter','none','FontSize',10);

fpath = fullfile('figures',fname(1:end-4));
try
	print(fpath,'-dpng')
catch
	disp('attempting to print figure one more time')
	print(fpath,'-dpng')
end
savefig(fig,[fpath '.fig'],'compact')


%-------------------------CODE GRAVEYARD-----------------------------------
%{
% distance histogram
nexttile(4)
ax1 = histogram(nndistList);
% ax1.BinLimits(1) = 0;
hold on
ax2 = histogram(nonintDists);

legend('intersecting','non-intersecting','Location','northoutside')
xlabel('multi-GB symmetrized \omega (rad)')
ylabel('counts')

%calculate SE and RMSE values
% ids = ~isnan(datainterp);
% interpSE = (data.props(ids)-datainterp(ids)).^2;
% nnSE = (data.props(ilist)-mesh.props(nnID)).^2;
% interpRMSE = sqrt(mean(interpSE));
% nnRMSE = sqrt(mean(nnSE));

% interpSE = [interpSE; NaN(length(nnSE),1)];
% histogram(interpSE,'Normalization','probability');
% hold on
% nnSE = [nnSE; NaN(length(interpSE),1)];
% histogram(nnSE,'Normalization','probability');

% use histcounts and specify your bins
% cntA = histcounts(interpSE,'BinEdges',edges);
% cntB = histcounts(nnSE,'BinEdges',edges);

% % convert bin edges into bin centers
% b = edges(1:end-1)+diff(edges)/2;
% % use bar
% bar(b,[cntA',cntB'],'stacked')

%}
