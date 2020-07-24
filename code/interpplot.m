function interpplot(fname)
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
% Usage:
%
% Dependencies:
%
% Notes:
%
%--------------------------------------------------------------------------

%load info
addpathdir({fname})
vars = {'mesh','data','datainterp','nnID','psdata','ilist',...
	'meshopts','dataopts','nndistList','nonintDists'};
load(fname,vars)

%figure setup
fig = figure;
fig.Position = [125,125,1000,600];

t = tiledlayout(2,3);

%
nexttile(1)

xmin = min([mesh.props(nnID);data.props],[],'all');
xmax = max([mesh.props(nnID);data.props],[],'all');

xlims = [xmin,xmax];

scatter(data.props,datainterp,5,'k','markerfacealpha',1);
hold on
% axis tight
scatter(data.props(ilist),mesh.props(nnID),5,'r','filled','markerfacealpha',1);

xlabel('BRK Energy @ datapoint')
ylabel('interpolated BRK Energy')

plot(xlims,xlims,'c')

legend({'interp data','NN data'},'Location','northoutside')


% distance histogram
nexttile(4)
ax1 = histogram(nndistList);
% ax1.BinLimits(1) = 0;
hold on
ax2 = histogram(nonintDists);

legend('intersecting','non-intersecting','Location','northoutside')
xlabel('7D Euclidean NN distance')
ylabel('counts')

% 5DOF plots
plot5DOF(mesh.five,'mesh',meshopts)
plot5DOF(data.five,'data',dataopts,ilist)

sgtitle({...
	['mesh == ' meshMethod '_octsubdiv' int2str(meshopts.octsubdiv)],...
	['data == ' dataMethod '_octsubdiv' int2str(dataopts.octsubdiv)]}, ...
	'Interpreter','none','FontSize',10);

print(meshdata.fname(1:end-4),'-dpng')
savefig(meshdata.fname(1:end-4))