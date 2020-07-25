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
	'meshopts','dataopts','nndistList','nonintDists','meshMethod',...
    'dataMethod'};
load(fname,vars{:})

%figure setup
fig = figure;
fig.Position = [125,125,1000,600];

t = tiledlayout(2,3);

%
nexttile(1)

xmin = min([mesh.props(nnID);data.props],[],'all');
xmax = max([mesh.props(nnID);data.props],[],'all');

xlims = [xmin,xmax];

%parity plot
scatter(data.props,datainterp,2,'k','filled','markerfacealpha',0.1);
hold on
% axis tight
scatter(data.props(ilist),mesh.props(nnID),2,'r','filled','markerfacealpha',0.1);

%calculate SE and RMSE values
ids = ~isnan(datainterp);
interpSE = (data.props(ids)-datainterp(ids)).^2;
nnSE = (data.props(ilist)-mesh.props(nnID)).^2;
interpRMSE = sqrt(mean(interpSE));
nnRMSE = sqrt(mean(nnSE));

xlabel('BRK Energy')
ylabel('interpolated BRK Energy')
title(['RMSE (rad): interp == ' num2str(interpRMSE) ', NN == ' num2str(nnRMSE)])

plot(xlims,xlims,'c')
axis tight

legend({'interp data','NN data','parity'},'Location','northoutside')


% distance histogram
nexttile(4)
histogram(interpSE);
hold on
histogram(nnSE);

legend('intersecting','non-intersecting','Location','northoutside')
xlabel('multi-GB symmetrized \omega (rad)')
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

plot5DOF(mesh.five,'mesh',meshopts)
plot5DOF(data.five,'data',dataopts,ilist)

nptsmesh = size(mesh.pts,1);
nptsdata = size(mesh.data,1);

nexttile(5)
legend('intersecting','non-intersecting','Location','northoutside')
nexttile(6)
legend('intersecting','non-intersecting','Location','northoutside')

sgtitle({...
	['mesh == ' meshMethod '_octsubdiv' int2str(meshopts.octsubdiv) ', npts == ' int2str(nptsmesh)],...
	['data == ' dataMethod '_octsubdiv' int2str(dataopts.octsubdiv) ', npts == ' int2str(nptsdata)]}, ...
	'Interpreter','none','FontSize',10);

fpath = fullfile('figures',fname(1:end-4));
print(fpath,'-dpng')
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
%}
