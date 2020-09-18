function parityplot(yactual,ypred,NV)
arguments
    yactual double
    ypred double
    NV.units char = 'J/m^2'
    NV.xname char = 'actual GBE'
    NV.yname char = 'predicted GBE'
    NV.title = '';
    NV.sz(1,1) int32 = 36
    NV.c double = [0 0 1]
    NV.mkr char = 'o'
    NV.fillQ(1,1) logical = false
    NV.scatterOpts(1,1) struct = struct()
    NV.reflineOpts(1,1) struct = struct()
    NV.Interpreter char = 'tex'
    NV.FontSize int32 = 11
    NV.legend = {'off'}
    NV.Location char = 'northwest'
end
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
%
% Date: 2020-09-03
%
% Description: Create a parity plot and pass options to scatter() and
% refline().
%
% Inputs:
%  y1 - x-axis data
%
%  y2 - y-axis data
%
%  NV - name-value pairs
%    scatter-specific
%      'sz', 'c', 'mkr', 'fillQ', 'scatterOpts'
%    refline-specific
%      'reflineOpts'
%    xlabel/ylabel/title/legend -specific
%      'Interpreter', 'FontSize', 'legend', 'Location'
%
% Outputs:
%  ax - handle to axis
%
% Usage:
%  parityplot(y1,y2);
%
% Dependencies:
%  *
%
% Notes:
%  See also: SCATTER, REFLINE
%
%  Look at arguments ... end in function to see classes and defaults for
%  name-value pairs
%--------------------------------------------------------------------------
assert(all(size(ypred)==size(yactual)),['y1 [' num2str(size(ypred)) '] and y2 [' num2str(size(yactual)) '] should be same size'])

%% scatter
if NV.fillQ
    ax1 = scatter(yactual,ypred,NV.sz,NV.c,NV.mkr,'filled');
else
    ax1 = scatter(yactual,ypred,NV.sz,NV.c,NV.mkr);
end

scatterNames = fields(NV.scatterOpts);
for i = 1:length(scatterNames)
    scatterName = scatterNames{i};
    ax1.(scatterName) = NV.scatterOpts.(scatterName);
end

%% refline
hold on
ax2 = refline(1,0);

reflineNames = fields(NV.scatterOpts);
for i = 1:length(reflineNames)
    reflineName = reflineNames{i};
    ax2.(reflineName) = NV.scatterOpts.(reflineName);
end

axis square

%% axes
xlbl = strjoin({NV.xname,NV.units},' ');
ylbl = strjoin({NV.yname,NV.units},' ');
xlabel(xlbl,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)
ylabel(ylbl,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)
title(NV.title,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)

lgd = legend(NV.legend{:});
if ~isempty(lgd)
    lgd.Location = NV.Location;
end
hold off
end