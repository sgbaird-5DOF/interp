function parityplot(yactual,ypred,plottype,NV)
arguments
    yactual double
    ypred double
    plottype char {mustBeMember(plottype,{'hex','scatter'})} = 'hex'
    NV.xunits char = 'J/m^2'
    NV.yunits char = 'J/m^2'
    NV.xname char = 'actual GBE'
    NV.yname char = 'predicted GBE'
    NV.title char = ''
    NV.charlbl char = ''
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
    NV.cblbl char = 'counts'
    NV.cbnds(1,2) double = [1 500]
    NV.cscale char {mustBeMember(NV.cscale,{'log','linear'})} = 'log'
    NV.res(1,1) double = 50
    NV.drawEdges(1,1) logical = 0
    NV.showZeros(1,1) logical = 0
    NV.xlim = [min([yactual(:);ypred(:)]) max([yactual(:);ypred(:)])]
    NV.ylim = [min([yactual(:);ypred(:)]) max([yactual(:);ypred(:)])]
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

switch plottype
    case 'hex'
        hexscatter(yactual,ypred,NV.xlim,NV.ylim,'cscale',NV.cscale,'cbnds',[1 500]);
    case 'scatter'
        %% scatter
        if NV.fillQ
            %     ax1 = plot(yactual,ypred,NV.mkr);
            ax1 = scatter(yactual,ypred,NV.sz,NV.c,NV.mkr,'filled');
        else
            %     ax1 = plot(yactual,ypred,NV.mkr);
            ax1 = scatter(yactual,ypred,NV.sz,NV.c,NV.mkr);
        end
        
        scatterNames = fields(NV.scatterOpts);
        for i = 1:length(scatterNames)
            scatterName = scatterNames{i};
            ax1.(scatterName) = NV.scatterOpts.(scatterName);
        end
end

%% refline
hold on
ax2 = refline(1,0);

reflineNames = fields(NV.reflineOpts);
for i = 1:length(reflineNames)
    reflineName = reflineNames{i};
    ax2.(reflineName) = NV.reflineOpts.(reflineName);
end

if ~isempty(NV.xlim)
    xlim(NV.xlim)
end
if ~isempty(NV.ylim)
    ylim(NV.ylim)
end

axis square tight

%% axes
if ~isempty(NV.xunits)
    xunits = ['(',NV.xunits,')'];
else
    xunits = '';
end
if ~isempty(NV.yunits)
    yunits = ['(',NV.yunits,')'];
else
    yunits = '';
end
xlbl = strjoin({NV.xname,xunits},' ');
ylbl = strjoin({NV.yname,yunits},' ');
xlabel(xlbl,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)
ylabel(ylbl,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)
title(NV.title,'Interpreter',NV.Interpreter,'FontSize',NV.FontSize)

% legend
lgd = legend(NV.legend{:});
if ~isempty(lgd)
    lgd.Location = NV.Location;
end

% label for figure tiles, e.g. '(a)', '(b)', '(c)', '(d)'
if ~isempty(NV.charlbl)
    text(0.025,0.95,NV.charlbl,'Units','normalized','FontSize',12)
end

hold off
end