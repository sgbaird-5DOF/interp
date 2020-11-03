function multiparity(parity,ID,plottype,NV)
arguments
    parity(:,1) cell
    ID
    plottype char {mustBeMember(plottype,{'hex','scatter'})} = 'hex'
    NV.charlblQ(1,1) logical = true
    NV.titleQ(1,1) logical = true
    NV.xlim = []
    NV.ylim = []
    NV.Interpreter char {mustBeMember(NV.Interpreter,{'tex','latex','none'})} = 'latex'
end

% fig=figure;
% fig.Position = [418.6 194.6 644.8 555.2];

% tiledlayout('flow','TileSpacing','compact','Padding','compact')
nIDs = length(ID);
switch nIDs
    case 1
        sz = [1 1];
    case 2
        sz = [2 1];
    case 3
        sz = [3 1];
    case 4
        sz = [2 2];
    otherwise
        sz = [];
end
if ~isempty(sz)
    paperfigure(sz(1),sz(2));
else
    tiledlayout('flow','TileSpacing','compact','Padding','compact');
end

for i = 1:nIDs
    nexttile
    ptmp = parity{i};
    if NV.charlblQ
        alphabet = ('a':'z').';
        chars = num2cell(alphabet(1:nIDs));
        chars = chars.';
        charlbl = strcat('(',chars,')');
    else
        charlbl = repelem({''},length(ID));
    end
    if NV.titleQ
        t = ID;
    else
        t = repelem({''},length(ID));
    end
    opts = struct('charlbl',charlbl{i},'title',t(i),'scatterOpts',...
        struct('MarkerEdgeAlpha',0.1),'xlim',NV.xlim,'ylim',NV.ylim);
    if isempty(NV.xlim)
        opts = rmfield(opts,'xlim');
    end
    if isempty(NV.ylim)
        opts = rmfield(opts,'ylim');
    end
    optpairs = namedargs2cell(opts);
    parityplot(ptmp.ytrue,ptmp.ypred,plottype,optpairs{:})
end