function multiparity(parity,ID,plottype,NV)
arguments
    parity(:,1) cell
    ID
    plottype char {mustBeMember(plottype,{'hex','scatter'})} = 'hex'
    NV.charlblQ(1,1) logical = true
    NV.autotitle(1,1) logical = false
    NV.titlelist = []
    NV.xlim = []
    NV.ylim = []
    NV.cbnds = []
    NV.Interpreter char {mustBeMember(NV.Interpreter,{'tex','latex','none'})} = 'latex'
end
% MULTIPARITY  create tiled parity plots using cell parity data

if ~isempty(NV.titlelist)
    assert(~NV.autotitle,'if titlelist is specified, autotitle should be set to false')
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
    switch sz(2)
        case 1
            ht = 14.5838333333333/2;
        case 2
            ht = 14.5838333333333;
        otherwise
            ht = [];
    end
    paperfigure(sz(1),sz(2),ht);
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
    if ~isempty(NV.titlelist)
        t = NV.titlelist;
    else
        if NV.autotitle
            t = ID;
        else
            t = repelem({''},length(ID));
        end
    end
    opts = struct('charlbl',charlbl{i},'title',t(i),'scatterOpts',...
        struct('MarkerEdgeAlpha',0.1),'xlim',NV.xlim,'ylim',NV.ylim,'cbnds',NV.cbnds);
    if isempty(NV.xlim)
        opts = rmfield(opts,'xlim');
    end
    if isempty(NV.ylim)
        opts = rmfield(opts,'ylim');
    end
    optpairs = namedargs2cell(opts);
    parityplot(ptmp.ytrue,ptmp.ypred,plottype,optpairs{:})
end