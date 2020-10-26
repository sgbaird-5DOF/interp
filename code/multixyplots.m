function multixyplots(tbl,xytypes,xtype,ytypes,nrows,ncols,NV)
arguments
    tbl
    xytypes cell = {'pbary','gpr','idw','nn'}
    xtype char = 'nmeshpts'
    ytypes(1,:) = {'rmse','mae'}
    nrows(1,1) double = []
    ncols(1,1) double = []
    NV.yunits char = 'J/m^2'
    NV.XScale char {mustBeMember(NV.XScale,{'log','linear'})} = 'log'
    NV.YScale char {mustBeMember(NV.YScale,{'log','linear'})} = 'linear'
    NV.xmin double = []
    NV.ymin double = []
    NV.lgdloc char = 'best'
    NV.charlblQ(1,1) logical = true
end

if iscell(xytypes{1})
    ntypes = length(xytypes);
    loopvar = 'xytypes';
elseif length(ytypes) > 1
    ntypes = length(ytypes);
    loopvar = 'ytypes';
else
    loopvar = 'ytypes';
    ntypes = 1;
end

if ~iscell(ytypes)
    ytypes = {ytypes};
end

fig = figure;
switch ntypes
    case 1
        fig.Position = [460.2 245 498.4 472.8];
    case 2
        if (nrows == 2) && (ncols == 1)
            fig.Position = [489 90.6 330.4 672.4];
        end
end

if ~isempty(nrows) && ~isempty(ncols)
    assert(nrows*ncols >= ntypes,['too many ytypes(' int2str(nrows*ncols) ')/(' int2str(ntypes) ')too few tiles'])
        tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact');
else
    tiledlayout('TileSpacing','compact','Padding','compact');
end

if NV.charlblQ
    alphabet = ('a':'z').';
    chars = num2cell(alphabet(1:ntypes));
    chars = chars.';
    charlbl = strcat('(',chars,')');
else
    charlbl = repelem({''},ntypes,1);
end

for i = 1:ntypes
    switch loopvar
        case 'ytypes'
            ytype = ytypes{i};
        case 'xytypes'
            subxytypes = xytypes{i};
    end
    
    ax1 = nexttile;
    hold(ax1,'on')
    ax1.XScale = NV.XScale;
    ax1.YScale = NV.YScale;
    NV.charlbl = charlbl{i};
    NV.charlblnum = i;
    NVpairs = namedargs2cell(NV);
    switch loopvar
        case 'ytypes'
            xyplots(tbl,xytypes,xtype,ytype,NVpairs{:})
        case 'xytypes'
            xyplots(tbl,subxytypes,xtype,ytypes{1},NVpairs{:})
    end
end

end