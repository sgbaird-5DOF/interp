function multixyplots(tbl,xytypes,xtype,ytypes,nrows,ncols,NV)
arguments
    tbl
    xytypes cell = {'pbary','gpr','idw','nn'}
    xtype char = 'nmeshpts'
    ytypes(1,:) cell = {'rmse','mae'}
    nrows(1,1) double = []
    ncols(1,1) double = []
    NV.yunits char = 'J/m^2'
    NV.XScale char {mustBeMember(NV.XScale,{'log','linear'})} = 'log'
    NV.YScale char {mustBeMember(NV.YScale,{'log','linear'})} = 'linear'
    NV.xmin double = []
    NV.ymin double = []
end

ntypes = length(ytypes);

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
    tiledlayout(nrows,ncols);
else
    tiledlayout;
end

for i = 1:ntypes
    ytype = ytypes{i};
    ax1 = nexttile;
    hold(ax1,'on')
    ax1.XScale = NV.XScale;
    ax1.YScale = NV.YScale;
    NVpairs = namedargs2cell(NV);
    xyplots(tbl,xytypes,xtype,ytype,NVpairs{:})
end

end