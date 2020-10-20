function xyplots(mdlparstbl,xytypes,xtype,ytype,NV)
arguments
   mdlparstbl table
   xytypes cell
   xtype char
   ytype char
   NV.yunits char = 'J/m^2'
   NV.XScale char {mustBeMember(NV.XScale,{'log','linear'})} = 'log'
   NV.YScale char {mustBeMember(NV.YScale,{'log','linear'})} = 'linear'
   NV.xmin double = []
   NV.ymin double = []
end

ax = gca;
ax.XScale = 'log';
xlabel(xtype)
if ~isempty(NV.yunits)
    ylabel([ytype ' (' NV.yunits ')'])
else
    ylabel(ytype)
end
axis square

for xytype = xytypes
    tbl = mdlparstbl(ismember(mdlparstbl.method,xytype),:);
    [G,ID] = findgroups(tbl.nmeshpts);
    
    xyplot(tbl,G,xtype,ytype)
end

if ~isempty(NV.xmin)
    ax.XLim(1) = NV.xmin;
end
if ~isempty(NV.ymin)
    ax.YLim(1) = NV.ymin;
end

legend(xytypes,'Location','best')

end

%% CODE GRAVEYARD
%{
%     [npts,rmse,mae,runtime] = deal(tbl.nmeshpts,tbl.rmse,tbl.mae,tbl.runtime);
    
    nmeshpts = splitapply(@(x)x(1),npts,G1);
    medRMSE = splitapply(@median,rmse,G1);
    stdRMSE = splitapply(@std,rmse,G1);
    medMAE = splitapply(@median,mae,G1);
    stdMAE = splitapply(@std,mae,G1);
    medruntime = splitapply(@median,runtime,G1);
    stdruntime = splitapply(@std,runtime,G1);
    
    semilogx(t1,nmeshpts,medRMSE,'-o')
    xlabel('nmeshpts')
    ylabel('RMSE (J/m^2)')
    axis square
    hold off
    
    semilogx(t2,nmeshpts,medMAE,'-o')
    xlabel('nmeshpts')
    ylabel('MAE (J/m^2)')
    axis square
    hold off
%     semilogx(nmeshpts,medruntime,'-o')

    figure(fig2)
    hold on
    semilogx(nmeshpts,medruntime,'-o')
    xlabel('nmeshpts')
    ylabel('runtime (s)')
    axis square
    hold off
%}