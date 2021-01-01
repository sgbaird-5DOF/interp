function xyplots(mastertbl,xytypes,xtype,ytype,NV)
arguments
   mastertbl table
   xytypes cell
   xtype char
   ytype char
   NV.xtypelbl char = xtype
   NV.ytypelbl char = ytype
   NV.xytypelbl char = xytypes
   NV.yunits char = 'J/m^2'
   NV.XScale char {mustBeMember(NV.XScale,{'log','linear'})} = 'log'
   NV.YScale char {mustBeMember(NV.YScale,{'log','linear'})} = 'linear'
   NV.xmin double = []
   NV.ymin double = []
   NV.xmax double = []
   NV.ymax double = []
   NV.lgdloc char = 'best'
   NV.charlblQ(1,1) logical = false
   NV.charlblnum double = []
   NV.Interpreter char {mustBeMember(NV.Interpreter,{'latex','tex'})} = 'latex'
end
% XYPLOTS  plot multiple datasets on the same axes using a "master" table and findgroups()
ax = gca;
ax.XScale = NV.XScale;
ax.YScale = NV.YScale;
xlabel(NV.xtypelbl)
if ~isempty(NV.yunits)
    ylabel([NV.ytypelbl ' (' NV.yunits ')'],'Interpreter',NV.Interpreter)
else
    ylabel(NV.ytypelbl,'Interpreter',NV.Interpreter)
end
axis square

for xytype = xytypes
    tbl = mastertbl(ismember(mastertbl.method,xytype),:);
    [G,ID] = findgroups(tbl.ninputpts);
    
    xyplot(tbl,G,xtype,ytype)
end

if ~isempty(NV.xmin)
    ax.XLim(1) = NV.xmin;
end
if ~isempty(NV.ymin)
    ax.YLim(1) = NV.ymin;
end
if ~isempty(NV.xmax)
    ax.XLim(1) = NV.xmax;
end
if ~isempty(NV.ymax)
    ax.YLim(2) = NV.ymax;
end

legend(NV.xytypelbl,'Location',NV.lgdloc,'FontSize',9,'Interpreter',NV.Interpreter)

% label for figure tiles, e.g. '(a)', '(b)', '(c)', '(d)'
if NV.charlblQ
    papertext(NV.charlblnum)
%     text(0.025,0.95,NV.charlbl,'Units','normalized','FontSize',12,'Interpreter',NV.Interpreter)
end

end

%% CODE GRAVEYARD
%{
%     [npts,rmse,mae,runtime] = deal(tbl.ninputpts,tbl.rmse,tbl.mae,tbl.runtime);
    
    ninputpts = splitapply(@(x)x(1),npts,G1);
    medRMSE = splitapply(@median,rmse,G1);
    stdRMSE = splitapply(@std,rmse,G1);
    medMAE = splitapply(@median,mae,G1);
    stdMAE = splitapply(@std,mae,G1);
    medruntime = splitapply(@median,runtime,G1);
    stdruntime = splitapply(@std,runtime,G1);
    
    semilogx(t1,ninputpts,medRMSE,'-o')
    xlabel('ninputpts')
    ylabel('RMSE (J/m^2)')
    axis square
    hold off
    
    semilogx(t2,ninputpts,medMAE,'-o')
    xlabel('ninputpts')
    ylabel('MAE (J/m^2)')
    axis square
    hold off
%     semilogx(ninputpts,medruntime,'-o')

    figure(fig2)
    hold on
    semilogx(ninputpts,medruntime,'-o')
    xlabel('ninputpts')
    ylabel('runtime (s)')
    axis square
    hold off

%   NV.charlbl char = ''
%}