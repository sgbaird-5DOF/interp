function xyplot(tbl,G,xtype,ytype)
arguments
    tbl table
    G double
    xtype char = 'nmeshpts'
    ytype char = 'rmse'
end
% XYPLOT  make an errorbar xyplot using a table (tbl) and output (G) from findgroups()

X = tbl.(xtype);
Y = tbl.(ytype);

x = splitapply(@(x)x(1),X,G);

ymean = splitapply(@mean,Y,G);
ystd = splitapply(@std,Y,G);

errorbar(x,ymean,ystd);
end


%% CODE GRAVEYARD
%{
%     nmeshpts = splitapply(@(x)x(1),npts,G1);
%     medRMSE = splitapply(@median,rmse,G1);
%     stdRMSE = splitapply(@std,rmse,G1);
%     medMAE = splitapply(@median,mae,G1);
%     stdMAE = splitapply(@std,mae,G1);
%     medruntime = splitapply(@median,runtime,G1);
%     stdruntime = splitapply(@std,runtime,G1);
    
%     semilogx(t1,nmeshpts,medRMSE,'-o')
%     xlabel('nmeshpts')
%     ylabel('RMSE (J/m^2)')
%     axis square
%     hold off
%
%     semilogx(t2,nmeshpts,medMAE,'-o')
%     xlabel('nmeshpts')
%     ylabel('MAE (J/m^2)')
%     axis square
%     hold off
    %     semilogx(nmeshpts,medruntime,'-o')
    
%     figure(fig2)
%     hold on
%     semilogx(nmeshpts,medruntime,'-o')
%     xlabel('nmeshpts')
%     ylabel('runtime (s)')
%     axis square
%     hold off

%     tbl = mdlparstbl(ismember(mdlparstbl.method,method),:);
%     [G1,ID1] = findgroups(tbl.nmeshpts);
    
%     [npts,rmse,mae,runtime] = deal(tbl.nmeshpts,tbl.rmse,tbl.mae,tbl.runtime);
%}