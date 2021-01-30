%% input parameters
rng(10)
ninputpts = 10000;
npredpts = 10000;
finaldimlist = fliplr(3:8);
nptstot = ninputpts+npredpts;
%% get VFZOs
ofull = get_ocubo(nptstot);
ofull = get_octpairs(ofull);
%% get true values
yfull = GB5DOF_setup(ofull(:,1:4),ofull(:,5:8));
pts = ofull;
%%
paperfigure(3,2,19.6);
i = 0;
for finaldim = finaldimlist
    i = i+1;
    nexttile
    % perform transformations
    ppts = proj_down(pts,'zeroQ',true);
    tripts = projfacet2hyperplane(mean(ppts),ppts);
    pppts = proj_down(tripts,'zeroQ',false);
    switch finaldim
        case 8
            Xfull = pts;
        case 7
            Xfull = ppts;
        case 6
            Xfull = pppts;
        otherwise
            nforce = size(pppts,2)-finaldim;
            downpts = proj_down(pppts,'zeroQ',false,'nforceQ',true,'nforce',nforce);
            Xfull = downpts;
    end
    
    X = Xfull(1:ninputpts,:);
    X2 = Xfull(ninputpts+1:end,:);
    y = yfull(1:ninputpts,:);
    ytrue = yfull(ninputpts+1:end,:);
    
    gprMdl = fitrgp(X,y);
    
    ypred = predict(gprMdl,X2);
    
    errmetrics(i) = get_errmetrics(ypred,ytrue); %#ok<SAGROW>
    parityplot(ytrue,ypred)
    papertext(i)
%     title([int2str(finaldim),'D'],'Interpreter','latex')
end
%%
% files = dir(fullfile('**','interp5DOF-paper','figures'));
% figfolder = files(1).folder;
figfolder = 'C:\Users\sterg\Documents\GitHub\posterior-sampling\figures';
finaldimtxt = arrayfun(@int2str,finaldimlist,'UniformOutput',false);
finaldimlbl = strjoin(finaldimtxt,'-');
fname = ['brkparity',int2str(ninputpts),'-finaldim-',finaldimlbl];
savefigpng(figfolder,fname)

%%
rmse = [errmetrics.rmse];
mae = [errmetrics.mae];
paperfigure();
plot(finaldimlist.',[rmse;mae].','*-')
set(gca,'xDir','reverse')
ticklbl = strcat(fliplr(finaldimtxt),'D');
xticklabels(ticklbl);
set(gca,'TickLabelInterpreter','latex');
xlabel('Final Dimension','Interpreter','latex')
ylabel('Error ($J/m^2$)','Interpreter','latex')
legend('RMSE','MAE','location','best')
fname = ['brk',int2str(ninputpts),'finaldim-',finaldimlbl,'-rmse-mae'];
savefigpng(figfolder,fname)

%%
paperfigure(3,2,10.1);
e = [errmetrics.e];
for i = 1:length(finaldimlist)
    finaldim = finaldimlist(i);
    nexttile
    etmp = e(:,i);
    histogram(etmp)
    axis tight
    xlabel([int2str(finaldim) 'D Error ($J/m^2$)'],'Interpreter','latex')
    ylabel('\# GBs','Interpreter','latex')
    papertext(i,'xypos',[0.025,0.85])
    
end
fname = 'finaldim-errhist';
savefigpng(figfolder,fname)
