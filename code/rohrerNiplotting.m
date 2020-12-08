%saving
files = dir(fullfile('**','interp5DOF-paper','figures'));
folder = files(1).folder;

[~,ids] = sort([mdlparscat.sigma]);
mdlparscat = mdlparscat(ids);
titlelist = strcat('$$\sigma_y$$ =',{' '},cellfun(@num2str,num2cell(vertcat(mdlparscat.sigma)),'UniformOutput',false),' $$J/m^2$$')

multiparity({mdlparscat.errmetrics},{mdlparscat.method},'autotitle',true)
multiparity({mdlparscat.errmetrics},titlelist,'autotitle',true)

mdlparstbl(:,{'datatype','sigma','method','starttime','rmse','mae'})

%%
paperfigure();
[~,D] = get_knn(mdl.data.pts,'omega',1);
histogram(D{1})
[~,lonelyID] = max(D{1});
pts = mdl.mesh.pts;
[~,D2] = get_knn(mdl.data.pts,'omega',10,'ID',lonelyID);
plot(1:10,[D2{:}],'o')
xlabel('n-th NN of max NN distance point','Interpreter','latex')
ylabel('$\omega (^{\circ})$','Interpreter','latex')

savefigpng(folder,'lonely-knn');

%% standard deviation color x-y plot and NN distance color x-y plot
% ids = randi(10000,10000,1);
ids = 1:10000;
paperfigure(2,2);
nexttile
ysd = mdl.ysd;
parityplot(mdl.errmetrics.ytrue(ids),mdl.errmetrics.ypred(ids));
nexttile
ysd2 = -rescale(ysd,-50,-5);
scatter(mdl.errmetrics.ytrue(ids),mdl.errmetrics.ypred(ids),ysd2(ids),ysd(ids),'o');
axis square
xlims = xlim;
ylims = ylim;
hold on
refline(1,0)
xlim(xlims);
ylim(ylims);
colorbar
colormap jet

nexttile
d = D{1};
d2 = rescale(d,5,50);
scatter(mdl.errmetrics.ytrue(ids),mdl.errmetrics.ypred(ids),d2(ids),d(ids),'.');
% set(gca,'ColorScale','log')
axis tight
xlims = xlim;
ylims = ylim;
hold on
refline(1,0)
xlim(xlims);
ylim(ylims);
axis square

colorbar
colormap jet

nexttile
xlims = [min(d(ids)),max(d(ids))];
ylims = [min(ysd(ids)),max(ysd(ids))];
hexscatter(d(ids),ysd(ids),xlims,ylims,'axis','normal');
xlabel('$\omega$ NN distance $(^\circ{})$','Interpreter','latex')
ylabel('$\sigma_{y} (J/m^{-2}$)','Interpreter','latex')
axis square
savefigpng(folder,'parity-ysd');

%%
paperfigure();
histogram(d)
xlabel('$\omega (^{\circ})$','Interpreter','latex')
ylabel('counts','Interpreter','latex')

savefigpng(folder,'rohrer-nn-hist');


%%
paperfigure();
plot([5000 10000 15000 20000],[76.73,110.32,157.67 189.39],'-o')
xlabel('total \# datapoints','Interpreter','latex')
ylabel('time (s)','Interpreter','latex')
title('PredictMethod=sd, ActiveSetMethod=sgma','Interpreter','latex')
savefigpng(folder,'sgma-linear-time')
%active set size is 2000
%estimated time for 180000 pts is 1704.5 s

%%
paperfigure()
mdl = mdlcat(1);
ids = randperm(10000,10000);
% ids2 = setdiff(1:10000,ids);
ypred = mdl.errmetrics.ypred;
ytrue = mdl.errmetrics.ytrue;
ytrue(ids) = normrnd(ytrue(ids),0.1);
parityplot(ytrue,ypred)

%%
paperfigure()
tunnelplot()
%saving
files = dir(fullfile('**','interp5DOF-paper','figures'));
folder = files(1).folder;
savefigpng(folder,'tunnelplot-lamb300m')


