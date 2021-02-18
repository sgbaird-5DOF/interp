% Zhang Fe Mobility Data Test
% Using data from Zhang, J.; Ludwig, W.; Zhang, Y.; Sørensen, H. H. B.;
% Rowenhorst, D. J.; Yamanaka, A.; Voorhees, P. W.; Poulsen, H. F. Grain
% Boundary Mobilities in Polycrystals. Acta Materialia 2020, 191, 211–220.
% https://doi.org/10.1016/j.actamat.2020.03.044.

clear; close all
%% Setup
initialepsijk = 1; %active (1) passive (-1)
files = dir(fullfile('**','interp5DOF-paper','figures'));
figfolder = files(1).folder;

%% Import
A = importdata('zhang-fe-mobility\mobility_5d.csv');
data = A.data;
eu = data(:,1:3);
nA = data(:,4:6);
y = data(:,10);

%% Property Histograms
setlatex()
xypos=[0.9,0.95];
paperfigure(1,2);
i = 0;
i = i+1;
nexttile(i)
histogram(y); axis tight
papertext(i,'xypos',xypos)
xylabel()
i = i+1;
nexttile(i);
histogram(y(y<0.4)); axis tight
papertext(i,'xypos',xypos)
xylabel()

savefigpng(figfolder,'zhang-mobility-hist')

%% Convert to Octonions
qm = eu2qu(eu,initialepsijk);
o = five2oct(qm,nA,epsijk);
o = get_octpairs(o);
[X,usv] = proj_down(o);

%% fit GPR Model (LOOCV)
gprMdl = fitrgp(X,y,'Leaveout','on');
% gprMdl = fitrgp(X,log(y+1),'Leaveout','on');

%% fit Lower Mobility Values
Xsub = X(y<0.4,:);
ysub = y(y<0.4);
gprMdlSub = fitrgp(Xsub,ysub,'Leaveout','on');

%% Results
ypred = kfoldPredict(gprMdl);
paperfigure(1,2,6.729);
i = 0;
i = i+1;
nexttile(i)
parityplot(y,ypred);
papertext(i)
i=i+1;
nexttile(i)
ypredsub = kfoldPredict(gprMdlSub);
parityplot(ysub,ypredsub);
papertext(i)
% parityplot(y,ypred,'xlim',[0,0.4],'ylim',[0,0.4]);

%% Helper Function(s)
function xylabel()
xlabel('Reduced Mobility ($\mu m^2/s$)')
ylabel('Number of GBs')
end
