% posterior_plotting
%% setup
% files = dir(fullfile('**','interp5DOF-paper','figures'));
% figfolder = files(1).folder;

addpath(genpath('.'))

set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultAxesFontSize',12)

figfolder = 'C:\Users\sterg\Documents\GitHub\posterior-sampling\figures';

brkstr = '\\acrfull{brk} validation function \cite{bulatovGrainBoundaryEnergy2014} (black)';
festr = '\ch{Fe} simulation dataset \cite{kimPhasefieldModeling3D2014}';

str1D = ['1D interpolation results for %s %s. The %s is used. The line \\overline{AB} is linearly extended '...
'\\SI{%d}{\\percent} in either direction before renormalizing and taking a uniform sampling of %d points.'];
%example: caption = sprintf(str1D,captionlist,'input points',brkstr,extend*100,n);

%% test function truncated multi-variate normal distribution sampling
testPosteriorCovarianceGP()
fname = 'testfn-mvn-tmvn';
savefigpng(figfolder,fname);

%% extended 1D arc, no ensemble or mixture, BRK validation
testnum = 2;
n = 100;
extend = 0.5; %percent to extend in either direction
ninputptslist = [1000 50000]; %1000 or 50000
paperfigure(2,2);

% load('tunnel-50000','A','B')
data = importdata(fullfile(figfolder,'AB-50000.csv')); %two points that are far apart in a particular VFZ
A = data(1,:);
B = data(2,:);

i = 0;
for ninputpts = ninputptslist
    nexttile
    i=i+1;
    tunnelplot_test(testnum,ninputpts,n,[],[],[],[],A,B,'extend',extend);
    papertext(i,'xypos',[0.1,0.95])
end
lbls = lblcat(ninputptslist);
savename = ['extended1Darc-gpr-brk-',lbls];
savefigpng(figfolder,savename)

captionlist = get_captionlist(ninputptslist);

caption = sprintf(str1D,captionlist,'input points',brkstr,extend*100,n);

savefigstr(caption,savename,figfolder);

%% Ensemble VFZO
for n = [1000 10000]
    seed = 10;
    rng(seed);
    n2 = min([n 10000]);
    K = 10;
    method = 'gpr';
    [ypred,ytrue,ypredlist,interpfnlist,mdllist,mdlparslist] = ensembleVFZO_test(n,n2,K,method);
    savename = ['ensemble-interp-',int2str(n)];
    savepath = fullfile(figfolder,savename);
    save(savepath,'mdllist','ypredlist')
end

%% four-panel extended arc (GPR, GPRM, EGPR, EGPRM)
seed = 10;
rng(seed);
n = 10000;
n2 = 10000;
thr = 1.1;
scl = 30;
epsijk = 1;

savename = ['panel-arc-' int2str(n)];

loadQ = true;
if loadQ
    load(savename,'mdls','ntunnelpts','extend','A','B')
else
    mdls = cell(1,4);
end

%%
paperfigure(2,2);
i = 0;
for K = [10 1]
    for mixQ = [true false]
        i = i+1;
        nexttile
        [mdls{i},ntunnelpts,extend,A,B] = ...
            egprm_test(n,n2,K,thr,scl,epsijk,'mixQ',mixQ,'lgdloc','south','mdl',mdls{i});
        switch n
            case 1000
                ylim([0.8,1.425])
            case 10000
                ylim([0.1,1.4])
        end
        papertext(i)
    end
end
%%
saveQ = false;
if saveQ
    save(savename,'mdls','ntunnelpts','extend','A','B','-v7.3') %#ok<UNRCH>
end
savefigpng(figfolder,savename);
% set(gca,'Legend','off')
%% 
captionlist = get_captionlist({...
    '10-component \\acrfull{egprm}',...
    '10-component \\acrfull{egpr}',...
    '\\acrfull{gprm}',...
    '\\acrfull{gpr}'});
str1Dci = [str1D '. \\SI{95}{\percent} confidence intervals (CI) are also shown.'];
caption = sprintf(str1D,captionlist,'models',brkstr,extend*100,ntunnelpts);
savefigstr(caption,savename,figfolder);


%% 2D visualization
% make a new function

% pick three points
% project to lower dimension (proj_down)
% get subdivided points (facet_subdiv)
% convert back to octonions (proj_up)
% sample new GBE (GB5DOF_setup)
% plot in ternary diagram (ternplot, add as submodule from GitHub)

%% 3D visualization
% make a new function 

% same steps as 2D visualization, except pick 4 points and plot in 3D

%% CODE GRAVEYARD
%{
nexttile
mdls{2} = egprm_test(n,n2,10,thr,scl,epsijk,'mixQ',false);
ylim([0.8,1.4])
nexttile
mdls{3} = egprm_test(n,n2,1,thr,scl,epsijk,'mixQ',true);
ylim([0.8,1.4])
nexttile
mdls{4} = egprm_test(n,n2,1,thr,scl,epsijk,'mixQ',false);
ylim([0.8,1.4])

% strArc = ['1D interpolation results for %s. Thresholds of %d and %d and ' ...
%     'a scale of %d were used for the mixture models.'];
% caption = sprintf(strArc,captionlist,'\\acrfull{brk} validation function',n,extend*100);

% caption = ['1D interpolation results for ' captionlist ' input points. ' ...
% 'The \acrfull{brk} validation function was used to sample ' int2str(n) ...
% ' equally spaced points. The endpoints A and B were extended ' ...
% '\SI{' num2str(extend*100) '}{\percent} in the respective directions.'];
%}