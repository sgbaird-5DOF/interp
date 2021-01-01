% RANDOCTPARITYDATA  submit sets of jobs to supercomputer or run locally for interp5DOF()
clear; close all

T = true;
F = false;

%% parameters
%loop through different combinations of parameters using random,
%octochorically sampled octonions
addpathdir({'var_names.m','writeparfile.m','walltimefns'})
runtype = 'full'; %'test','full'
nreps = 10; % number of runs or repetitions

% job submission environment
env = 'slurm'; %'slurm', 'local'
%whether to skip running the jobs and just compile results
dryrunQ = F;
metaQ = F; %whether to load full model or only metay-data at end

%make sure the parameters here correspond with the input to "pars" below
switch runtype
    case 'test'
        ninputpts = floor(58604*0.8); %56442; %floor(67886*0.8); %floor(264276*.8); %17176; %1893*2; %[2366]; %[1893*1]; % 5000 10000 20000 50000];
        npredpts = ceil(58604*0.2); %11443; %floor(67886*0.2); %ceil(264276*0.2); %67886-17176; %67886-1893*2; %65520; %473*1;
        method = {'gpr'}; % 'sphbary', 'pbary', 'gpr', 'sphgpr', 'nn', 'avg'
        datatype = {'kim'}; % 'brk', 'kim', 'rohrer-Ni', 'rohrer-test', 'rohrer-brk-test'
        pgnum = 32; %m-3m (i.e. m\overbar{3}m) FCC symmetry default for e.g. Ni
        sigma = [0]; %J/m^2, standard deviation, added to "y"
        genseed = 11;
        brkQ = false;
        
    case 'full'
        ninputpts = [100 388 500 1000 5000 10000 20000 50000]; % 388, 500, 1000, 2000, 5000, 10000, 20000, 50000
        npredpts = 10000;
        method = {'sphgpr','gpr','sphbary','pbary','nn','avg','idw'}; % 'sphbary', 'pbary', 'gpr', 'sphgpr', 'nn', 'avg'
        datatype = {'brk'};
        pgnum = 32; %m-3m (i.e. m\overbar{3}m) FCC symmetry default for e.g. Ni
        sigma = [0]; %mJ/m^2, standard deviation, added to "y"
        genseed = []; %set to empty if you want to vary for each run
        brkQ = false;
end

%comment (no spaces, used in filename)
% comment = 'paper-data2';
% comment = 'rohrer-Ni-test5';
% comment = 'rohrer-Ni-regularization';
% comment = 'rohrer-Ni-lamb500m';
% comment = 'rohrer-Ni-lamb300m';
% comment = 'rohrer-Ni-lamb300m';
% comment = 'rohrer-test';
% comment = 'kim-test5-neg-epsijk';
% comment = 'kim-test6-pos-epsijk';
% comment = 'kim-paper-data2';
% comment = 'kim-paper-data-exact';
% comment ='kim-mech-spec-equal-ard';
% comment = 'kim-mech-spec-equal-hyperopt';
% comment = 'kim-mech-spec-equal-exact-exact';
% comment = 'kim-mech-train-spec-test-fic-fic';
% comment = 'kim-mech-fic-fic';
% comment = 'kim-brk-minrepeats2';
% comment = 'kim-minrepeats';
% comment = 'kim-minrepeats3';
% comment = 'kim-minrepeats2-trainsigma0.2';
% comment = 'kim-minrepeats2-traintestsigma0.2';
% comment = 'kim-minrepeats5-trainsigma0.2-posnoise';
% comment = 'kim-rng11';
% comment = 'paper-data3';
comment = 'test';
% comment = 'kim-trainRepeat-testNoRepeat';
% comment = 'rohrer';
% comment = 'kim-paper-data-equal-sig0.005-exact-exact';
% comment = 'kim-brk';
% comment = 'kim-brk-sig-0.1-test2';
% comment = 'kim-brk-sig-0.1-exp-cstsig';
% comment = 'test';
% comment = 'brk-test-kim-compare-sig-0';
% comment = 'kim-brk-sig-0-fic';
% comment = 'brk-test';
% comment = 'rohrer-brk-test';
% comment = 'brk-test-uniform';
% comment = 'idw-test-3pt5deg';

disp(['env = ' env])

if strcmp(env,'slurm') && dryrunQ
    error('did you mean to change dryrunQ to false?')
end

disp(['dryrunQ = ' int2str(dryrunQ)])
if strcmp(env,'local')
    savecatQ = T; % whether to save the catenated model and/or parameters (depends on metaQ)
    disp(['savecatQ = ' int2str(savecatQ)])
    disp(['metaQ = ' int2str(metaQ)])
end

m = input(['default comment: ' comment '. Continue (y) or override (n)? '],'s');
if ~strcmp(m,'y') && ~strcmp(m,'Y')
    disp('enter a comment without spaces for use in a filename')
    comment = input('new comment: ','s');
end

% # cores
switch env
    case 'slurm'
        cores = 12;
    case 'local'
        if ~dryrunQ
            p = gcp;
            cores = p.NumWorkers;
        else
            cores = 1;
        end
end

%% functions to generate save filepaths
%diary
files = dir(fullfile('**','data','randOctParity','diary'));
diaryfolder = files(1).folder;
diarynamefn = @(method,ninputpts,gitcommit,puuid) [method int2str(ninputpts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '_' comment '.txt'];
diarypathfn = @(method,ninputpts,gitcommit,puuid) fullfile(diaryfolder,diarynamefn(method,ninputpts,gitcommit,puuid));
%data
files = dir(fullfile('**','data','randOctParity','pcombs'));
savefolder = files(1).folder;
savenamefn = @(method,ninputpts,gitcommit,puuid) [method int2str(ninputpts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '_' comment '.mat'];

%for use with dir
savepathgen = fullfile(savefolder,'*gitID-*puuID*.mat');

savenamematch = [ ...
    ['(' strjoin(method,'|') ')'] ... match (exactly) any of the method options
    ['*'] ... %     ['(' strjoin(cellfun(@num2str,num2cell(ninputpts),'UniformOutput',false),'|') ')'] ... match (exactly) any of the ninputpts options
    '(_gitID-[a-z0-9]*)' ... match any combination of 0 or more lowercase alphabetic or numeric characters (for git hash)
    '(_puuID-[a-z0-9]+)' ... match any combination of 1 or more lowercase alphabetic or numeric characters (for param combo uuid)
    ['(_' comment ')'] ...
    ]; %e.g. '(sphbary|gpr|nn)(50)(_gitID-)[a-z0-9]*)(_puuID-[a-z0-9]+)' matches sphbary50_gitID-abcd123_puuID-abcd123.mat
if metaQ
    savenamematch = [ savenamematch '(_meta.mat)'];
else
    savenamematch = [savenamematch '(.mat)'];
end
    
savepathfn = @(method,ninputpts,gitcommit,puuid) fullfile(savefolder,savenamefn(method,ninputpts,gitcommit,puuid));

%parameters
%**ADD ALL PARAMETERS HERE** (see runtype switch statement)
pars = var_names(ninputpts,npredpts,method,cores,datatype,pgnum,sigma,genseed,brkQ);
if ~dryrunQ
    %% parameter file setup
    %function to execute and output arguments from function
    execfn = @(ninputpts,npredpts,method,datatype,pgnum,sigma,genseed,brkQ) ... **NAMES NEED TO MATCH PARS FIELDS** (see above)
        interp5DOF_setup(ninputpts,npredpts,method,datatype,...
        'pgnum',pgnum,'sigma',sigma,'genseed',genseed,'brkQ',brkQ); %**NAMES NEED TO MATCH PARS FIELDS AND EXECFN ARGUMENTS**
    argoutnames = {'ypred','interpfn','mdl','mdlpars'};
    %i.e. [ypred,interpfn,mdl,mdlpars] = interp5DOF_setup(ninputpts,npredpts,method,datatype,...);
    
    % walltimefn = @() 300; %can set to constant or to depend on parameters, probably fine when using standby queue
    walltimefn = @(ninputpts,npredpts,method,cores,datatype) get_walltimefn(ninputpts,npredpts,method,cores,datatype);
    
    %% parameter file
    [parpath, parcombsets, Ntrim, jobwalltimes] = ...
        writeparfile(pars,execfn,argoutnames,walltimefn,'diarypathfn',diarypathfn,'savepathfn',savepathfn,'nreps',nreps);
end
%% job submission
switch env
    case 'slurm'
        %setup
        mem = 1024*12*cores; %total memory of job, MB
        qosopt = 'standby'; %'', 'test', 'standby'
        scriptfpath = fullfile('MATslurm','code','submit.sh');
        %submission
        submit_sbatch(parpath,cores,mem,qosopt,scriptfpath); %submit.sh --> exec_combs.m --> execfn (defined above)
        
    case 'local'
        if ~dryrunQ
            %nested loop through jobs and tasks to generate results
            for jid = 1:length(Ntrim)
                N = Ntrim(jid);
                for tid = 1:N
                    outtmp = exec_combs(parpath, jid, tid);
                end
            end
        end
        %get names and folders from files that match general savepath
        files = dir(savepathgen);
        namestmp = {files.name};
        folderstmp = {files.folder};
        
        %get names and folders of filenames that match savenamematch
        startIDs = regexp(namestmp,savenamematch);
        matchIDs = ~cellfun(@isempty,startIDs);
        names = namestmp(matchIDs);
        folders = folderstmp(matchIDs);
        nfiles = length(names);
        assert(nfiles ~= 0,'no files matched. Verify files exist and check regexp savenamematch.')
        %initialize
        init1 = cell(nfiles,1);
        [Slist,ypredlist,mdllist,mdlparslist,interpfnlist] = deal(init1);
        
        %extract results from files
        for i = 1:nfiles
            if mod(i,50) == 0
                disp(i)
            end
            folder = folders{i};
            name = names{i};
            fpath = fullfile(folder,name);
            
            if metaQ
                S = load(fpath);
                mdlparslist{i} = S;
            else
                loadvars = {'ypred','mdl','mdlpars','interpfn'};
                S = load(fpath,loadvars{:});
                % if you get a "variable not loaded here" warning, you may need to
                % delete a previous erroneous data file you generated with the
                % wrong variable names. Typically shouldn't be an issue
                % though..
                [ypredlist{i},mdllist{i},mdlparslist{i},interpfnlist{i}] = ...
                deal(S.ypred, S.mdl, S.mdlpars, S.interpfn);
            end
            %             Slist{i} = S;
        end
        
        if ~metaQ
            %concatenate models and parameters
            mdlcat = structvertcat(mdllist{:});
            mdltbl = struct2table(mdlcat,'AsArray',true);
            mdltbl = tblfilt(mdltbl,pars);
            mdlcat = table2struct(mdltbl);
            clear mdllist
            %         mdltbl = struct2table(mdlcat,'AsArray',true);
        end
        
        mdlparscat = structvertcat(mdlparslist{:});
        mdlparstbl = struct2table(mdlparscat,'AsArray',true);
        
        mdlparstbltmp = tblfilt(mdlparstbl,pars);
        if isempty(mdlparstbltmp)
            error('mdlparstbltmp was was empty, check tblfilt()')
        else
            mdlparstbl = mdlparstbltmp;
        end
        mdlparscat = table2struct(mdlparstbl);
        
        %         gitcommit = get_gitcommit();
        %save models and parameters
        fpath = fullfile(savefolder,['gitID-' get_gitcommit '_uuID-' get_uuid() '_' comment '.mat']);
        if savecatQ
            if metaQ
                savevars = {'mdlparscat','mdlparstbl'};
            else
                savevars = {'ypredlist','interpfnlist','mdlcat',...
                    'mdlparscat','mdlparstbl'};
            end
            save(fpath,savevars{:},'-v7.3')
        end
        writetable(mdlparstbl,[fpath(1:end-4) '.csv'],'WriteVariableNames',true)
        disp(mdlparstbl)
        
        %         %nested loop through jobs and tasks to load results
        %         S = load(parpath);
        %         for jid = 1:length(Ntrim)
        %             for tid = 1:jid
        %                 S.
        %             end
        %         end
        mdltbl(:,{'method','ninputpts','npredpts','rmse','mae'})

        %% plotting
        mdlnum = 1;
        mdl = mdlcat(mdlnum);
        Kpars = mdl.KernelInformation.KernelParameters;
        Lval = Kpars(1);
        sigval = Kpars(2);
        paperfigure(1,2);
        nexttile
        t1 = ['$L_{\mathrm{kernel}}$: ' num2str(rad2deg(Lval)*2) ' ($^\circ{}$), ' ...
        '$\sigma_{\mathrm{kernel}}$: ' num2str(sigval) ' ($J m^{-2}$)'];
        brkQ = true;
        tunnelplot(mdl,'brkQ',brkQ)
        title(t1,'Interpreter','latex')
        % multiparity({mdl.errmetrics},'charlblQ',false);
        nexttile
        % parityplot(mdl.errmetrics.ytrue,mdl.errmetrics.ypred,'scatter',...
        %     'mkr','o','fillQ',true,'scatterOpts',struct('MarkerFaceAlpha',0.005));

        parityplot(mdl.errmetrics.ytrue,mdl.errmetrics.ypred);

        t2 = ['RMSE: ' num2str(mdl.rmse) ' ($J m^{-2}$), MAE: ' num2str(mdl.mae) ' ($J m^{-2}$)'];
        title(t2,'Interpreter','latex')

        t3 = ['method: ' upper(char(mdl.method)) ', datatype: ' upper(char(mdl.datatype)),...
        ', ninputpts: ' num2str(mdl.ninputpts) ', npredpts: ' num2str(mdl.npredpts)];
        sgtitle(t3,'Interpreter','latex')
end

disp('end randOctParityData.m')
disp(' ')


%% CODE GRAVEYARD
%{

% sigma plotting
[~,ids] = sort([mdlparscat.sigma]);
mdlparscat = mdlparscat(ids);
titlelist = strcat('$$\sigma_y$$ =',{' '},cellfun(@num2str,num2cell(vertcat(mdlparscat.sigma)),'UniformOutput',false),' $$J/m^2$$')
multiparity({mdlparscat.errmetrics},[2 3 4 1],'titlelist',titlelist)


for i = 1:length(ninputptsList)
    %unpack
    ninputpts = ninputptsList(i);
    for j = 1:length(methods)
        method = methods{j};
        fname = ['randOctParityData_' method int2str(ninputpts) '.mat'];
        S = load(
    end
end

% execfn = @interp5DOF_setup;

%initialize
init1 = cell(length(ninputpts),length(method));
propOutlist = init1;
mdllist = init1;
mdlparslist = init1;
interpfnlist = init1;


% %nested loop through ninputptsList and methodlist
% for i = 1:length(ninputpts)
%     %unpack
%     ninputpts = ninputpts(i);
%     for j = 1:length(method)
%         %% setup
%         %unpack
%         method = method{j};
%         disp(method)
%
%         %get a unique ID
%         uuid = get_uuid();
%
%         %% interpolation
%         [propOutlist{i,j},interpfnlist{i,j},mdllist{i,j},mdlparslist{i,j}]=interp5DOF_setup(ninputpts,npredpts,method,uuid);
%
%         %% output
%         %package
%         propOut = propOutlist{i,j};
%         interpfn = interpfnlist{i,j};
%         mdl = mdllist{i,j};
%         mdlpars = mdlparslist{i,j};
%         %save
%         save(savepathfn(method,ninputpts,mdl.gitcommit,uuid),...
%             'propOut','interpfn','mdl','mdlpars','ninputpts','method','uuID',...
%             '-v7.3','-nocompression')
%         disp(' ')
%     end
% end

%git commit
[status,cmdout] = system('git log -1 --oneline');
if status == 0
    gitcommit = cmdout;
else
    gitcommit = '';
end

%job universal unique ID
uuid = {get_uuid()};

pars = var_names(ninputpts,npredpts,method.uuid);

%         S(nfiles) = struct();
%         S(1).out = struct();


memfn = @(method,ninputpts) ...
    (strcmp(method,'gpr') && ninputpts > 10000) *1024*8 + ...
    (~strcmp(method,'gpr') && ninputpts <= 10000) * 1024*4;

%         k = 0; %initialize counter
        %nested loop through jobs and tasks to generate results
        for jid = 1:length(Ntrim)
            for tid = 1:jid
%                 k = k+1;
                outtmp = exec_combs(parpath, jid, tid);
            end
        end

diaryfolder = fullfile('data','randOctParity','diary');
savefolder = fullfile('data','randOctParity','pcombs');


% job unique identifier
% juuid = get_uuid();


        emptyIDs = cellfun(@isempty,mdlparslist);
        if metaQ
            [mdlparslist{emptyIDs}] = deal([]);
        else
            [ypredlist{emptyIDs},mdllist{emptyIDs},mdlparslist{emptyIDs},interpfnlist{emptyIDs}] = ...
            deal([]);
        end


                if any(strcmp(S.datatype,datatype))
                    mdlparslist{i} = S;
                end


                if any(strcmp(S.mdl.datatype,datatype))
                    [ypredlist{i},mdllist{i},mdlparslist{i},interpfnlist{i}] = ...
                    deal(S.ypred, S.mdl, S.mdlpars, S.interpfn);
                end

%}
