clear; close all
%% parameters
%loop through different combinations of parameters using random,
%octochorically sampled octonions
addpathdir({'var_names.m','writeparfile.m','walltimefns'})
runtype = 'full'; %'test','full'
nreps = 10; % number of runs or repetitions

%make sure the parameters here correspond with the input to "pars" below
switch runtype
    case 'test'
        ndatapts = [1000]; % 5000 10000 20000 50000];
        npredpts = 1000;
        method = {'idw','gpr'}; % 'sphbary', 'pbary', 'gpr', 'sphgpr', 'nn', 'avg'
        datatype = {'brk','kim'};
        pgnum = 32; %m-3m (i.e. m\overbar{3}m) FCC symmetry default for e.g. Ni
        
    case 'full'
        ndatapts = [100 388 500 1000 5000 10000 20000 50000];
        npredpts = 10000;
        method = {'sphgpr','gpr','sphbary','pbary','nn','avg','idw'}; % 'sphbary', 'pbary', 'gpr', 'sphgpr', 'nn', 'avg'
        datatype = {'brk','kim'};
        pgnum = 32; %m-3m (i.e. m\overbar{3}m) FCC symmetry default for e.g. Ni
end

%comment (no spaces, used in filename)
% comment = 'paper-data2';
comment = 'set4';
% comment = 'idw-test-3pt5deg';

% job submission environment
env = 'slurm'; %'slurm', 'local'
T = true;
F = false;
%whether to skip running the jobs and just compile results
dryrunQ = F;
disp(['env = ' env])

if strcmp(env,'slurm') && dryrunQ
    error('did you mean to change dryrunQ to false?')
end

metaQ = T; %whether to load full model or only meta-data at end
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
        end
end

%% functions to generate save filepaths
%diary
files = dir(fullfile('**','data','randOctParity','diary'));
diaryfolder = files(1).folder;
diarynamefn = @(method,ndatapts,gitcommit,puuid) [method int2str(ndatapts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '_' comment '.txt'];
diarypathfn = @(method,ndatapts,gitcommit,puuid) fullfile(diaryfolder,diarynamefn(method,ndatapts,gitcommit,puuid));
%data
files = dir(fullfile('**','data','randOctParity','pcombs'));
savefolder = files(1).folder;
savenamefn = @(method,ndatapts,gitcommit,puuid) [method int2str(ndatapts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '_' comment '.mat'];

%for use with dir
savepathgen = fullfile(savefolder,'*gitID-*puuID*.mat');

savenamematch = [ ...
    ['(' strjoin(method,'|') ')'] ... match (exactly) any of the method options
    ['(' strjoin(cellfun(@num2str,num2cell(ndatapts),'UniformOutput',false),'|') ')'] ... match (exactly) any of the ndatapts options
    '(_gitID-[a-z0-9]*)' ... match any combination of 0 or more lowercase alphabetic or numeric characters (for git hash)
    '(_puuID-[a-z0-9]+)' ... match any combination of 1 or more lowercase alphabetic or numeric characters (for param combo uuid)
    ['(_' comment ')'] ...
    ]; %e.g. '(sphbary|gpr|nn)(50)(_gitID-)[a-z0-9]*)(_puuID-[a-z0-9]+)' matches sphbary50_gitID-abcd123_puuID-abcd123.mat
if metaQ
    savenamematch = [ savenamematch '(_meta.mat)'];
else
    savenamematch = [savenamematch '(.mat)'];
end
    
savepathfn = @(method,ndatapts,gitcommit,puuid) fullfile(savefolder,savenamefn(method,ndatapts,gitcommit,puuid));

if ~dryrunQ
    %% parameter file setup
    %parameters
    pars = var_names(ndatapts,npredpts,method,cores,datatype,pgnum); %**ADD ALL PARAMETERS HERE** (see runtype switch statement)
    %function to execute and output arguments from function
    execfn = @(ndatapts,npredpts,method,datatype,pgnum) ... **NAMES NEED TO MATCH PARS FIELDS** (see above)
        interp5DOF_setup(ndatapts,npredpts,method,datatype,'pgnum',pgnum); %**NAMES NEED TO MATCH PARS FIELDS AND EXECFN ARGUMENTS**
    argoutnames = {'ypred','interpfn','mdl','mdlpars'};
    %i.e. [propOut,interpfn,mdl,mdlpars] = interp5DOF_setup(ndatapts,npredpts,method,'datatype',datatype);
    
    % walltimefn = @() 300; %can set to constant or to depend on parameters, probably fine when using standby queue
    walltimefn = @(ndatapts,npredpts,method,cores) get_walltimefn(ndatapts,npredpts,method,cores);
    
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
            clear mdllist
            %         mdltbl = struct2table(mdlcat,'AsArray',true);
        end
        
        mdlparscat = structvertcat(mdlparslist{:});
        mdlparstbl = struct2table(mdlparscat,'AsArray',true);
        
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
        writetable(mdlparstbl,[fpath(1:end-4) '.csv'])
        disp(mdlparstbl)
        
        %         %nested loop through jobs and tasks to load results
        %         S = load(parpath);
        %         for jid = 1:length(Ntrim)
        %             for tid = 1:jid
        %                 S.
        %             end
        %         end
end

disp('end randOctParityData.m')
disp(' ')
%-----------------------------CODE GRAVEYARD-------------------------------
%{


for i = 1:length(ndataptsList)
    %unpack
    ndatapts = ndataptsList(i);
    for j = 1:length(methods)
        method = methods{j};
        fname = ['randOctParityData_' method int2str(ndatapts) '.mat'];
        S = load(
    end
end

% execfn = @interp5DOF_setup;

%initialize
init1 = cell(length(ndatapts),length(method));
propOutlist = init1;
mdllist = init1;
mdlparslist = init1;
interpfnlist = init1;


% %nested loop through ndataptsList and methodlist
% for i = 1:length(ndatapts)
%     %unpack
%     ndatapts = ndatapts(i);
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
%         [propOutlist{i,j},interpfnlist{i,j},mdllist{i,j},mdlparslist{i,j}]=interp5DOF_setup(ndatapts,npredpts,method,uuid);
%
%         %% output
%         %package
%         propOut = propOutlist{i,j};
%         interpfn = interpfnlist{i,j};
%         mdl = mdllist{i,j};
%         mdlpars = mdlparslist{i,j};
%         %save
%         save(savepathfn(method,ndatapts,mdl.gitcommit,uuid),...
%             'propOut','interpfn','mdl','mdlpars','ndatapts','method','uuID',...
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

pars = var_names(ndatapts,npredpts,method.uuid);

%         S(nfiles) = struct();
%         S(1).out = struct();


memfn = @(method,ndatapts) ...
    (strcmp(method,'gpr') && ndatapts > 10000) *1024*8 + ...
    (~strcmp(method,'gpr') && ndatapts <= 10000) * 1024*4;

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

%}
