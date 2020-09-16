clear; close all

addpathdir({'var_names.m','writeparfile.m'})
runtype = 'full'; %'test','full'
switch runtype
    case 'test'
        ndatapts = 388;
        npredpts = 500;
        method = {'sphbary','gpr'};
        inputtype = {'5dof','octonion'}; %'5dof','octonion'
    case 'full'
        ndatapts = [10000 20000];
        npredpts = 10000;
        method = {'gpr'}; %'sphbary','pbary','gpr','nn'
        inputtype = {'5dof'};
end

%% functions to generate save filepaths
%diary
diaryfolder = fullfile('data','randOctParity','diary');
diarynamefn = @(method,ndatapts,gitcommit,puuid) [method int2str(ndatapts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '.txt'];
diarypathfn = @(method,ndatapts,gitcommit,puuid) fullfile(diaryfolder,diarynamefn(method,ndatapts,gitcommit,puuid));
%data
savefolder = fullfile('data','randOctParity','pcombs');
savenamefn = @(method,ndatapts,gitcommit,puuid) [method int2str(ndatapts) '_gitID-' gitcommit(1:7) '_puuID-' puuid '.mat'];

savepathgen = fullfile(savefolder,'*gitID-*puuID*.mat'); %for use with dir
savenamematch = [... %for use with regexp
    ['(' strjoin(method,'|') ')'] ... match (exactly) any of the method options
    ['(' strjoin(num2cell(num2str(ndatapts),2),'|') ')'] ... match (exactly) any of the ndatapts options
    '(_gitID-[a-z0-9]*)' ... match any combination of 0 or more lowercase alphabetic or numeric characters (for git hash)
    '(_puuID-[a-z0-9]+)' ... match any combination of 1 or more lowercase alphabetic or numeric characters (for param combo uuid)
    '[.mat]' ... %match (exactly) the file-extension
    ]; %e.g. '(sphbary|gpr)(50)(_gitID-)[a-z0-9]*)(_puuID-[a-z0-9]+).mat'

savepathfn = @(method,ndatapts,gitcommit,puuid) fullfile(savefolder,savenamefn(method,ndatapts,gitcommit,puuid));

%% parameter file setup
%parameters
pars = var_names(ndatapts,npredpts,method,inputtype); %add all parameters here (see runtype switch statement)
%function to execute and output arguments from function
execfn = @(ndatapts,npredpts,method,inputtype) ...
    interp5DOF_setup(ndatapts,npredpts,method,get_uuid(),inputtype); %names need to match pars fields
argoutnames = {'propOut','interpfn','mdl','mdlpars'};
%i.e. [propOut,interpfn,mdl,mdlpars] = interp5DOF_setup(ndatapts,npredpts,method);

walltimefn = @() 600; %can set to constant or to depend on parameters

%% parameter file
[parpath, parcombsets, Ntrim, jobwalltimes] = ...
    writeparfile(pars,execfn,argoutnames,walltimefn,'diarypathfn',diarypathfn,'savepathfn',savepathfn);

%% job submission
env = 'slurm'; %'slurm', 'local'
switch env
    case 'slurm'
        %setup
        cores = 6;
        mem = 1024*8*cores; %total memory of job, MB
        qosopt = 'standby'; %'', 'test', 'standby'
        scriptfpath = fullfile('MATslurm','code','submit.sh');
        %submission
        submit_sbatch(parpath,cores,mem,qosopt,scriptfpath); %submit.sh --> exec_combs.m --> execfn (defined above)
        
    case 'local'
        %nested loop through jobs and tasks to generate results
        for jid = 1:length(Ntrim)
            N = Ntrim(jid);
            for tid = 1:N
                outtmp = exec_combs(parpath, jid, tid);
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
        assert(nfiles ~= 0,'no files matched')
        %initialize
        init1 = cell(nfiles,1);
        [outlist,propOutlist,mdllist,mdlparslist,interpfnlist] = deal(init1);
        
        %extract results from files
        for i = 1:nfiles
            folder = folders{i};
            name = names{i};
            fpath = fullfile(folder,name);
            clear S out
            S = load(fpath,'out');
            out = S.out;
            [propOutlist{i},mdllist{i},mdlparslist{i},interpfnlist{i}] = ...
                deal(out.propOut, out.mdl, out.mdlpars, out.interpfn);
            outlist{i} = out;
        end
        
        %concatenate models and parameters
        mdlcat = structvertcat(mdllist{:});
        mdlparscat = structvertcat(mdlparslist{:});
        mdlparstbl = struct2table(mdlparscat,'AsArray',true);
        
        gitcommit = get_gitcommit();
        %save models and parameters
        fpath = fullfile(savefolder,['gitID-' get_gitcommit '_uuID-' get_uuid() '.mat']);
        save(fpath,'propOutlist','interpfnlist','mdllist',...
            'mdlparslist','mdlcat','mdlparscat','mdlparstbl','outlist',...
            '-v7.3','-nocompression')
        writetable(mdlparstbl,[fpath(1:end-4) '.xlsx'])
        
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
%}
