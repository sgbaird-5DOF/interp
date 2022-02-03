function [parpath, parcombsets, Ntrim, jobwalltimes] = ...
    writeparfile(pars,execfn,argoutnames,walltimefn,uuid,excludefpath,NV)
arguments
    pars struct
    execfn function_handle %#ok<INUSA> %gets saved in parameter file
    argoutnames cell %#ok<INUSA> %gets saved in parameter file
    walltimefn function_handle = @() 60
    uuid char = get_uuid()
    excludefpath char = ''
    NV.maxwalltime(1,1) double = 59
    NV.maxtasks(1,1) double = 500
    NV.maxtaskstot(1,1) double = 5000
    NV.walltimecushion(1,1) double = 20
    NV.nreps(1,1) double = 1
    NV.diarypathfn function_handle = @(puuid) fullfile('data','diary',[puuid '.txt']); %#ok<FVUBD>
    NV.savepathfn function_handle = @(puuid) fullfile('data','pcombs',[puuid '.txt']); %#ok<FVUBD>
end
% WRITEPARFILE generate a parameters file for an sbatch submission
%--------------------------------------------------------------------------
% Inputs:
%  pars - struct containing parameter lists, where the lists can be cells
%  or numeric, and a single row corresponds to a single parameter out of
%  one of the parameter lists.
%
%  NV - name-value pairs
%   nreps - number of repeat runs to perform
%
%   walltimefn - function handle to compute a walltime based on a parameter
%   combination. Make sure that the names of the parameters in the
%   functional handle are consistent with the names of the parameters in
%   pars. The function handle can include a call to a user-defined .m file
%   if needed.
%
% Outputs:
%  parpath - relative filepath to parameter file
%
% Usage:
%  parpath = writeparfile(pars,execfn,argoutnames,walltimefn)
%
% Dependencies:
%  allcomb.m
%
% Notes:
%  Next is submit_sbatches.m
%
% Author(s): Sterling Baird
%
% Date: 2020-09-09
%--------------------------------------------------------------------------
diarypathfn = NV.diarypathfn; %#ok<NASGU>
savepathfn = NV.savepathfn; %#ok<NASGU>

%% parameter combinations setup

%convert struct to cell
parcell = namedargs2cell(pars);
parnames = parcell(1:2:end);
parvals = parcell(2:2:end);

%convert numerics and chars to cells
fns = {@isnumeric,@ischar,@islogical};
changefns = {@(x) num2cell(x,1),@(x) num2cell(x,2),@(x) num2cell(x,1)};
for i = 1:length(fns)
    %unpack
    changefn = changefns{i};
    fn = fns{i};
    
    numIDs = cellfun(fn,parvals);
    nums = cellfun(@(parval) changefn(parval),parvals,'UniformOutput',false);
    tmp = nums(numIDs);
    [parvals{numIDs}] = tmp{:};
end

%check for non-cell parameter lists
if any(~cellfun(@iscell,parvals))
    error(['Parameter lists should be in cell format. Inputs can only be cell or numeric vectors. Noncell names: ' ...
        strcat(parnames{~cellfun(@iscell,parvals)})])
end

%% parameter combinations
parcombstmp = allcomb(parvals{:});
parcombstmp = repmat(parcombstmp,NV.nreps,1); %repeat inferences
ncombs = size(parcombstmp,1);

parcombs = cell2struct(parcombstmp,parnames,2);

%% excluded parameter combinations
if ~isempty(excludefpath) && (exist(excludefpath,'file') == 2)
    S = load(excludefpath,'parcombs');
    excludeparcombs = S.('parcombs');
    T1 = struct2table(parcombs);
    T2 = struct2table(excludeparcombs);
%     T1 = cell2table(parcombs);
%     T2 = cell2table(excludeparcombs);
%     table1 = cell2table(parcombs,'VariableNames',{'sel','ng','nObs','L'});
%     table2 = cell2table(excludeparcombs,'VariableNames',{'sel','ng','nObs','L'});
    uniquesets = setdiff(T1,T2);
    parcombs = table2cell(uniquesets);
elseif ~isempty(excludefpath) && (exist(excludefpath,'file') ~= 2)
    disp(['exclude file ' excludefpath ' not found'])
end

%% walltimes
%compute
walltimes = get_walltimes(parcombs,walltimefn);
%package
parcombs(1).walltime = [];
tmp = num2cell(walltimes);
[parcombs.walltime] = tmp{:};

% parcombs3 = parcombs2;
% parcombs3(:,walltimes_idx) = num2cell(walltimes);

%% random seed for each parameter combinations
masterseed = rng('shuffle');
%seeds
seeds = randi(100000,ncombs,1);
%package
parcombs(1).seed = [];
tmp = num2cell(seeds);
[parcombs.seed] = tmp{:};

%% git commit hash for each parameter combination
[gitcommittmp,~,warnedQ] = get_gitcommit();
if warnedQ
    m = input(['Warning issued from get_gitcommit. Commit/push as needed. Continue (y) or abort (n)? '],'s');
    if ~strcmp(m,'y') && ~strcmp(m,'Y')
        error('User chose to abort')
    end
    gitcommittmp = get_gitcommit();
end
gitcommit = repelem(gitcommittmp,ncombs,1);
%package
parcombs(1).gitcommit = '';
tmp = num2cell(gitcommit,2);
[parcombs.gitcommit] = tmp{:};

%% random unique ID for each parameter combination
puuid = cell(ncombs,1);
for i = 1:ncombs
    puuid{i} = get_uuid();
end
%package
parcombs(1).puuid = '';
[parcombs.puuid] = puuid{:};

%% sort parameter combinations into task bins
%get appropriate name-value pairs in struct format (pare down from NV)
sortbinsNV = struct('maxwalltime',NV.maxwalltime,'maxtasks',NV.maxtasks,...
    'maxtaskstot',NV.maxtaskstot,'walltimecushion',NV.walltimecushion);
%convert to cell name value pairs
NVpairs = namedargs2cell(sortbinsNV);
%sort into bins
[Ntrim, parcombsets, jobwalltimes] = sortbins(parcombs, walltimes, NVpairs{:});

%% saving
%data filepath
datafolder = fullfile('templates');
datafilename = @(ng) ['.mat'];
datafilepath = @(ng) fullfile(datafolder,datafilename(ng));

%parameter
disp('saving parameter file')
parfolder = fullfile('parameters');
if exist(parfolder,'dir') ~= 7
    mkdir(parfolder)
end
parname = ['uuID-',uuid,'.mat'];
parpath = fullfile(parfolder,parname);
save(parpath)
% parpath = ['"',savefilepath,'"'];

end

%-------------------------helper functions---------------------------------
%%
%^^^^^^^^^^^^^^^^^^^^^^^^helper functions^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%%
%{
% --------------Code Graveyard----------------
%  tic;
% 	load(data_filepath(ng))
% 	time.load = toc;
% 	disp([int2str(time.load),' s to load data'])
%


%---------export filepath to shell--------
% PATH = getenv('PATH');
% setenv('PATH', [PATH ':/lustre/scratch/usr/sbaird9/1DOF_functional']);
% unix('export ninf=10; matlab sse_inference')
%
% setenv('save_filepath',save_filepath)


%binrange{1} = edges(1):edges(2);

% 	[edge_min,edge_max] = bounds(ngnss_tot);
% 	comp_binsize = (edge_max-edge_min)/(N);
% 	edges = edge_min:comp_binsize:edge_max;

%parcombs_sets{i} = sortrows(parcombs_sets{i},ngnsstot_idx);

% 	%----separate for all but last bin------
%
% 	for i = 1:length(edges)-2
% 		temp_idx = edges(i) <= temp & temp < edges(i+1);
% 		parcombs_sets{i} = parcombs2_sorted(temp_idx,:);
% 		parcombs_sets{i} = sortrows(parcombs_sets{i},ngnsstot_idx);
% 		if ~isempty(parcombs_sets{i})
% 			job_walltimes(i) = 9/5000*cell2mat(parcombs_sets{i}(end,ngnsstot_idx))*walltime_safetyfactor+10;
% 		end
% 		njobs(i) = size(parcombs_sets{i},1);
% 	end
%
% 	%-------separate for last bin-----------
% 	i = length(edges)-1;
% 	temp_idx = edges(i) <= temp & temp < edges(i+1);
% 	parcombs_sets{i} = parcombs2_sorted(temp_idx,:);
% 	if ~isempty(parcombs_sets{i})
% 		job_walltimes(i) = 9/5000*cell2mat(parcombs_sets{i}(end,ngnsstot_idx))*walltime_safetyfactor+20;
% 	end
% 	njobs(i) = size(parcombs_sets{i},1);


% 	for i = 1:Nidx_length - 1
% 		walltimes2(i) = cell2mat(parcombs_sets{i}(:,walltimes_idx));
% 		nbins2(i) = histcounts(walltimes2(i));
% 	end


% 	while Ntot > max_taskcount
% 		for i = 1:Nidx_length-1
% 			while max(Nsub{i}) > max_binsize || walltimes2(i) > max_walltime
% 				nbins2(i) = nbins2(i) + 1;
% 				[Nsub{i},edges] = histcounts(walltimes_temp{i},nbins2(i));
%
% 				Nidx2 = find(Nsub{i}~=0);
% 				Nidx_length2 = length(Nidx2);
%
% 				for j = 1:Nidx_length2-1
% 					%find elements corresponding to each sub-bin
% 					edge1 = edges(Nidx2(j));
% 					edge2 = edges(Nidx2(j)+1);
% 					temp_idx = edge1 <= temp & temp < edge2;
%
% 					%check to see if max_walltime is violated
% 					walltimes_temp2(j) = sum(walltimes_temp{i}(temp_idx));
% 				end
%
% 				walltimes2(i) = max(walltimes_temp2);
%
% 				Ntot_temp(i) = sum(Nsub{i});
% 				if sum(Ntot_temp) <= max_taskcount ...
% 						&& max(Nsub{i}) <= max_binsize && walltimes2(i) <= max_walltime
% 					break
% 				end
% 			end
% 		end
% 		Ntot = sum(Ntot_temp);
% 	end


% 		max_binsize = 1000; %maximum number of elements in a sub-bin (i.e. # loop elements in sse_inference.m)


% %--sort into bins with equal width--
% 		N = histcounts(timetot_list);
% 		iter_factor = 0.5;
% 		N_temp = ceil(max(N)/iter_factor);
% 		nbins = length(N);
% 		nbins_factor = 2;
%
% 		disp('sorting into bins')
% 		max_taskcount = 5000; %maximum number of total tasks to submit (limit on BYU supercomputer is 5000 as of 2019-07-08)
% 		Ntot = sum(N); %initialize
% 		max_bins = 500;
% 		ktot = max_bins+1;
% 		while max(ktot) > max_bins
% 			N_temp = ceil(N_temp*iter_factor);
% 			disp(['N_temp = ',int2str(N_temp)])
% 			while(max(N) > N_temp)
% 				nbins = nbins*nbins_factor;
% 				N = histcounts(walltimes,nbins);
% 			end
%
% 			%while(sum(N) > 5000)
% 			[N, edges] = histcounts(walltimes,nbins);
%
% 			Nidx = find(N~=0);
% 			Nidx_length = length(Nidx);
%
% 			N_trim = N(Nidx);
% 			%end
%
%
% 			%----separate for all but last bin------
% 			for i = 1:Nidx_length-1
% 				edge1 = edges(Nidx(i));
% 				edge2 = edges(Nidx(i)+1);
% 				temp_idx = edge1 <= temp & temp < edge2;
% 				parcombs_sets{i} = parcombs3_sorted(temp_idx,:);
% 				if ~isempty(parcombs_sets{i})
% 					job_walltimes(i) = 8/5000*cell2mat(parcombs_sets{i}(end,ngnsstot_idx))*walltime_safetyfactor;
% 				end
% 			end
%
% 			%-------separate for last bin-----------
% 			i = Nidx_length;
% 			edge1 = edges(Nidx(i));
% 			edge2 = edges(Nidx(i)+1);
% 			temp_idx = edge1 <= temp & temp <= edge2;
% 			parcombs_sets{i} = parcombs3_sorted(temp_idx,:);
% 			if ~isempty(parcombs_sets{i})
% 				job_walltimes(i) = 8/5000*cell2mat(parcombs_sets{i}(end,ngnsstot_idx))*walltime_safetyfactor;
% 			end
%
% 			%-----separate into sub-bins-------
%
% 			%-----------job loop------------
% 			for i = 1:Nidx_length
% 				walltimes_temp = cell2mat(parcombs_sets{i}(:,walltimes_idx));
% 				%initialize counters
% 				j = 0; %nth parameter combo within job
% 				k = 0; %number of tasks (each task is associated with a for loop containing sublist_length(i,k) elements)
% 				timecounter_full = []; %running list of times for each parameter grouping
% 				%--sort parameter groups--
% 				while j < length(walltimes_temp)
% 					j0 = j; %initialize starting j value for each group
%
% 					time_counter = 0;
% 					if walltimes_temp(j+1) > max_walltime || length(walltimes_temp) == 1
% 						j = j+1;
% 						time_counter = walltimes_temp(j);
% 					end
% 					%--create parameter group--
% 					while time_counter <= max_walltime
% 						%--exit loop if all par combos in job done--
% 						if j == length(walltimes_temp) || length(walltimes_temp) == 1
% 							break
% 						end
% 						j = j+1;
%
% 						time_counter = time_counter + walltimes_temp(j);
% 						%--reverse last iteration--
% 						if time_counter > max_walltime
% 							time_counter = time_counter - walltimes_temp(j);
% 							j = j-1;
% 							break
% 						end
% 					end
% 					k = k+1;
% 					timecounter_full = [timecounter_full time_counter];
% 					jf = j; %final parameter combo within group
% 					subidx_list = j0+1:jf;
% 					sublist_length{i,k} = jf-j0; %number of parameter combos within parameter group
% 					parcombs_sets{i}(subidx_list,subbin_idx) = {k}; %take ###(#==k) in sse_inference.m
% 				end
% 				ktot(i) = k; %total number of parameter groups
% 				walltimes2(i) = max(timecounter_full)+20; %final walltime for job
% 			end
% 			disp(['max(ktot) = ',int2str(max(ktot))])
% 		end

% 					if i == length(id_list)-1
% 						break
% 					end

% 			i = 0;
% 			j = 0;
% 			k = 0;
% 			tid_list = [];
% 			while i < length(id_list)
% 				tidlist_temp = [];
% 				k = k+1;
% 				while j < max_tasks
% 					i = i+1;
% 					j = j+1;
% 					tidlist_temp = [tidlist_temp i];
% 					tid_list{k} = idlist_temp2; %#ok<*AGROW>
% 				end
% 			end

% subbin_idx = 8;


% 				if ~isempty(parcombs_sets{i})
% 					job_walltimes(i) = 8/5000*cell2mat(parcombs_sets{i}(end,ngnsstot_idx))*walltime_safetyfactor;
% 				end

seedlisttype = 'master' %'master' (faster) or 'individual' (longer)
%note: I haven't noticed any qualitative difference between the two.
switch seedlisttype
    case 'master'
        masterseed = rng('shuffle');
        %seeds
        seeds = randi(10000000,ncombs,1);
        %package
        parstruct.seeds = seeds;
%         parcombs3(:,seed_idx) = num2cell(seedlist);
        
    case 'individual'
        f = waitbar(0);
        seeds{size(parcombs3,1)} = rng('shuffle');
        for i = 1:size(parcombs3,1)
            waitbar(i/size(parcombs3,1),f,[int2str(i),'/',int2str(size(parcombs3,1))])
            seeds{i} = rng('shuffle');
        end
        close(f)
        parcombs3(:,seed_idx) = seeds;
end

%% sort the full cell array by walltime
%find sorted indices
[sortedwalltimes,idx] = sortrows(parstruct.walltimes);
%rearrange parameter structure
sortedparstruct = parstruct(idx);

% parcombs3_sorted = sortrows(parcombs3,walltimes_idx);
% temp = cell2mat(parcombs3_sorted(:,walltimes_idx)); %list of sorted
% walltimes?



%parcombs = allcomb(SelectionMethod_list,num2cell(ng_list),num2cell(nObs_list),num2cell(L_list));
%}


%{
	USER: make sure the idx's are correct (i.e. they correspond to the
	correct column in the parameter combination array - par_combs, some
	assigned before here, and some assigned later. These variables are
	loaded from the saved parameter file in run.m 2020-05-05 SB)
%}
% select_idx = 1;
% ng_idx = 2;
% nObs_idx = 3;
% L_idx = 4;
% ngnsstot_idx = 5;
% walltimes_idx = 6;
% seed_idx = 7;

%---------calculate est. time----------
% disp('estimate walltimes')
% ngnsstot = cell2mat(parcombs(:,ng_idx)).*cell2mat(parcombs(:,nObs_idx));
% 
% parcombs2 = parcombs;
% parcombs2(:,ngnsstot_idx) = num2cell(ngnsstot);
% 
% walltime_safetyfactor = 1.2;
% ngnss_timefactor = 20/5000;
% walltimes = ngnss_timefactor*ngnsstot*walltime_safetyfactor;

% walltimefntxt = func2str(walltimefn);
% startID = 3;
% endID = strfind(walltimefntxt,')')-1;
% walltimefnvars = strsplit(walltimefntxt(startID:endID),',');


% maxwalltime = 1; %minutes (1 hr jobs on BYU supercomputer run quickly, cushion can be added as well)
% maxtasks = 500;
% maxtaskstot = 5000;
% walltimecushion = 20; %minutes

%{
In the following example, ng_idx = 3, and nss_idx = 4 (i.e. positions that
they appear in in the allcomb() function)
e.g. par_combs = allcomb(selectionMethod_list,ninf_list,ng_list,nObs_list);
%}

% ^^^^^^^^^^^^^^Code Graveyard^^^^^^^^^^^^^^^^
