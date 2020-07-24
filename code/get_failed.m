%find failed jobs

%% determine failed jobs (if applicable)
if strcmp(seed_type,'previous')
	if failedjobQ
		%variables to save and load for 'find' and 'from_file', respectively
		failed_vars = {'N_trim','Nidx_length','parcombs_sets','ninf',...
			'job_walltimes','select_idx','ng_idx','nObs_idx','L_idx',...
			'seed_idx','rand_id'};
		
		failedjob_type = 'find'; %'find', 'from_file'
		
		switch failedjob_type
			case 'find'
				ngnss_timefactor = 20/5000;
				%get full list
				varid_list = [select_idx ng_idx nObs_idx L_idx seed_idx];
				table_all = parcombs3_sorted(:,varid_list);
				
				%find successful result files
				result_folder = fullfile('result');
				%result_filenamegen = ['result_*_ng*_nObs*_L*_ninf',int2str(ninf),'_jid*_tid*_id*_seed*_rid',int2str(rand_id),'.mat'];
				result_filenamegen = ['result_*_ng*_nObs*_L*_ninf',int2str(ninf),'_jid*_tid*_id*_rid',int2str(rand_id),'.mat'];
				result_filepath = fullfile(result_folder,result_filenamegen);
				
				result_files = dir(result_filepath);
				nfiles = length(result_files);
				
				result_folders = {result_files.folder};
				result_filenames = {result_files.name};
				
				f = @(i) fullfile(result_folders{i},result_filenames{i});
				result_filepaths = arrayfun(f,1:nfiles,'UniformOutput',false); %compile filepaths
				
				table_success = cell(nfiles,5); %initialize
				loadparfor = true;
				if ~loadparfor
					msg = 'loading successful parameter combinations';
					wb = waitbar(0,msg);
					for i = 1:nfiles
						temp=load(result_filepaths{i},'table_full');
						temp = temp.('table_full');
						table_success(i,:) = ...
							{temp.SelectionMethod{1}, temp.ng, temp.nObs, temp.L, temp.seed};
						waitbar(i/nfiles,wb,msg);
					end
					close(wb);
				else
					parfor i = 1:nfiles
						temp=load(result_filepaths{i},'table_full');
						temp = temp.('table_full');
						table_success(i,:) = ...
							{temp.SelectionMethod{1}, temp.ng, temp.nObs, temp.L, temp.seed};
					end
				end
				
				%determine failed parameter combinations
				[~,failedid_list] = setdiff(string(table_all),string(table_success),'rows');
				table_failed = parcombs3_sorted(failedid_list,:);
				nfailed = size(table_failed,1);
				disp(['nfailed = ',int2str(nfailed)])
				
				% Prompt to continue
				m = input('Does this seem correct? y/n:','s');
				if ~strcmp(m,'y') && ~strcmp(m,'Y')
					return
				end
				
				%---create new parameter sets---
				ngnss_tot = ([table_failed{:,2}].*[table_failed{:,3}])';
				[~,sortid_list] = sort(ngnss_tot);
				tablefailed_sorted = table_failed(sortid_list,:);
				walltimes = ngnss_timefactor*ngnss_tot*walltime_safetyfactor;
				
				%--------sorting parameters------
				max_walltime = 1; %minutes
				max_tasks = 500;
				max_taskstot = 5000;
				walltime_cushion = 60; %minutes
				
				%------------sorting-------------
				[N_trim, Nidx_length, parcombs_sets, job_walltimes] = ...
					sortbins(walltimes, max_walltime, max_tasks, max_taskstot,...
					tablefailed_sorted, walltime_cushion);
				
				%---save failed parameter combos--
				%--find existing files--
				failed_folder = fullfile('parameters');
				failed_filenamegen = ['failed_rid',int2str(rand_id),'*.mat'];
				failed_filepathgen = fullfile(failed_folder,failed_filenamegen);
				
				nfailed_files = length(dir(failed_filepathgen));
				%^^find existing files^^
				
				%{
							--------filenaming--------
						If 'failed_rid<rand_id>.mat' exists, then save a file named
							'failed_rid<rand_id>_1.mat'.
						If 'failed_rid<rand_id>.mat' and 'failed_rid<rand_id>_1.mat'
							both exist, then save 'failed_rid<rand_id>_2.mat'.
						etc.
				%}
				if nfailed_files == 0
					failed_filename = ['failed_rid',int2str(rand_id),'.mat'];
					failed_filepath = fullfile(failed_folder,failed_filename);
				else
					failed_filename = ['failed_rid',int2str(rand_id),'_',int2str(nfailed_files),'.mat'];
					failed_filepath = fullfile(failed_folder,failed_filename);
				end
				
				save(failed_filepath,failed_vars{:});
				%^^^save failed parameter combos^^
				
			case 'from_file'
				%----general filepath----
				failed_folder = fullfile('parameters');
				failed_filenamegen = ['failed_rid',int2str(rand_id),'*.mat'];
				failed_filepathgen = fullfile(failed_folder,failed_filenamegen);
				
				%-----take last file-----
				failed_files = dir(failed_filepathgen);
				failed_names = {failed_files.name};
				failednames_sorted = sort(failed_names);
				failed_filename = failednames_sorted{end};
				
				%--------filepath--------
				failed_filepath = fullfile(failed_folder,failed_filename);
				
				disp('loading failed job information')
				load(failed_filepath,failed_vars{:});
				
				nfailed = size(vertcat(parcombs_sets{:}),1);
				disp(['nfailed = ',int2str(nfailed)])
				
				% Prompt to continue
				m = input('Does this seem correct? y/n:','s');
				if ~strcmp(m,'y') && ~strcmp(m,'Y')
					return
				end
		end
	end
end