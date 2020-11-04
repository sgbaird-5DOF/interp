% SBATCH_SETUP  setup for SLURM sbatch submissions (deprecated)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
% 
% Inputs:
%
% Outputs:
%
% Dependencies:
%
%--------------------------------------------------------------------------

%load failed jobs

rng('shuffle')
rand_id = randi(100000);
%%
%--------options and setup------------
%{
runtype_list can take on individual cell values of: 'sse_dataprep',
'sse_inference', 'sse_compile', 'sse_plot' 'sse_compile', 'sse_plot',
'sse_timingfn', 'jspace_plot'
as in: runtype_list = {'sse_inference','sse_compile','sse_plot'};
%}
runtype = 'infonly'; 


dataset_type = 'full'; %'test' or 'full'
environment = 'slurm'; %'slurm', 'local'
failedjobQ = false; %true if repeating failed jobs, false if otherwise
seed_type = 'previous'; %'new' or 'previous', set to 'previous' if failedjobQ == true

if failedjobQ
	Qstr = 'FAILED dataset from a '; %#ok<*UNRCH>
else
	Qstr = '';
end

disp('-----------Options-----------')
disp(['Do a(n) ',upper(char(join(runtype))),' using a ',Qstr,upper(dataset_type),' dataset on ',upper(environment)])


%% find or load failed jobs
disp(['seed_type: ',upper(seed_type)])
if strcmp(seed_type,'previous')
	rand_id = 51848; %specify rand_id if rerunning failed jobs
end

%% files to exclude
exclude_folder = ''; %'parameters';
exclude_filename = ''; %'parameters_ninf10_rid2933.mat';
exclude_filepath = fullfile(exclude_folder,exclude_filename);

%% retrieve job info
switch seed_type
	case 'new'
		[load_filepath, Nidx_length, N_trim, job_walltimes] = ...
			sse_parameters(SelectionMethod_list,ninf,ng_list,nObs_list,L_list,rand_id,exclude_filepath);
	case 'previous'
		load_folder = ''; %'parameters'
		load_filename = ['parameters.mat'];
		load_filepath = [fullfile(load_folder,load_filename)];
		
		varname_list = {'Nidx_length','N_trim','job_walltimes',...
			'select_idx','ng_idx','nObs_idx','L_idx','seed_idx',...
			'parcombs3_sorted','walltime_safetyfactor'}; %,'ngnss_timefactor'
		load(load_filepath,varname_list{:}) %'Nidx_length','N_trim','job_walltimes')
end

%% find failed jobs



%% submit sbatch
qosopt = ''; %'', '--qos=test'
dirpath = './'; % directory to change to for submission/execution
%package output into structs

jbinfo = var_names(N_trim,walltimes,failedjobQ,env,for_type,load_filepath,timefactor);
sbopts = var_names(cores,mem,qosopt,dirpath,script_fpath);

submit_sbatch(jbinfo,sbopts)


%-----------------------------HELPER FUNCTIONS-----------------------------
function out = var_names(varargin)
%--------------------------------------------------------------------------
% https://www.mathworks.com/matlabcentral/answers/79281#answer_89015 
%--------------------------------------------------------------------------
for n = 1:nargin
	out.(inputname(n)) = varargin{n};
end

end %var_names



%---------------------------CODE GRAVEYARD---------------------------------
%{
switch runtypelist_type
	case 'all'
		runtype_list = {'sse_inference','sse_compile','sse_plot'};
	case 'infonly'
		runtype_list = {'sse_inference'};
	case 'compplot'
		runtype_list = {'sse_compile','sse_plot'};
	case 'plotonly'
		runtype_list = {'sse_plot'};
	case 'componly'
		runtype_list = {'sse_compile'};
end


% jb_vnames = {'N_trim','walltimes','failedjobQ','env','for_type','load_filepath','timefactor'};
% sb_vnames = {'cores','mem','qosopt','dirpath','script_fpath'};
% % jbinfo struct
% for vname = jb_vnames
% 	evalc([jbinfo.(vname) '= vname']);
% end
% 
% % sbopts struct
% for vname = sb_vnames
% 	evalc([sbopts.(vname) '= vname']);
% end

%}