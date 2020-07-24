
seed_type = 'previous'; %'new', 'previous'
failedjobQ = true;

%variables to save and load for 'find' and 'from_file', respectively
failed_vars = {'N_trim','Nidx_length','parcombs_sets','ninf',...
	'job_walltimes','select_idx','ng_idx','nObs_idx','L_idx',...
	'seed_idx','rand_id'};