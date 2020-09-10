function sbatch_submit(jbinfo,sbopts)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description:
%
% Inputs:
%		load_filepath		===	parameter file
%
%		failed_filepath	===	failed parameter filepath
%
%		sbopts				===	sbatch options
%
% Outputs:
%
% Dependencies:
%		get_cmd.m
%
% Notes:
%		Depends on having a parameter file for the jobs and a single result
%		file for each job. The existence of result files is used to determine
%		which jobs failed.
%
%--------------------------------------------------------------------------
disp(' ')
%% Package info & options into struct
S.jbinfo = jbinfo;
S.sbopts = sbopts;

%% Assign variables from struct for brevity
opt2 = 'eval'; %'eval' or 'manual'
switch opt2
	case 'eval'
		%see footnote [1]
		snames = fields(S); %struct names
		
		%load field names
		for i = 1:length(snames)
			vnames{i} = fields(S.(snames{i})); %variable names
		end
		
		%check if structures share field names
		if ~isempty(intersect(vnames{:}))
			warning('overlapping fieldnames betweeen structs..')
		end
		% unpack fields from struct
		for i = 1:length(snames)
			sname = snames{i};
			for j = 1:length(vnames{i})
				vname = vnames{i}{j};
				temp = S.(sname).(vname); %temporary value of vName
				evalc([vname '= temp']); %assign temp value to a short name
			end
		end
		
	case 'manual'
		%see footnote [2]
		N_trim = jbinfo.N_trim;
		walltimes = jbinfo.walltimes;
		failedjobQ = jbinfo.failedjobQ;
		failed_filepath = jbinfo.failed_filepath;
		env = jbinfo.env;
		
		cores=sbopts.cores;
		mem=sbopts.mem; %MB
		qosopt= sbopts.qosopt; %e.g. ' --qos=test ';
end

%% ------------run/submit jobs------------
Nidx_length = length(N_trim);

setenv('for_type',for_type)
disp(['for_type = ' getenv('for_type')])

setenv('cores',sprintf(int2str(cores)))
disp(['cores = ' getenv('cores')])

if ~failedjobQ
	setenv('load_filepath',['"',load_filepath,'"'])
else
	walltimes = int2str(timefactor*walltimes);
	setenv('load_filepath',['"',failed_filepath,'"'])
end

disp(['load_filepath = ' getenv('load_filepath')])
disp(' ')

switch env
	case 'slurm'
		%-------submit sbatch scripts--------
		for jid = 1:Nidx_length
			%get sbatch command
			cmd = get_cmd(N_trim,jid,cores,walltimes,mem,qosopt,dirpath,script_fpath);
			disp(cmd)
			%send command to shell
			[status,cmdout] = system(sprintf(cmd),'-echo');
			disp(' ');
		end
		
	case 'local'
		if ~failedjobQ
			filepath = load_filepath(2:end-1);
		else
			filepath = failed_filepath(2:end-1);
		end
		for jid = 1:Nidx_length
			time2 = [];
			disp(['jid = ',int2str(jid)]);
			for tid = 1:N_trim(jid) %[119 120 121 122 140 146 154]
				disp(['tid = ',int2str(tid)])
				tic
				sse_inference(filepath,jid,tid)
				time2 = [time2 toc]
			end
		end
		
end
end

%-----------------------------HELPER FUNCTIONS-----------------------------

%------------------------------FOOTNOTES-----------------------------------
%{

[1] normally I avoid eval statements AT ALL COSTS. In this case, it may
actually make sense to use because it allows me to add more variables
without having cumbersome naming assignments all throughout (e.g.
T(k).ncycles) or having to manually update the new variable assignment in
more than one place. When using eval, only the statement above that assigns
to vnames and the corresponding lists would need to be changed. That is, by
adding the new variable name to vnames, and the corresponding list to the
'T' struct, everything else should go through.

[2] manual alternative to the eval statement above. Can copy vnames as
comment somewhere above for reference or display it on the command line to
avoid excessive scrolling back and forth when double checking all variables
are accounted for
%}


%-------------------------------CODE GRAVEYARD-----------------------------
%{
				%find and format corresponding tid's
				%tid_list = cellfun(@(num) int2str(num),...
				%	num2cell(sort(failed.tid(failed.jid==jid))),...
				%	'UniformOutput',false);	%e.g. {'28', '36', '40', ...}
				%tidstr_list = strjoin(tid_list{jid},','); %e.g. '28,36,40,...'
				
				
				setenv('jid',sprintf(int2str(jid)))
				walltimes = int2str(timefactor*walltimes(jid));
				setenv(sprintf(['walltime',int2str(jid)]),sprintf(walltime));
				
				%----------sbatch options-----------
				arrayopt = ['sbatch -a ',tidstr_list];
				coreopt = [' -c ',int2str(cores)];
				qosopt = [ qos ];
				timeopt = [' --time=',walltime,':00'];
				imageopt = [' -C rhel7 '];
				jobnameopt = [' --job-name=batch',int2str(jid)];
				memopt = [' --mem-per-cpu=',int2str(mem)];
				diropt = [' -D ~/compute/1DOF_functional'];
				scriptopt = [ script_fpath ];

				bckgdopt = ''; %[' &'];

%}