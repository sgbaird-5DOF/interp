function cmd = get_cmd(N_trim,jid,cores,nworkers,walltimes,mem,qosopt,scriptfpath,dirpath)
arguments
    N_trim = 1
    jid = 1
    cores = 24
    nworkers = cores
    walltimes = 60
    mem double = 1024*4*cores
    qosopt char = 'standby'
    scriptfpath char = 'submit.sh'
    dirpath char = '../../' %just above MATslurm
end
% GET_CMD  generate the sbatch command (in char format) to submit a SLURM
% job with specified parameters.
%--------------------------------------------------------------------------
% Inputs:
%  N_trim - vector of # of tasks in each job, needs a better name
%  jid - job ID number
%  cores - number of cores or processors to request from the job
%  nworkers - number of workers to request in the parallel pool. Keep in
%  mind that SPMD is disabled in submit.sh.
%  walltimes - list of walltimes to be indexed into by jid
%  mem - total memory allocated for the job in MB
%  qosopt - "Quality of Service" option: 'test', 'standby', ''
%  scriptfpath - filepath of the script to be run
%  dirpath - working directory to change to before executing the script.
%
% Outputs:
%  cmd - the sbatch command in char format
%
% Usage:
%  cmd = get_cmd(N_trim,jid,cores,nworkers,walltimes,mem,qosopt,scriptfpath,dirpath)
%
% Dependencies:
%  SLURM workload manager
%
% Author(s): Sterling Baird
%
% Date: 2020-09-15
%--------------------------------------------------------------------------
%set environment variables for UNIX shell
setenv('jid',sprintf(int2str(jid)))
disp(['jid = ' getenv('jid')])

walltimestr = ['walltime',int2str(jid)];
setenv(sprintf(walltimestr),sprintf(int2str(walltimes(jid))))
disp([walltimestr ' = ' getenv(walltimestr)])

walltime = walltimes(jid);
tidstr_list = ['1-',int2str(N_trim(jid))];

setenv('nworkers',sprintf(int2str(nworkers)))
%----------sbatch options-----------
arrayopt = ['sbatch -a ' tidstr_list]; %consider implementing tid_list batching
coreopt = ['-c ',int2str(cores)];
timeopt = ['--time=',int2str(walltime),':00']; %minutes
if ~isempty(qosopt)
    qosopt = ['--qos=' qosopt];
else
    qosopt = '';
end
imageopt = ['']; %[' -C rhel7 '];
jobnameopt = ['--job-name=batch',int2str(jid)];
memopt = ['--mem=',int2str(mem)];
diropt = ['-D ' dirpath];
scriptopt = [ scriptfpath ];
bckgdopt = ''; %[' &'];

%-------submit sbatch script--------
cmd = strjoin({arrayopt,coreopt,qosopt,timeopt,...
	imageopt,jobnameopt,memopt,diropt,scriptopt,bckgdopt},' ');
end %get_cmd
