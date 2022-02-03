function cmd = get_cmd(Ntrim,jid,cores,jobwalltimes,mem,qosopt,scriptfpath,dirpath,NV)
arguments
   Ntrim double = 1
   jid double = 1
   cores = 24
   jobwalltimes double = 60
   mem(1,1) double = 1024*4 %MB
   qosopt char = 'test'
   scriptfpath char = fullfile('MATslurm','code','submit.sh')
   dirpath char = './' %default to just above MATslurm
   NV.requeueQ logical = true
end
% GET_CMD  get slurm command to submit an sbatch script
%--------------------------------------------------------------------------
% Inputs:
%  N_trim
%  jid
%  cores
%  walltimes
%  mem
%  qosopt
%  script_fpath
%  dirpath
%
% Outputs:
%
% Author: Sterling Baird
%
% Date: 2020-09-07
%--------------------------------------------------------------------------

%set environment variables for UNIX shell
setenv('jid',sprintf(int2str(jid)))
disp(['jid = ' getenv('jid')])

walltimestr = ['walltime',int2str(jid)];
setenv(sprintf(walltimestr),sprintf(int2str(jobwalltimes(jid))))
disp([walltimestr ' = ' getenv(walltimestr)])

walltime = jobwalltimes(jid);
tidstr_list = ['1-',int2str(Ntrim(jid))];

% setenv('cores',sprintf(int2str(cores)))
setenv('cores',sprintf(int2str(cores)))
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
memopt = ['--mem=',int2str(mem)]; %changed from --mem-per-cpu to --mem 2020-09-11 SGB
if NV.requeueQ
    reQopt = ['--requeue']; %'', '--requeue', e.g. requeue if preempted
else
    reQopt = [''];
end
diropt = ['-D ' dirpath];
scriptopt = [ scriptfpath ];
bckgdopt = ''; %[' &'];

%-------submit sbatch script--------
cmd = strjoin({arrayopt,coreopt,qosopt,timeopt,...
	imageopt,jobnameopt,memopt,reQopt,diropt,scriptopt,bckgdopt},' ');
end %get_cmd
