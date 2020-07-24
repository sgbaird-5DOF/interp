function cmd = get_cmd(N_trim,jid,cores,walltimes,mem,qosopt,script_fpath,dirpath)
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
%set environment variables for UNIX shell
setenv('jid',sprintf(int2str(jid)))
disp(['jid = ' getenv('jid')])

walltimestr = ['walltime',int2str(jid)];
setenv(sprintf(walltimestr),sprintf(int2str(walltimes(jid))))
disp([walltimestr ' = ' getenv(walltimestr)])

walltime = walltimes(jid);
tidstr_list = ['1-',int2str(N_trim(jid))];

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
memopt = ['--mem-per-cpu=',int2str(mem)];
diropt = ['-D ' dirpath];
scriptopt = [ script_fpath ];
bckgdopt = ''; %[' &'];

%-------submit sbatch script--------
cmd = strjoin({arrayopt,coreopt,qosopt,timeopt,...
	imageopt,jobnameopt,memopt,diropt,scriptopt,bckgdopt},' ');
end %get_cmd
