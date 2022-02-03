function submit_sbatch(parpath,cores,mem,qosopt,scriptfpath,NV)
arguments
    parpath = []
    cores = 24
    mem = 1024*4*cores %MB
    qosopt = 'test'
    scriptfpath = fullfile('MATslurm','code','submit.sh')
    NV.Ntrim = 1
    NV.jobwalltimes = 60
    NV.spmdQ = true
end
% SUBMIT_SBATCH submit multiple sbatch jobs that depend on a parameter file via writeparfile.m
%--------------------------------------------------------------------------
% Inputs:
%  parpath -
%  cores -
%  mem -
%  qosopt - 
%  scriptfpath - 
%  NV.Ntrim - 
%  NV.jobwalltimes - 
%
% Outputs:
%  b - b
%
% Usage:
%  a = b(a);
%
% Dependencies:
%  *
%
% Notes:
%  *
%
% see also WRITEPARFILE
%
% Author(s): Sterling Baird
%
% Date: 2020-09-08
%--------------------------------------------------------------------------
if ~isempty(parpath)
    load(parpath,'Ntrim','jobwalltimes')
else
    Ntrim = NV.Ntrim;
    jobwalltimes = NV.jobwalltimes;
end
Nidxlength = length(Ntrim); %number of jobs

% general environment variables
setenv('cores',int2str(cores))
setenv('spmdQ',int2str(NV.spmdQ))
setenv('parpath',['"',parpath,'"']);

for jid = 1:Nidxlength
    %set job-specific environment variables
%     setenv('jid',sprintf(int2str(jid)))
%     setenv(sprintf(['walltime',int2str(jid)]),sprintf(int2str(jobwalltimes(jid)))) %only for displaying on cmd line
    
    cmd = get_cmd(Ntrim,jid,cores,jobwalltimes,mem,qosopt,scriptfpath);
    disp(cmd)
    [status,cmdout] = system(cmd);
    disp(' ')
    disp(['status = ' int2str(status)])
    disp(['cmdout = ' cmdout])
end

end
