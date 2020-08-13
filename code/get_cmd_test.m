%get_cmd test
%N_trim = [10 20 15 11];
N_trim = [1];
jid = 1;
cores = 8;
setenv('cores',int2str(cores))
walltimes = repelem(60,length(N_trim));
mem = 1024*6; %MB
qosopt = 'test'; %'', 'test'
script_fpath = 'submit.sh';
dirpath = './';

cmd = get_cmd(N_trim,jid,cores,walltimes,mem,qosopt,script_fpath,dirpath)

system(cmd)
