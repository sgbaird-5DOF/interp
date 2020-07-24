%% Variable Names

ninf = 1;
rand_id = 1348;
load_folder = 'parameters';
load_filename = ['parameters_ninf',int2str(ninf),'_rid',int2str(rand_id),'.mat'];
load_filepath = ['"',fullfile(load_folder,load_filename),'"'];

njobs = 3;
jbinfo.N_trim = [randi(100,1,njobs)];
jbinfo.walltimes = repelem(10,njobs); %minutes
jbinfo.failedjobQ = false;
jbinfo.failed_filepath = '';
jbinfo.env = 'slurm';
jbinfo.for_type = 'single';
jbinfo.load_filepath = load_filepath;
jbinfo.timefactor = 3; % number of times to multiply previous timefactor by (to help ensure that job finishes)

sbopts.cores = 2;
sbopts.mem = 1024*4;
sbopts.qosopt = '';
sbopts.dirpath = '~/compute/1DOF_functional';
sbopts.script_fpath = ' ~/compute/1DOF_functional/infer.sh';

submit_sbatch(jbinfo,sbopts)


