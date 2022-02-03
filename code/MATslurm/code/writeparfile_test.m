%writeparfile_test
nparvars = 10;
par1 = randi(10000,nparvars,1);
par2 = randi(10000,nparvars,1);
pars = struct('par1',par1,'par2',par2);

walltimefn = @(par1,par2) 20*rand();

addpathdir('allcomb.m')

[fpath, Nidxlength, Ntrim, jobwalltimes] = writeparfile(pars,walltimefn,'1');