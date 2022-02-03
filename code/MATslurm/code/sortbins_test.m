%sortbins_test
npts = 1000;
walltimes = 10*rand(npts,1);
par1 = randi(10000,npts,1);
parstruct = struct('walltimes',num2cell(walltimes),'par1',num2cell(par1));

[Ntrim, Nidxlength, parcombsets, jobwalltimes] = ...
    sortbins(parstruct, walltimes);