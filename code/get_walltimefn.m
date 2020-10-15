function walltimefn = get_walltimefn(ndatapts,npredpts,method,cores)

walltimeBuffer = 5;

switch method
    case {'sphgpr','gpr'}
        walltimefn = @(ndatapts,cores) gprwalltimefn(ndatapts,cores);
    case {'pbary','sphbary'}
        walltimefn = @(ndatapts,npredpts,cores) barywalltimefn(ndatapts,npredpts,cores);
    case {'nn'}
        walltimefn = @() walltimeBuffer;
    case {'avg'}
        walltimefn = @() walltimeBuffer;
end

end

function walltime = gprwalltimefn(ndatapts,cores)
switch ndatapts
    case num2cell(1:1000)
        walltime = 5;
    case num2cell(1001:5000)
        walltime = 8;
    case num2cell(5001:10000)
        walltime = 12;
    case num2cell(10001:50000)
        walltime = 20;
    otherwise
        error('too many datapts, update gprwalltimefn')
end
walltime = walltime/cores + walltimeBuffer;
end

function walltime = barywalltimefn(ndatapts,npredpts,cores)
switch ndatapts
    case num2cell(1:1000)
        walltime0 = 5;
    case num2cell(1001:5000)
        walltime0 = 8;
    case num2cell(5001:10000)
        walltime0 = 12;
    case num2cell(10001:50000)
        walltime0 = 20;
    otherwise
        error('too many datapts, update barywalltimefn')
end
walltime0 = walltime0 + walltimeBuffer;

switch npredpts
    case num2cell(1:1000)
        walltime1 = 5;
    case num2cell(1001:5000)
        walltime1 = 8;
    case num2cell(5001:10000)
        walltime1 = 12;
    case num2cell(10001:50000)
        walltime1 = 20;
    otherwise
        error('too many predpts, update barywalltimefn')
end
end
