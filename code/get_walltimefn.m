function walltimefn = get_walltimefn(ndatapts,npredpts,method,cores) %#ok<*INUSD,*INUSL>
arguments
    ndatapts(1,1) double
    npredpts(1,1) double
    method char {mustBeMember(method,{'sphgpr','gpr','pbary','sphbary','nn','avg'})}
    cores(1,1) double
end
walltimeBuffer = 5;

switch method
    case {'sphgpr','gpr'}
        walltime = gprwalltimefn(ndatapts,cores) + walltimeBuffer;
    case {'pbary','sphbary'}
        walltime = barywalltimefn(ndatapts,npredpts,cores) + walltimeBuffer;
    case {'nn'}
        walltime = walltimeBuffer;
    case {'avg'}
        walltime = walltimeBuffer;
end

end