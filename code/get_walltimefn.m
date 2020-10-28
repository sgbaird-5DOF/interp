function walltimefn = get_walltimefn(ndatapts,npredpts,method,cores,datatype) %#ok<*INUSD,*INUSL>
arguments
    ndatapts(1,1) double
    npredpts(1,1) double
    method char {mustBeMember(method,{'sphgpr','gpr','pbary','sphbary','idw','nn','avg'})}
    cores(1,1) double
    datatype char {mustBeMember(datatype,{'brk','kim','rohrer-Ni'})} = 'brk'
end
walltimeBuffer = 5;

switch method
    case {'sphgpr','gpr'}
        walltimefn = gprwalltimefn(ndatapts,cores) + walltimeBuffer;
    case {'pbary','sphbary'}
        walltimefn = barywalltimefn(ndatapts,npredpts,cores,datatype) + walltimeBuffer;
    case {'idw'}
        walltimefn = walltimeBuffer*2;
    case {'nn'}
        walltimefn = walltimeBuffer;
    case {'avg'}
        walltimefn = walltimeBuffer;
end

end