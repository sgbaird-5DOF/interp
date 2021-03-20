function walltimefn = get_walltimefn(ndatapts,npredpts,method,cores,datatype,K) %#ok<*INUSD,*INUSL>
% GET_WALLTIMEFN  get a function handle for computing interp5DOF SLURM walltimes (e.g. 'gpr')
arguments
    ndatapts(1,1) double
    npredpts(1,1) double
    method char {mustBeMember(method,{'sphgpr','gpr','pbary','sphbary',...
        'idw','nn','avg','egprm','egpr','gprm','gprmix'})}
    cores(1,1) double
    datatype char %{mustBeMember(datatype,{'brk','kim','rohrer-Ni',...
        %'rohrer-test','rohrer-brk-test','olmsted-Ni','olmsted-MgO'})} = 'brk'
    K(1,1) double = 1
end
walltimeBuffer = 5;

switch method
    case {'sphgpr','gpr','egprm','egpr','gprm','gprmix'}
        walltimefn = K*gprwalltimefn(ndatapts,cores) + walltimeBuffer;
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