function walltime = gprwalltimefn(ndatapts,cores)
switch ndatapts
    case num2cell(1:1000)
        walltime = 5*cores;
    case num2cell(1001:5000)
        walltime = 8*cores;
    case num2cell(5001:10000)
        walltime = 12*cores;
    case num2cell(10001:50000)
        walltime = 60*cores;
    otherwise
        error('too many datapts, update gprwalltimefn')
end
walltime = walltime/cores;
end
