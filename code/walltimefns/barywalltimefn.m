function walltime = barywalltimefn(ndatapts,npredpts,cores,datatype)
arguments
   ndatapts(1,1) double
   npredpts(1,1) double
   cores(1,1) double
   datatype char {mustBeMember(datatype,{'brk','kim','rohrer-Ni','rohrer-test','olmsted-Ni'})} = 'brk'
end
switch ndatapts
    case num2cell(1:1000)
        walltime0 = 5*cores;
    case num2cell(1001:5000)
        walltime0 = 15*cores;
    case num2cell(5001:10000)
        walltime0 = 30*cores;
    case num2cell(10001:50000)
        walltime0 = 90*cores;
    case num2cell(50001:100000)
        walltime0 = 300*cores;
    otherwise
        error('too many datapts, update barywalltimefn')
end

scale = 1;
walltime1 = scale*npredpts/cores; % scale was ~0.7 for 1000 mesh points, 500 data points, 6 cores

walltime = (walltime0 + walltime1)/cores;

if strcmp(datatype,'kim')
    walltime = 2*walltime;
end

end

%% CODE GRAVEYARD
%{
% switch npredpts
%     case num2cell(1:1000)
%         walltime1 = 100;
%     case num2cell(1001:5000)
%         walltime1 = 200;
%     case num2cell(5001:10000)
%         walltime1 = 400;
%     case num2cell(10001:50000)
%         walltime1 = 800;
%     otherwise
%         error('too many predpts, update barywalltimefn')
% end
%}
