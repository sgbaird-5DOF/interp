function walltime = barywalltimefn(ndatapts,npredpts,cores)
switch ndatapts
    case num2cell(1:1000)
        walltime0 = 5*cores;
    case num2cell(1001:5000)
        walltime0 = 10*cores;
    case num2cell(5001:10000)
        walltime0 = 15*cores;
    case num2cell(10001:50000)
        walltime0 = 30*cores;
    otherwise
        error('too many datapts, update barywalltimefn')
end

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

scale = 2;
walltime1 = scale*cores/npredpts; % scale was ~0.7 for 1000 mesh points, 500 data points, 6 cores

walltime = (walltime0 + walltime1)/cores;

end
