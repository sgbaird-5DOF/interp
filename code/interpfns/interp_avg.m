function [ypred,yavg] = interp_avg(propList,npredpts)
% INTERP_AVG assign average of input values as constant output model
yavg = mean(propList);
ypred = repelem(yavg,npredpts,1);
end