function errout = get_errmetrics(predicted,measured,type)
arguments
    predicted double {mustBeReal,mustBeFinite}
    measured double {mustBeReal,mustBeFinite}
    type char {mustBeMember(type,{'e','ae','mae','se','rmse','all'})} = 'all'
end
% GET_ERRMETRICS  get various error metrics for measured data relative to
% true data, and specify types based on the lowercase symbol. e.g. 'rmse',
% or 'all' to output a struct of error metrics. e = predicted - measured
%--------------------------------------------------------------------------
% Inputs:
%  predicted - scalar, vector, or matrix of predicted values
%
%  measured - scalar, vector, or matrix of measured (i.e. true or lookup)
%  values, matching size of predicted
%
%  type - type of error metric to output
%
% Outputs:
%  errout - the error metric (or a structure of error metrics) specified by
%  type, where e, ae, and se, are column vectors,
%
% Usage:
%  errmetrics = get_errmetrics(datatrue,dataout); %outputs struct
%
%  rmse = get_errmetrics(datatrue,dataout,'rmse'); %outputs rmse value
%
% Dependencies:
%  var_names.m
%
% Notes:
%  See https://en.wikipedia.org/wiki/Root-mean-square_deviation
%
% Author(s): Sterling Baird
%
% Date: 2020-09-17
%--------------------------------------------------------------------------

%additional argument validation
szpred = size(predicted);
szmeas = size(measured);
assert(all(szpred==szmeas),['predicted size: ' num2str(szpred), ', measured size: ' num2str(szmeas)])

%error metrics
e = predicted-measured; %error
ae = abs(e); %absolute error
mae = mean(ae,'all'); %mean absolute error
se = e.^2; %square error
mse = mean(se,'all'); %mean square error
rmse = sqrt(mse); %root mean square error

%compile into struct
errmetrics = var_names(e,ae,mae,se,mse,rmse);

%assign output error value(s)
switch type
    case 'all'
        errout = errmetrics;
    otherwise
        errout = errmetrics.(type); %dynamic field reference
end

end

%---------------------------CODE GRAVEYARD---------------------------------
%{
%assign output error value(s)
switch type
    case 'e'
        errout = e;
    case 'ae'
        errout = ae;
    case 'se'
        errout = se;
    case 'mse'
        errout = mse;
    case 'rmse'
        errout = rmse;
    case 'all'
        
end

%}