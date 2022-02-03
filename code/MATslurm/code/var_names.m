function out = var_names(varargin)
% VAR_NAMES  take variables and output a combined struct with each of the variable names
%--------------------------------------------------------------------------
% https://www.mathworks.com/matlabcentral/answers/79281#answer_89015 
%--------------------------------------------------------------------------
for n = 1:nargin
	out.(inputname(n)) = varargin{n};
end
end