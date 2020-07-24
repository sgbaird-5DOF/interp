function out = var_names(varargin)
%--------------------------------------------------------------------------
% https://www.mathworks.com/matlabcentral/answers/79281#answer_89015 
%--------------------------------------------------------------------------
for n = 1:nargin
	out.(inputname(n)) = varargin{n};
end
end