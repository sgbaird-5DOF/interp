function fname = get_fname(sampleMethod,opts)
arguments
	sampleMethod char
	opts struct
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-09
%
% Description: Get a filename to load or save
%
% Inputs:
%		sampleMethod	=== the sampling method used in datagen
%
%		res				=== as in resDegrees for meshFZ.m
%
%		nint				=== for meshBP.m
%
% Outputs:
%		fname				=== filename
%
% Notes:
%		example outputs
%			'5DOF_exterior_res15_nint2_octsubdiv1.mat'
%
%			'5DOF_vtx_octsubdiv2.mat'
%
%			'Rohrer2009_octsubdiv2.mat'
%
% Dependencies:
%
%--------------------------------------------------------------------------


if any(cellfun(@(pattern) contains(sampleMethod,pattern),{'resolution','interior','exterior'}))
	
	fname = [sampleMethod '_res' int2str(opts.res) ...
		'_nint' int2str(opts.nint) '_octsubdiv' int2str(opts.octsubdiv) '.mat'];
	
elseif contains(sampleMethod,'ocubo')
	% ocuboOpts fields: n, method, sidelength, seed (2020-07-27)
	varnames = fields(opts.ocuboOpts);
	fname = sampleMethod;
	for i = 1:length(varnames)
		varname = varnames{i};
		var = opts.ocuboOpts.(varname);
		if ~isempty(var)
			if ischar(var)
				fname = [fname '_' varname var];
			else
				fname = [fname '_' varname int2str(var)];
			end
		end
	end
	fname = [fname '_octsubdiv' int2str(opts.octsubdiv) '.mat'];
	
else
	fname = [sampleMethod '_octsubdiv' int2str(opts.octsubdiv) '.mat'];
end
end


%----------------------------CODE GRAVEYARD--------------------------------
%{

switch sampleMethod
	case {'5DOF_resolution','5DOF_interior','5DOF_exterior'}
		fname = [sampleMethod '_res' int2str(res) ...
			'_nint' int2str(nint) '_octsubdiv' int2str(octsubdiv) '.mat'];
	otherwise
		fname = [sampleMethod '_octsubdiv' int2str(octsubdiv) '.mat'];
end

% if nargin > 1
% 	res = varargin{1};
% 	nint = varargin{2};
% 	octsubdiv = varargin{3};
% 	ocuboOpts = varargin{4};
% end

	n = opts.ocuboOpts.n;
	method = opts.ocuboOpts.method;
	sidelength = opts.ocuboOpts.sidelength;
	seed = opts.ocubo

%}