function newpts = proj_up(pts,usv)
arguments
	pts double
	usv struct {mustContainFields(usv,{'V'})}
end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date:
%
% Description: project up (restore null dimensions)
% 
% Inputs:
%
% Outputs:
%
% Dependencies:
%
%--------------------------------------------------------------------------
% account for "zeropt" being in usv (see proj_down.m)
if isfield(usv,'zeropt')
	if ~isempty(usv.zeropt)
		pts = pts + usv.zeropt;
	end
end

V = usv.V;
avg = usv.avg;

%lower dimension
d1 = size(pts,2);

%higher dimension
d2 = size(V,2);

ndegdim = d2-d1; %number of degenerate dimensions
% ndegdim = sum(abs(diag(S)) < 1E-6);

newpts = padarray(pts,[0 ndegdim],'post')*V'+avg;
% newpts = padarray(pts,[0 ndegdim],'post')*V';


end