function newpts = proj_up(pts,usv)
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
V = usv.V;
avg = usv.avg;

%lower dimension
d1 = size(pts,2);

%higher dimension
d2 = size(V,1);

ndegdim = d2-d1; %number of degenerate dimensions
% ndegdim = sum(abs(diag(S)) < 1E-6);

newpts = padarray(pts,[0 ndegdim],'post')*V'+avg;
% newpts = padarray(pts,[0 ndegdim],'post')*V';


end