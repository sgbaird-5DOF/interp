function nrlist = normr(rlist)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
%
% Description: normalizes vectors row-by-row. Outputs a zero vector if a
%					zero vector is input.
% 
% Inputs:
%
%		rlist		===	rows of vectors to be normalized
%
% Outputs:
%
%		nrlist	===	normalized rows of vectors
%
% Dependencies:
%
% Note:	normr.m shadows a built-in function that's part of the Computer
%			Vision Toolbox
%
%--------------------------------------------------------------------------

% determine size & initalize array
if ~any(size(rlist) == 1)
	npts = size(rlist,1); %assume rows of points
	nr = zeros(npts,size(rlist,2));
else
	npts = 1;
	nr = zeros(1,npts);
end

%loop through points
nrlist = zeros(size(rlist));
for i = 1:npts
	r = rlist(i,:);
	
	nm = norm(r);
	if nm ~= 0
		nr = r./nm;
	else
		nr = zeros(size(r));
	end
	
	nrlist(i,:) = nr;
	
end

end %normr