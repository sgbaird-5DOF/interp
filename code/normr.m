function nrlist = normr(rlist)
% NORMR  normalizes vectors row-by-row. Outputs zero vector if zero vector input (shadowed by Computer Vision Toolbox)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-03
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

%determine size
[n,d] = size(rlist);

% shortcut if one row
if n == 1
	nm = norm(rlist);
	if nm ~= 0
		nrlist = rlist./norm(rlist);
	else
		nrlist = zeros(n,d);
	end
else
	%initialize
	nrlist = zeros(size(rlist));
	
	%compute norms
	nmlist = vecnorm(rlist,2,2);
	
	%get indices of non-zero elements
	ids = find(nmlist);
	
	%normalize only rows where norm is non-zero
	nrlist(ids,:) = rlist(ids,:)./nmlist(ids);
	
	%note: when nm(id) == 0, nrlist(id) == zeros(1,d)
end

end %normr


%-----------------------------CODE GRAVEYARD-------------------------------
%{
for i = 1:npts
		r = rlist(i,:);
		
		nm = nmlist(i);
		if nm ~= 0
			nr = r./nm;
		else
			nr = zeros(size(r));
		end
		
		nrlist(i,:) = nr;
		
	end

% 	ids = nmlist ~= 0;
%}