function repsets = get_repsets(pts)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-06-30
%
% Description: Find sets of non-unique (i.e. degenerate) points 
%
% Inputs:
%		pts				===	set of points (rows) that may or may not contain
%									duplicate values
% 
% Outputs:
%		irepset			===	cell array that contains sets of indices
%									of a unique point in pts. 
%
% % Example:
% pts = [...
% 	0.7071    0.7071    0.0000
% 	0.8629    0.5054    0.0000
% 	0.9675    0.2527    0.0000
% 	1.0000   -0.0000    0.0000
% 	0.6325    0.6325    0.4472
% 	0.8125    0.3366    0.4760
% 	0.8944   -0.0000    0.4472
% 	0.3162    0.3162    0.8944
% 	0.4472   -0.0000    0.8944
% 	-0.0000        0    1.0000
% 	0.7071    0.7071    0.0000
% 	0.8629    0.5054    0.0000
% 	0.9675    0.2527    0.0000 ];
%
% output: irepset{:} = [1 11], [2 12], [3 13].
%
% Note that pts(irepset{i},:) gives the duplicate values. Set uniquetype ==
% 'tol' if you want to avoid issues with numerical precision.
%--------------------------------------------------------------------------

%find index lists from unique()
uniquetype = 'tol'; %'tol','std' as in standard

switch uniquetype
	case 'std'
		[~,ia,ic] = unique(pts,'rows');
	case 'tol'
		tol = 1e-6;
		[~,ia,ic] = uniquetol(pts,tol,'ByRows',true);
end

iaNum = numel(ia);
[N, ~, icount] = histcounts(ic,iaNum);

irep = find(N(icount) >= 1);

% extract the cases
ic_rep = icount(irep);

ic_unique = sort(unique(ic_rep),1);
iculength = length(ic_unique);
icset = cell(1,iculength);
repsets = icset;

for i = 1:iculength	
	%take a single microstructure that corresponds to a non-unique
	%microstructure set
	ic_val = ic_unique(i);
	
	% find the index of that microstructure in the set of all non-unique
	% microstructures
	icset{i} = find(ic_rep == ic_val);
	
	% correlate that index back with the master list of microstructures
	repsets{i} = irep(icset{i});
end

%-----------------------------CODE GRAVEYARD-------------------------------
%{

if ~includeUniqueQ
	irep = find(N(icount) > 1);
else
	irep = find(N(icount) >= 1);
end
%}

%{
includeUniqueQ	===	whether or not to include unique sets as 1x1 %
numerics within an irepset cell.

Also includes indices of unique points as their own cell if includeUniqueQ
== true
%}