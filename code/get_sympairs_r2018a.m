function Spairs = get_sympairs_r2018a(pgnum,dispnameQ)
% arguments
% 	pgnum(1,1) int8 {mustBeInteger} = 32 %default == cubic
% 	dispnameQ(1,1) logical = false
% end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: get all combinations (pairs) of operators for a point group
% 
% Inputs:
%		pgnum - point group number, 32 == Oh (cubic)
%
%		dispnameQ - whether or not to display the name of the point group
%
% Outputs:
%		Spairs - combinations of symmetry operators, including repeats
%
% Usage:
%		Spairs = get_sympairs(32); % combinations of cubic symmetry operators
%
% Dependencies:
%		Pgsymops.mat, PGnames.mat (optional, default is not required)
%
%--------------------------------------------------------------------------

%load operators
symops = load('PGsymops.mat');

if dispnameQ
	%display point group name
	symnames = load('PGnames.mat');
	pgname = symnames.PG_names{pgnum};
	disp(pgname)
end

%unpack point group
qpt = symops.Q{pgnum};

%get all combinations of symmetry operators
qpttmp = num2cell(qpt,2);
qptpairs = allcomb(qpttmp,qpttmp);

% unpack symmetry operator combinations
SAlist = vertcat(qptpairs{:,1});
SBlist = vertcat(qptpairs{:,2});

% catenate combinations
Spairs = [SAlist SBlist];

end %get_sympairs