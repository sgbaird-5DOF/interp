function olist = osymsets(data,pgnum,varargin)
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-15
%
% Description: Get symmetrically equivalent octonions (osymsets) to a list
%					of input octonions
% 
% Inputs:
%		data	===	rows of octonions
%
%		pgnum ===	point group number (e.g. 32 == cubic)
%
% Outputs:
%		olist ===	1*size(data,1) cell array containing rows of
%						symmetrically equivalent octonions
%
% Dependencies:
%
%		Pgnames.mat, PGsymops.mat
%
% Notes:
%		Adapted a portion from Grain Boundary Octonion function GBdist.m from
%		Elizabeth Holm's CMU group.
%--------------------------------------------------------------------------
if nargin - 2 == 1
	usv = varargin{1};
else
	usv = [];
end

%% load crystal symmetry
symnames = load('PGnames.mat'); %need to add crystal_symmetry_ops to path in order for this to work
symops = load('PGsymops.mat');

pgname = symnames.PG_names{pgnum};

%unpack point group
qpt = symops.Q{pgnum}; %choose point group symmetry
npt = size(qpt,1);

%% perform distance calculation

ndatapts = size(data,1);
% disp('Performing distance calculation ...')

if size(data,2) == 7 && ~isempty(usv)
	data = proj_up(data,usv);
else
	data = [data zeros(size(data,1),1)];
end

olist = cell(1,ndatapts);
for i = 1:ndatapts
	%unpack quaternions
	qA = normr(data(i,1:4));
	qB = normr(data(i,5:8));
	
	%loop through symmetries
	%initialize
	olist{i} = zeros(npt^2*8,8);
	m = 0;
	for j = 1:npt
		for k = 1:npt
			m = m+1; %counter
			nsyms = 8;
			n = nsyms*m;
			
			%define symmetry operators
			Si = qpt(j,:);
			Sj = qpt(k,:);
			
			%apply symmetry operators
			qSA = qmult(Si,qA);
			qSB = qmult(Sj,qB);
			
			%catenate grain exchange & double cover
			octs = [...
				 qSA	 qSB
				 qSA	-qSB
				-qSA	 qSB
				-qSA	-qSB
				 qSB	 qSA
				 qSB	-qSA
				-qSB	 qSA
				-qSB	-qSA];
			
			olist{i}(n-nsyms+1:n,:) = octs;
		end
	end
	olist{i} = uniquetol(round(olist{i},12),'ByRows',true);
end


%-----------------------CODE GRAVEYARD-------------------------------------
%{
			% 			%double cover
			% 			qSC = qmult(Si,-qA);
			% 			qSD = qmult(Sj,-qB);

%}