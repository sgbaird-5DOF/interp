function osets = osymsets_r2018a(oct,pgnum,usv,grainexchangeQ,doublecoverQ)
% arguments
% 	oct(:,8) double {mustBeFinite,mustBeReal,mustBeSqrt2Norm}
% 	pgnum(1,1) double {mustBeInteger} = 32 %default to Oh cubic point group
% 	usv(1,1) struct = struct
% 	grainexchangeQ logical = false
% 	doublecoverQ logical = false
% end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-15
%
% Description: Get symmetrically equivalent octonions (osymsets) for each
%					octonion in a list of octonions
%
% Inputs:
%		data	===	rows of octonions
%
%		pgnum ===	point group number (e.g. 32 == cubic)
%
%		usv	===	optional, use if "data" is 7D and needs to be projected
%						back to 8D
%
% Outputs:
%		olist ===	1*size(data,1) cell array containing rows of
%						unique symmetrically equivalent octonions
%
% Dependencies:
%		Pgnames.mat, PGsymops.mat
%
%		osymset.m
%			--qmult.m
%
% Notes:
%		Adapted a portion from Grain Boundary Octonion function GBdist.m from
%		Elizabeth Holm's CMU group github page.
%--------------------------------------------------------------------------

%% load symmetry operator combinations
Spairs = get_sympairs_r2018a(pgnum,false);

%% reformat data (if applicable)
ndatapts = size(oct,1);
if size(oct,2) == 7 && ~isempty(fields(usv))
	oct = proj_up(oct,usv);
elseif size(oct,2) == 7
	oct = [oct zeros(size(oct,1),1)];
end

%% get symmetrically equivalent octonions
%initialize
osets = cell(1,ndatapts);

%normalize quaternions
qAlist = normr(oct(:,1:4));
qBlist = normr(oct(:,5:8));

%loop through quaternion pairs
for i = 1:ndatapts
	%unpack quaternions
	qA = qAlist(i,:);
	qB = qBlist(i,:);
	
	%get symmetrically equivalent octonions
	osets{i} = osymset_r2018a(qA,qB,Spairs,grainexchangeQ,doublecoverQ);
end

end %osymsets

%-----------------------CODE GRAVEYARD-------------------------------------
%{
			% 			%double cover
			% 			qSC = qmult(Si,-qA);
			% 			qSD = qmult(Sj,-qB);



for i = 1:ndatapts
	%unpack quaternions
	qA = qAlist(i,:);
	qB = qBlist(i,:);
	
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




	%unpack quaternions
	qA = qAlist(i,:);
	qB = qBlist(i,:);
	
	%vertically stack copies of quaternions
	qArep = repmat(qA,nsyms);
	qBrep = repmat(qB,nsyms);
	
	%apply symmetry operators
	qSA = qmult(SA,qArep);
	qSB = qmult(SB,qBrep);
	
	% apply grain exchange & double cover
	symocts = [...
		qSA	 qSB
		qSA	-qSB
		-qSA	 qSB
		-qSA	-qSB
		qSB	 qSA
		qSB	-qSA
		-qSB	 qSA
		-qSB	-qSA];
	
	olist{i} = uniquetol(round(symocts,12),'ByRows',true);


symnames = load('PGnames.mat'); %need to add crystal_symmetry_ops to path in order for this to work
symops = load('PGsymops.mat');

pgname = symnames.PG_names{pgnum};

%unpack point group
qpt = symops.Q{pgnum}; %choose point group symmetry
npt = size(qpt,1);

%get all combinations of symmetry operators
qpttmp = num2cell(qpt,2);
qptpairs = allcomb(qpttmp,qpttmp);

% unpack symmetry operator combinations
SAlist = vertcat(qptpairs{:,1});
SBlist = vertcat(qptpairs{:,2});

%}