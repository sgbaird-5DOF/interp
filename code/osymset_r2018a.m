function symocts = osymset_r2018a(qA,qB,Spairs,grainexchangeQ,doublecoverQ)
% arguments
% 	qA(1,4) double {mustBeNumeric,mustBeFinite}
% 	qB(1,4) double {mustBeNumeric,mustBeFinite}
% 	Spairs(:,8) double {mustBeNumeric} = get_sympairs(qA,qB)
% 	grainexchangeQ logical = false
% 	doublecoverQ logical = false
% end
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
%
% Description: get symmetrically equivalent octonions
%
% Inputs:
%		(qA, qB) - quaternions
%
%		Slist - list of pairs of symmetry operators to be applied to qA and qB
%
% Outputs: rows of symmetrically equivalent octonions
%
% Usage:
%			symocts = osymset(qA,qB); (calls get_sympairs once per function call)
%
%			symocts = osymset(qA,qB,Slist);
%
% Dependencies:
%		get_sympairs.m (optional, required if Slist not supplied)
%			--allcomb.m
%
%		qmult.m
%
%--------------------------------------------------------------------------
%number of symmetry operator pairs
nsyms = size(Spairs,1);

%vertically stack copies of quaternions
qArep = repmat(qA,nsyms);
qBrep = repmat(qB,nsyms);

%unpack pairs
SAlist = Spairs(:,1:4);
SBlist = Spairs(:,5:8);

%apply symmetry operators
qSA = qmult(SAlist,qArep);
qSB = qmult(SBlist,qBrep);

if grainexchangeQ && doublecoverQ
	%apply grain exchange & double cover
	symocts = [...
		qSA	 qSB
		qSA	-qSB
		-qSA	 qSB
		-qSA	-qSB
		qSB	 qSA
		qSB	-qSA
		-qSB	 qSA
		-qSB	-qSA];
	
elseif grainexchangeQ && ~doublecoverQ
	symocts = [...
		qSA qSB
		qSB qSA];
	
elseif ~grainexchangeQ && doublecoverQ
	symocts = [...
		qSA qSB
		-qSA qSB
		qSA -qSB
		-qSA -qSB];
	
elseif ~(grainexchangeQ || doublecoverQ)
	symocts = [...
		qSA qSB];
end

%reduce to unique set of octonions
symocts = uniquetol(round(symocts,12),'ByRows',true);

end %osymset