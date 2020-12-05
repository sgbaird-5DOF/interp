function symocts = osymset(qA,qB,Spairs,grainexchangeQ,doublecoverQ,epsijk)
arguments
    qA(1,4) double {mustBeNumeric,mustBeFinite}
    qB(1,4) double {mustBeNumeric,mustBeFinite}
    Spairs(:,8) double {mustBeNumeric} = get_sympairs(qA,qB)
    grainexchangeQ logical = false
    doublecoverQ logical = false
    epsijk = 1
end
% OSYMSET  get symmetrically equivalent octonions
%--------------------------------------------------------------------------
% Author: Sterling Baird
%
% Date: 2020-07-27
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
qArep = repmat(qA,nsyms,1);
qBrep = repmat(qB,nsyms,1);

%unpack pairs
SAlist = Spairs(:,1:4);
SBlist = Spairs(:,5:8);

%apply symmetry operators
qSA = qmult(SAlist,qArep,epsijk);
qSB = qmult(SBlist,qBrep,epsijk);

if grainexchangeQ && doublecoverQ
    %apply grain exchange & double cover
    symocts = [...
        qSA     qSB
        qSA     -qSB
        -qSA    qSB
        -qSA    -qSB
        qSB     qSA
        qSB     -qSA
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