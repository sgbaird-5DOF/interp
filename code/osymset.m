function symocts = osymset(qA,qB,Spairs,pgnum,grainexchangeQ,doublecoverQ,uniqueQ,epsijk)
arguments
    qA(1,4) double {mustBeNumeric,mustBeFinite}
    qB(1,4) double {mustBeNumeric,mustBeFinite}
    Spairs(:,8) double {mustBeNumeric} = get_sympairs(32,false) %default to cubic Oh symmetry
    pgnum(1,1) double {mustBeInteger} = 32 % default == cubic Oh point group
    grainexchangeQ(1,1) logical {mustBeLogical} = true
    doublecoverQ(1,1) logical {mustBeLogical} = true
    uniqueQ(1,1) logical {mustBeLogical} = false
    epsijk(1,1) double {mustBeInteger} = 1
end
% OSYMSET  get symmetrically equivalent octonions
%--------------------------------------------------------------------------
% Author: Sterling Baird / Oliver Johnson
%
% Date: 2020-07-27
%
% Inputs:
%		(qA, qB) - quaternions
%		Spairs - list of pairs of symmetry operators to be applied to qA and qB
%       pgnum - point group number (1-32)
%       grainexchangeQ - logical indicating whether or not to apply grain
%                        exchange symmetry
%       doublecoverQ - logical indicating whether or not to apply double
%                      cover symmetry
%       uniqueQ - logical indicating whether or not to reduce output to
%                 only the unique octonions
%       epsijk - (1) or (-1) indicating the active or passive convention,
%                respectively
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
% Notes:
%  Could be sped up by doing multiple qA/qB pairs at a time instead of a
%  single qA/qB pair (i.e. batching/vectorizing approach). Would need to
%  pay attention to stacking order and perhaps better to output as a cell
%  instead of an array.
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
qSA = qmult(qArep,SAlist,epsijk);
qSB = qmult(qBrep,SBlist,epsijk);

% apply grain exchange and/or double cover symmetries
qXpi = repmat([0 1 0 0],nsyms,1); % rotation by pi around the x axis used for grain exchange symmetry
if grainexchangeQ && doublecoverQ
    symocts = [...
        qSA     qSB
        qSA     -qSB
        -qSA    qSB
        -qSA    -qSB
        qmult(qXpi,qSB,epsijk)     qmult(qXpi,qSA,epsijk)
        qmult(qXpi,qSB,epsijk)     qmult(qXpi,-qSA,epsijk)
        qmult(qXpi,-qSB,epsijk)	   qmult(qXpi,qSA,epsijk)
        qmult(qXpi,-qSB,epsijk)	   qmult(qXpi,-qSA,epsijk)];
    
elseif grainexchangeQ && ~doublecoverQ
    symocts = [...
        qSA qSB
        qmult(qXpi,qSB,epsijk) qmult(qXpi,qSA,epsijk)];
    
elseif ~grainexchangeQ && doublecoverQ
    symocts = [...
        qSA qSB
        qSA -qSB
        -qSA qSB
        -qSA -qSB];
    
elseif ~(grainexchangeQ || doublecoverQ)
    symocts = [...
        qSA qSB];
end

% apply inversion symmetry if point group is a Laue (centrosymmetric) group
% Laue Groups | pgnum  | name
LaueGroupNum = [2,...  % -1
                5,...  % 2/m
                8,...  % mmm
                11,... % 4/m
                15,... % 4/mmm
                17,... % -3
                20,... % -3m
                23,... % 6/m
                27,... % 6/mmm
                29,... % m-3
                32];   % m-3m
if ismember(pgnum,LaueGroupNum)
    qXpi = repmat([0 1 0 0],size(symocts,1),1); % rotation by pi around the x axis
    symocts = [...
        symocts
        qmult(qXpi,symocts(:,1:4),epsijk)    qmult(qXpi,symocts(:,5:8),epsijk)]; % this enforces [M,n]~[M,-n]
end

%reduce to unique set of octonions
if uniqueQ
    symocts = uniquetol(round(symocts,12),'ByRows',true);
end

end %osymset

%% CODE GRAVEYARD
%{
%following seems to produce inconsistent results in VFZ workflow:
% qSA = qmult(SAlist,qArep,epsijk);
% qSB = qmult(SBlist,qBrep,epsijk);
%}
