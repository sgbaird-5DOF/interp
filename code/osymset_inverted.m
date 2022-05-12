function symocts = osymset(qA,qB,Spairs,grainexchangeQ,doublecoverQ,uniqueQ,epsijk)
arguments
    qA(1,4) double {mustBeNumeric,mustBeFinite}
    qB(1,4) double {mustBeNumeric,mustBeFinite}
    Spairs(:,8) double {mustBeNumeric} = get_sympairs(32,false) %default to cubic Oh symmetry
    grainexchangeQ(1,1) logical {mustBeLogical} = false
    doublecoverQ(1,1) logical {mustBeLogical} = false
    uniqueQ(1,1) logical {mustBeLogical} = false
    epsijk(1,1) double {mustBeInteger} = 1
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

if grainexchangeQ && doublecoverQ
    %apply grain exchange & double cover
    symocts_tmp = [...
        qA     qB
        qA     -qB
        -qA    qB
        -qA    -qB
        qB     qA
        qB     -qA
        -qB	 qA
        -qB	-qA];
    
symocts = cell(1,8);
for i = 1:size(symocts_tmp,1)
    symoct = symocts_tmp(i,:);
    qSA = qmult(qA, symoct(1:4), epsijk);
    qSB = qmult(qB, symoct(5:8), epsijk);
    symocts{i} = qmult(
end

% if grainexchangeQ && doublecoverQ
%     %apply grain exchange & double cover
%     symocts = [...
%         qSA     qSB
%         qSA     -qSB
%         -qSA    qSB
%         -qSA    -qSB
%         qSB     qSA
%         qSB     -qSA
%         -qSB	 qSA
%         -qSB	-qSA];
%     
% elseif grainexchangeQ && ~doublecoverQ
%     symocts = [...
%         qSA qSB
%         qSB qSA];
%     
% elseif ~grainexchangeQ && doublecoverQ
%     symocts = [...
%         qSA qSB
%         -qSA qSB
%         qSA -qSB
%         -qSA -qSB];
%     
% elseif ~(grainexchangeQ || doublecoverQ)
%     symocts = [...
%         qSA qSB];
% end

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
