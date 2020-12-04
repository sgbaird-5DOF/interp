function [qm,nA] = eumA2five(eA,eB,mA,epsijk)
arguments
    eA(:,3) double {mustBeFinite,mustBeReal}
    eB(:,3) double {mustBeFinite,mustBeReal}
    mA(:,3) double {mustBeFinite,mustBeReal}
    epsijk(1,1) double = 1
end
% EUMA2FIVE Convert sample frame euler angles of grain A and grain B and
% sample frame boundary plane normal pointing outward from grain A towards
% grain B to misorientation quaternion and boundary plane normal (crystal
% reference frame of grain A).
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-08-22
% 
% Inputs:
%  eA, eB - Euler angles of grains A and B in sample reference frame,
%    resp.
%  mA - boundary plane normal (sample reference frame) pointing from grain
%  A to grain B
%
% Outputs:
%  five -struct containing qm and nA
%   --qm - misorientation quaternion (qB qA*)
%
%   --nA - boundary plane normal in grain A crystal frame vec[qA mA qA*]
%
% Usage:
%  five = eumA2qnA(euA,euB,mA)
%
% Dependencies:
%  qmult.m
%  qinv_francis.m
%  eu2qu.m
%
% References:
%  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
%
%  [2]
%  Seita, M.; Volpi, M.; Patala, S.; McCue, I.; Schuh, C. A.; Diamanti, M.
%  V.; Erlebacher, J.; Demkowicz, M. J. A High-Throughput Technique for
%  Determining Grain Boundary Character Non-Destructively in
%  Microstructures with through-Thickness Grains. Npj Computational
%  Materials 2016, 2, 16016.

%--------------------------------------------------------------------------
qA = eu2qu(eA,epsijk); % grain A, sample frame
qB = eu2qu(eB,epsijk); % grain B, sample frame

[qm,nA] = qmA2five(qA,qB,mA,epsijk);

%----------------------------CODE GRAVEYARD--------------------------------
%{
setGlobal_epsijk(1); %Morawiec convention seems to be consistent

    qmconvention char {mustBeMember(qmconvention,{'francis','johnson'})} = 'johnson'
qm = qlab2qm(qA,qB,qmconvention);

% qm = qmult(qB,qinv_francis(qA)); %ref [1], eqn. S-4, francis convention
% qm = qmult(qinv_francis(qA),qB); %johnson convention

%covert BP normal from sample frame to grain 1 crystal frame: undo equation
%S-5 of [1], such that mA = vec[pA mAhat pA*], except in the notation I'm
%using, nA = vec[qA mA qA*]
% Zero = zeros(size(mA,1),1); %column of zeros
% nAtmp = qmult(qA,qmult([Zero mA],qinv(qA)));
% nA = nAtmp(:,2:4);


% oA = eu2om(eA);
[nA,nB] = deal(zeros(npts,3));
for i = 1:npts
%     nA(i,:) = Lpr(qinv(qA(i,:)),mA(i,:),epsijk);
    nA(i,:) = Lpr(qA(i,:),mA(i,:),epsijk); %sample frame to crystal frame?
    nB(i,:) = Lpr(qB(i,:),mA(i,:),epsijk); %sample frame to crystal frame?
%     nA(i,:) = (oA{i}.'\(mA(i,:).')).'; % eq(3) from [2]
%     nA(i,:) = oA{i}(:,3).';
end

% nA = Lpr(qA,mA,epsijk);


npts = size(eA,1);
qm = qlab2qm(qA,qB,epsijk);

%boundary plane normal
[nA,nB] = deal(zeros(npts,3));
for i = 1:npts
    nA(i,:) = Lpr(qA(i,:),mA(i,:),epsijk); %sample frame to crystal frame?
%     nB(i,:) = Lpr(qB(i,:),mA(i,:),epsijk); %sample frame to crystal frame?
end

% nA = Lpr(qA,mA,epsijk);

%}