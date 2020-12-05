function nA = qmA2nA(qA,mA,epsijk)
arguments
   qA(:,4) double
   mA(:,3) double
   epsijk(1,1) double
end
% QLAB2FIVE Convert lab/sample frame quaternions of grain A and grain B and
% and compute misorientation quaternion according to Toby Francis's or
% Oliver Johnson's convention
%--------------------------------------------------------------------------
% Author(s): Sterling Baird
% Date: 2020-08-22
% 
% Inputs:
%  qA - orientation of a grain in sample reference frame
%  mA - boundary plane normal in sample reference frame
%  convention - francis or johnson convention
%
% Outputs:
%  nA - boundary plane normal in crystal reference frame of q
%
% Usage:
%  nA = qmA2nA(q,mA,convention)
%
% Dependencies:
%  qmult.m
%  qinv_francis.m
%
% References:
%  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
%--------------------------------------------------------------------------

nA = Lpr(qA,mA,epsijk);

% npts = size(qA,1);
% Zero = zeros(npts,1);
% mA0 = [Zero mA];
% 
% nA = qmult(qinv(qA),qmult(mA0,qA,epsijk),epsijk);
% nA = nA(:,2:4);

end
%% CODE GRAVEYARD
%{
%    convention char {mustBeMember(convention,{'francis','johnson'})} = 'johnson'


npts = size(qA,1);
Zero = zeros(npts,1);
mA0 = [Zero mA];
switch convention
    case 'francis'
        nA = qmult(qinv_francis(qA),qmult(mA0,qA));
    case 'johnson'
        nA = qmult(qA,qmult(mA0,qinv_francis(qA)));
end
nA = nA(:,2:4);
%}