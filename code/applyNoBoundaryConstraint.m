function [ypred,ysd,yint] = applyNoBoundaryConstraint(ypred,ysd,yint,qm2,mdl,nv)
arguments
    ypred(:,1) double = []
    ysd(:,1) double = []
    yint(:,2) double = []
    qm2(:,4) = [] % query misorientations
    mdl = [] % trained model 
    nv.weights char {mustBeMember(nv.weights,{'rsw','kernel'})} = 'rsw'
    nv.wthreshold(1,1) double = deg2rad(5)
end

% APPLYNOBOUNDARYCONSTRAINT  Adjusts model predictions to enforce the
% no-boundary constraint by applying a weighting function eta(w) in [0,1],
% which is a function of the disorientation angle, w.
%--------------------------------------------------------------------------
% Inputs:
%  ypred - interpolated property values of queried grain boundaries (the
%          mean model of the posterior).
%
%  ysd - prediction uncertainty specified as one standard deviation for the
%        mean model of the posterior.
%
%  yint - prediction uncertainty specified as the 95% confidence interval
%         around the mean model of the posterior.
%
%  qm2 - list of misorientation quaternions of the query points, as in 
%        qmult(qB,qinv(qA)) for grains A and B, where qA and qB are in the 
%        sample frame
%
%  mdl - trained model (a GPR model or compact GPR model), or []
%
%  nv - method-specific name-value pairs
%       'weights' - string ('rsw', or 'kernel') indicating which weighting
%                   scheme to use. 'rsw' weighting uses a 
%                   Read-Shockley-Wolf (RSW) type function [1-2]. The 
%                   'kernel'  weighting scheme uses the same correlation 
%                   function  contained in the trained model (mdl). The 
%                   'kernel' weighting scheme is NOT recommended as it 
%                   incorrectly suppresses predicted values far from the 
%                   no-boundary  domain. A warning is issued if the user 
%                   specifies the 'kernel' weighting scheme.
%
%       'wthreshold' - a scalar indicating the disorientation angle
%                      threshold (in radians) to use with the 'rsw'
%                      weighting scheme. When using 'rsw' weights, the
%                      weighting function, eta(w), will go smoothly from 0
%                      at w = 0 to 1 at w = wthreshold. NOTE: the parameter
%                      'a' in the RSW function controls the slope. A value
%                      of a = 1 (the default) will make the slope equal to
%                      0 at w = wthreshold.
%
% Outputs:
%  ypred - interpolated property values of queried grain boundaries (the
%          mean model of the posterior), updated after applying the
%          no-boundary constraint.
%
%  ysd - prediction uncertainty specified as one standard deviation for the
%        mean model of the posterior, updated after applying the
%          no-boundary constraint.
%
%  yint - prediction uncertainty specified as the 95% confidence interval
%         around the mean model of the posterior, updated after applying the
%          no-boundary constraint.
%
% Usage:
%  [ypred,ysd,yint] = applyNoBoundaryConstraint(ypred,ysd,yint,qm2,mdl);
%  [ypred,ysd,yint] = applyNoBoundaryConstraint(ypred,ysd,yint,qm2,mdl,'weights','rsw');
%  [ypred,ysd,yint] = applyNoBoundaryConstraint(ypred,ysd,yint,qm2,mdl,'weights','rsw','wthreshold',deg2rad(5));
%
% Dependencies:
%  MATLAB 2019b or higher (mainly for the "arguments" syntax checking at
%  the beginning of functions, which is used extensively throughout)
%
% Author(s): Oliver Johnson
%
% Date: 2022-03-23
%
% [1] Bulatov, V. V, Reed, B. W., & Kumar, M. (2014). Grain boundary energy
%     function for fcc metals. Acta Materialia, 65, 161–175. 
%     https://doi.org/10.1016/j.actamat.2013.10.057
% [2] Wolf, D. (1989). A read-shockley model for high-angle grain 
%     boundaries. Scripta Metallurgica, 23(10), 1713–1718. 
%     https://doi.org/10.1016/0036-9748(89)90348-7
%--------------------------------------------------------------------------

% compute disorientation angles
qd = disorientation(qm2,'cubic'); % there may be a faster way to do this since we just want the angle
[w,~,~] = q2rot(qd);

% assign weights
eta = ones(size(w));
switch nv.weights
    case 'rsw'

        eta(w <= nv.wthreshold) = RSWfun(w(w <= nv.wthreshold),0,nv.wthreshold,0,1,1);

    case 'kernel'
        
        str = input('The ''kernel'' weighting scheme is not recommended as it excessively suppresses the predictions far from the no-boundary domain. Use at your own risk. Continue? (y/n):','s');
        
        switch lower(str)
            case 'y'
        
                dE = w/4; % need distance in dE rad = 0.25*w, since model uses dE rad (and hence kernel parameters are in dE rad)
                switch lower(mdl.KernelInformation.Name)
                    case 'exponential'
                        eta = 1-exp(-dE/mdl.KernelInformation.KernelParameters(1));
                    case 'squaredexponential'
                        eta = 1-exp(-0.5*(dE/mdl.KernelInformation.KernelParameters(1)).^2);
                end

            case 'n'

                disp('applyNoBoundaryConstraint cancelled.')
                return

            otherwise

                error('Unrecognized input.')

        end
end

% update predictions
ypred = eta.*ypred;
if nargout > 1
    ysd = eta.*ysd;
    yint = eta.*yint;
end

end

%--------------------------Helper Functions-------------------------------%

function fout = RSWfun(w,wmin,wmax,ymin,ymax,a)

fout = zeros(size(w));
fout(w ~= wmin) = ((ymax-ymin)*sin((pi/2)*(w(w~=wmin)-wmin)/(wmax-wmin)).*(1-a*log(sin((pi/2)*(w(w~=wmin)-wmin)/(wmax-wmin))))+ymin);
fout(w == wmin) = ymin; % when w == wmin you get 0*inf = nan, so we have to correct for this

end