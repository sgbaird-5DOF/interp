function [ypred,ysd,yint] = predict5DOF(qm2,nA2,mdl,nv)
arguments
    qm2(:,4) = [] % query misorientations
    nA2(:,3) = [] % query BP normals
    mdl = [] % trained model 
    nv.epsijk(1,1) double = 1
    nv.noboundaryQ(1,1) logical {mustBeLogical} = true % enforce no-boundary constraint
    nv.weights char {mustBeMember(nv.weights,{'rsw','kernel'})} = 'rsw'
    nv.wthreshold = 'min'
    nv.qm(:,4) = [] % original input misorientation quaternions
end

% PREDICT5DOF  Predict GB property values from a trained model for GBs
% specified by the ordered pair of misorientations and GB normals (qm2,nA2)
%--------------------------------------------------------------------------
% Inputs:
%  qm2 - list of misorientation quaternions of the query points, as in 
%        qmult(qB,qinv(qA)) for grains A and B, where qA and qB are in the 
%        sample frame
%
%  nA2 - list of boundary plane Cartesian unit normals (grain A frame) of
%        the query points
%
%  mdl - trained model (a GPR model or compact GPR model)
%
%  nv - method-specific name-value pairs
%       'noboundarQ' - logical, whether or not to enforce the no boundary
%                      constraint (i.e. that the properties a strictly 
%                      constant when the disorientation angle is equal to 
%                      zero)
%       'weights' - string ('rsw', or 'kernel') indicating which weighting
%                   scheme to use. See applyNoBoundaryConstraint.m
%
%       'wthreshold' - a scalar indicating the disorientation angle
%                      threshold (in radians) to use with the 'rsw'
%                      weighting scheme; OR the string 'min', in which case
%                      wthreshold is set to the minimum of 
%                      [min(w),deg2rad(5)], where w contains the 
%                      disorientation angles of all of the input GBs
%                      (data). See applyNoBoundaryConstraint.m
%
% Outputs:
%  ypred - interpolated property values of queried grain boundaries (the
%          mean model of the posterior).
%
%  ysd - prediction uncertainty specified as one standard deviation for the
%        mean model of the posterior.
%
%  yint - prediction uncertainty specified as the 95% confidence interval
%         around the mean model of the posterior.
%
% Usage:
%  ypred = predict5DOF(qm2,nA2,mdl);
%  ypred = predict5DOF(qm2,nA2,mdl,'noboundaryQ',false);
%
%  [ypred,ysd,yint] = predict5DOF(qm2,nA2,mdl,'noboundaryQ',true);
%
% Dependencies:
%  MATLAB 2019b or higher (mainly for the "arguments" syntax checking at
%  the beginning of functions, which is used extensively throughout)
%
%  Toolboxes
%  -Statistics and Machine Learning Toolbox
%
% Author(s): Oliver Johnson, Sterling Baird
%
% Date: 2022-03-22
%--------------------------------------------------------------------------

%% Convert query (qm2,nA2) to proper coordinates

% convert query points to octonions
o = five2oct(qm2,nA2,nv.epsijk); % will have norm of sqrt(2)

% map to same VFZ as the model
o = get_octpairs(o,nv.epsijk,'oref',mdl.oref,'dispQ',false);

% convert to unit norm
o = normr(o);

% project to 7D if necessary
if mdl.projQ
    pptspred = proj_down(o,mdl.projtol,mdl.usv,'zero',mdl.zeroQ);
else
    pptspred = o;
end
ptspred = o;

%% Make predictions

switch mdl.method
    case 'pbary'
        mdl.data.pts = ptspred;
        mdl.data.ppts = pptspred;
        mdl.data.npts = size(ptspred,1);
        mdl.data.props = GB5DOF_setup(ptspred(:,1:4),ptspred(:,5:8),[0 0 1],'Ni',1);
        [intfacetIDs,databary,klist] = intersect_facet(mdl.mesh.ppts,mdl.mesh.sphK,pptspred,mdl.inttol,'nnMax',mdl.nnMax);
        ypred = get_interp(mdl.mesh,mdl.data,intfacetIDs,mdl.barytype,mdl.barytol);
    case 'gpr'
        if isfield(mdl,'gprMdl')
            gprMdl = mdl.gprMdl;
        else
            gprMdl = mdl.cgprMdl;
        end
        if nargout > 1
            [ypred,ysd,yint] = predict(gprMdl,pptspred);
        else
            ypred = predict(gprMdl,pptspred);
            [ysd,yint] = deal([]);
        end

    case 'gprmix'
        alpha = 0.05;
        [ypred,ysd,yint] = gprmix(mdl,[],[],pptspred,'gprMdl2',nv.gprMdl2);
        kfn = kmix(nv.kfntmp,nv.kfntmp2);
        a = sigfn(ypred{1});
        covmat = kfn(pptspred,pptspred,a);
        covmat = nearestSPD(covmat);
        if ~isempty(nv.nsamp)
            ypred = mvnrnd(ypred.',covmat,nv.nsamp);
        end

    case 'idw'
        ypred = idw(mdl.mesh.ppts,pptspred,mdl.mesh.props,mdl.r,mdl.L);

    case 'nn'
        ypred = mdl.mesh.props(dsearchn(mdl.mesh.ppts,pptspred));

    case {'egprm','egpr','gprm'}
        K = mdl.K;
        mixQ = mdl.mixQ;
        [~,ypred,ysd,~,yint] = ...
            egprm([],[],[],[],[],K,'mixQ',mixQ,...
            'o2',sqrt2norm(ptspred,'quat'),'egprmMdl',mdl,...
            'dispQ',false,'KdispQ',false,'egprmDispQ',false);
end

%% Apply no-boundary constraint

if nv.noboundaryQ

    % determine threshold
    if ischar(nv.wthreshold) && strcmpi(nv.wthreshold,'min') % use the minimum of the input disorientation angles and the constant 5 deg

        % compute input disorientation angles
        qd = disorientation(nv.qm,'cubic'); % there may be a faster way to do this since we just want the angle
        [w,~,~] = q2rot(qd); % disorientation angle

        % take the minimum
        wthreshold = min([min(w),deg2rad(5)]);

    elseif isscalar(nv.wthreshold)

        % use the supplied value
        wthreshold = nv.wthreshold;

    else

        % throw error
        error('Unrecognized wthreshold value.')

    end

    % apply constraint
    [ypred,ysd,yint] = applyNoBoundaryConstraint(ypred,ysd,yint,qm2,mdl,'weights',nv.weights,'wthreshold',wthreshold);

%     % compute disorientation angles
%     qd = disorientation(qm2,'cubic'); % there may be a faster way to do this since we just want the angle
%     [w,~,~] = q2rot(qd);
% 
%     % assign weights
%     eta = ones(size(w));
%     switch nv.weights
%         case 'rsw'
%             eta(w <= nv.wthreshold) = RSWfun(w(w <= nv.wthreshold),0,nv.wthreshold,0,1,1); % may want to use exponential correlation function instead?
%         case 'kernel'
%             dE = w/4; % need distance in dE rad = 0.25*w, since model uses dE rad (and hence kernel parameters are in dE rad)
%             switch lower(mdl.KernelInformation.Name)
%                 case 'exponential'
%                     eta = 1-exp(-dE/mdl.KernelInformation.KernelParameters(1));
%                 case 'squaredexponential'
%                     eta = 1-exp(-0.5*(dE/mdl.KernelInformation.KernelParameters(1)).^2);
%             end
%     end
% 
%     % update predictions
%     ypred = eta.*ypred;
%     if nargout > 1
%         ysd = eta.*ysd;
%         yint = eta.*yint;
%     end

end

end

% %--------------------------Helper Functions--------------------------------
% function fout = RSWfun(w,wmin,wmax,ymin,ymax,a)
% 
% fout = zeros(size(w));
% fout(w ~= wmin) = ((ymax-ymin)*sin((pi/2)*(w(w~=wmin)-wmin)/(wmax-wmin)).*(1-a*log(sin((pi/2)*(w(w~=wmin)-wmin)/(wmax-wmin))))+ymin);
% fout(w == wmin) = ymin; % when w == wmin you get 0*inf = nan, so we have to correct for this
% 
% end