function T = predict5DOF_timing(qm2,nA2,mdl,nv)
arguments
    qm2(:,4) double % query misorientations
    nA2(:,3) double % query BP normals
    mdl struct % trained 'gpr' model from interp5DOF
    nv.epsijk(1,1) double = 1
    nv.dispQ(1,1) logical = true
end
% PREDICT5DOF_TIMING  Time each stage of predict5DOF for a trained GPR model
%--------------------------------------------------------------------------
% Runs the same steps as predict5DOF (without the no-boundary constraint)
% and reports how long each one takes. Symmetrizing the query GBs into the
% VFZ costs the same per query regardless of the number of training points
% (N). The GPR mean costs O(N) per query. With the exact method, the
% standard deviation and confidence interval cost O(N^2) per query.
%
% Inputs:
%  qm2 - list of misorientation quaternions of the query points, as in
%        predict5DOF
%
%  nA2 - list of boundary plane Cartesian unit normals (grain A frame) of
%        the query points
%
%  mdl - trained model, as returned by interp5DOF(...,'gpr')
%
% Outputs:
%  T - table with the time spent in each stage, in total and per query
%
% Usage:
%  [~,~,mdl] = interp5DOF(qm,nA,y,qm2,nA2,'gpr');
%  T = predict5DOF_timing(qm2,nA2,mdl);
%
% Date: 2026-10-02
%--------------------------------------------------------------------------

npts = size(qm2,1);
if isfield(mdl,'gprMdl')
    gprMdl = mdl.gprMdl;
else
    gprMdl = mdl.cgprMdl;
end

Stage = {'five2oct';'get_octpairs (symmetrize into VFZ)';'normr + proj_down';...
    'predict: mean only';'predict: mean, sd and interval'};
Seconds = nan(numel(Stage),1);

tic
o = five2oct(qm2,nA2,nv.epsijk);
Seconds(1) = toc;

tic
o = get_octpairs(o,nv.epsijk,'oref',mdl.oref,'dispQ',false,'IncludeTies',false);
Seconds(2) = toc;

tic
o = normr(o);
if mdl.projQ
    X = proj_down(o,mdl.projtol,mdl.usv,'zero',mdl.zeroQ);
else
    X = o;
end
Seconds(3) = toc;

tic
ypred = predict(gprMdl,X); %#ok<NASGU>
Seconds(4) = toc;

% 'bcd' predictions don't provide sd or intervals
if ~strcmp(gprMdl.PredictMethod,'bcd')
    tic
    [ypred,ysd,yint] = predict(gprMdl,X); %#ok<ASGLU>
    Seconds(5) = toc;
end

MillisecondsPerQuery = 1e3*Seconds/npts;
T = table(Stage,Seconds,MillisecondsPerQuery);

if nv.dispQ
    fprintf('predict5DOF stages for %d query GBs (%s prediction, %d active set points):\n',...
        npts,gprMdl.PredictMethod,size(gprMdl.ActiveSetVectors,1))
    disp(T)
end

end
