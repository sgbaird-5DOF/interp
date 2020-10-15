function [ypred,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,propList,qm2,nA2,method,NV)
arguments
    qm
    nA
    propList(:,1)
    qm2
    nA2
    method = 'gpr'
    NV.databary = []
    NV.facetIDs = []
    NV.dataprops = []
    NV.modelparsspec = struct()
    NV.brkQ(1,1) logical = true
    NV.gpropts = struct.empty
    NV.uuid = get_uuid()
    NV.o = []
    NV.o2 = []
end
% INTERP5DOF  Convert misorientation and boundary plane normal 5DOF input
% data to a closed, octonion, hyperspherical mesh and interpolate property
% values for arbitrary grain boundaries using spherical barycentric
% coordinates, planar barycentric coordinates, or a Gaussian process
% regression.
%--------------------------------------------------------------------------
% Inputs:
%  qm - list of misorientation quaternions (data), as in qmult(qB,qinv(qA))
%  for grains A and B, where qA and qB are in the sample frame (same for
%  qm2)
%
%  nA - list of boundary plane Cartesian unit normals (grain A frame)
%
%  propList - property value for each grain boundary (GB)
%
%  qm2 - list of misorientation quaternions for query GBs
%
%  nA2 - list of boundary plane Cartesian unit normals for query GBs (grain
%   A frame)
%
%  method - interpolation scheme to use. Possible methods are:
%    'sphbary' - spherical barycentric coordinates
%    'pbary' - planar barycentric coordinates
%    'gpr' - Gaussian process regression
%    'nn' - nearest neighbor interpolation
%    'insertnamehere' - no functionality, but contains instructions for
%      implementing your own interpolation scheme
%
%  NV - method-specific name-value pairs
%     method == 'sphbary' or 'pbary'
%       'databary' - supply barycentric coordinates to reduce redundant
%       computation if repeating interp5DOF calls for same list of qm/nA
%       pairs, even if property values are different. I.e. the interpolation
%       is calculated via a simple table lookup and a dot product.
%       facetprops must also be supplied.
%
%       'facetprops' - supply properties of each facet vertex (i.e. same
%       size as databary) for barycentric data interpolation. databary must
%       also be supplied
%
%       'brkQ' - logical, whether or not to calculate BRK energy values for
%       the query points to use for error calculations. If false and
%       dataprops is not supplied, then data.props is assigned a NaN vector
%       of the same size as mesh.props
%
%       'dataprops' - user supplied properties for query points for error
%       calculations. If not supplied, data.props depends on brkQ
%
%       'modelparsspec' - user supplied struct of model-specific
%       parameters, only gets used if databary and facetprops also supplied
%       via name-value pairs
%
%     method == 'gpr'
%       'gpropts' - options structure that will be passed to fitrgp() in place
%       of the defaults that are supplied later in this function.
%
% Outputs:
%  propOut - interpolated property values of queried grain boundaries
%
%  varargout
%   method == 'gpr'
%    varargout{1} - gprMdl, a Guassian Process Regression Object
%    varargout{2} - ysd, standard deviations of predicted values (propOut)
%    varargout{3} - yint, 95% confidence intervals of predicted values
%
%   method == 'sphbary' or 'pbary'
%    varargout{1} - databary, barycentric coordinates
%    varargout{2} - fname, which contains get_interp() workspace
%
% Usage:
%  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
%  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'sphbary')
%  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'pbary')
%
%  [propOut,gprMdl] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
%  [propOut,gprMdl,ysd,yint] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
%
%  [propOut] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr','gpropts',gpropts);
%
%  [propOut,databary,fname] = interp5DOF(qm,nA,propList,qm2,nA2,'sphbary')
%  [propOut,databary,fname] = interp5DOF(qm,nA,propList,qm2,nA2,'pbary')
%
% Simple Example Data
%  npts = 100;
%  qm = get_cubo(npts); nA = normr(rand(npts,3)); %random (qm,nA) pairs
%  propList = rand(npts); %random property values
%  qm2 = get_cubo(npts); nA2 = normr(rand(npts,3)); %random (qm,nA) pairs
%  %(I suggest you look at interp5DOF_test.m instead)
%
% Dependencies:
%  MATLAB 2019b or higher (mainly for the "arguments" syntax checking at
%  the beginning of functions, which is used extensively throughout)
%
%  Toolboxes
%  -Statistics and Machine Learning Toolbox
%  -Symbolic Math Toolbox (optional, for numStabBary.m)
%
%  Functions
%   see "FUNCTION DEPENDENCIES" section at end of this file (prior to "CODE
%   GRAVEYARD" section)
%
% Notes:
%  Simpler, plug & play, input/output version of run.m (fewer options)
%
%  Dependencies updated 2020-09-03 SGB
%
%  *If you have the Computer Vision toolbox, normr.m will be shadowed by
%  the corresponding function in the toolbox. However, both should produce
%  the same results for the purposes here.
%
%  **If you have addpathdir.m available on your path, all the other
%  dependencies should be added as long as the functions are in a
%  sub-folder of your current working directory. Alternatively, call the
%  following line of code while in the "code" directory, and then move
%  to your directory of choice:
%  addpathdir({'normr.m','GB5DOF_setup.m','cu2qu.m','q2rod.m','GBfive2oct.m','correctdis.m'})
%
%  In the context of this function, mesh is equivalent to predictors &
%  predictor responses, and data is equivalent to the query points you want
%  to interpolate or predict at.
%
%  Test functions (e.g. GBfive2oct_test.m, get_octpairs_test.m) are
%  available for many of the functions listed in "FUNCTION DEPENDENCIES" at
%  the end which can help with understanding and visualizing what each
%  function does. These are also very useful for debugging. The test
%  functions may have other dependencies than the ones listed in this file,
%  such as plot5DOF.m. There is also a "test" function for this function
%  (interp5DOF.m), namely interp5DOF_test.m
%
%  If you want to implement a custom interpolation scheme other than the
%  three available here, see the instructions beneath case 'insertnamehere'
%  in the switch method statement.
%
%  You can minimize this preamble text in MATLAB by clicking the "minus"
%  symbol at the top-left.
%
% Author(s): Sterling Baird
%
% Date: 2020-09-03
%--------------------------------------------------------------------------

%unpack (some) name-value pairs
brkQ = NV.brkQ;
uuid = NV.uuid;

%display method
disp(['method = ' method])

% add relevant folders to path (by searching subfolders for functions)
addpathdir({'normr.m','GB5DOF_setup.m','cu2qu.m','q2rod.m','GBfive2oct.m','correctdis.m','interp_gpr.m'})

%% convert to octonions & symmetrize
%predictor points
if isempty(qm) && isempty(nA) && ~isempty(NV.o)
    predinput = 'octonion';
    otmp = NV.o;
else
    predinput = '5dof';
    otmp = GBfive2oct(qm,nA);
end
wtol = 1e-6;
o = get_octpairs(otmp,'wtol',wtol);
nmeshpts = size(o,1);

%query points
if isempty(qm2) && isempty(nA2) && ~isempty(NV.o2)
    queryinput = 'octonion';
    otmp2 = NV.o2;
else
    queryinput = '5dof';
    otmp2 = GBfive2oct(qm2,nA2);
end
o2 = get_octpairs(otmp2,'wtol',wtol);
ndatapts = size(o2,1);

disp(['nmeshpts = ' int2str(nmeshpts) ', ndatapts = ' int2str(ndatapts)])

%% projection
%important to project both sets together so they use same SVD decomposition

projtol = 1e-4;
switch method
    case {'sphgpr','gpr'}
        projQ = false;
    otherwise
        projQ = true;
end
zeroQ = false;
o = normr(o);
o2 = normr(o2);
[a,usv] = proj_down([o;o2],projtol,'zeroQ',zeroQ);

%projected points
if size(a,2) == 7
    ppts = proj_down(o,projtol,usv,'zeroQ',zeroQ);
	ppts2 = proj_down(o2,projtol,usv,'zeroQ',zeroQ);
else
    error("Input doesn't have degenerate dimension or has too few (i.e. check input data), or try reducing proj_down tolerance input (tol)")
end

%% mesh and data struct setup
%mesh
mesh.pts = o;
mesh.ppts = ppts;
mesh.props = propList;
mesh.npts = nmeshpts;

%data
data.pts = o2;
data.ppts = ppts2;
data.npts = ndatapts;

if isempty(NV.dataprops)
    if NV.brkQ
        for i = 1:data.npts
            om1 = qu2om(o2(i,1:4));
            om2 = qu2om(o2(i,5:8));
            data.props(i) = GB5DOF(om1,om2,'Ni');
        end
    else
        data.props = nan(size(ppts2,1),1);
    end
else
    data.props = NV.dataprops;
end
ytrue = data.props;

%% additional variables
% current date and time
starttime = datetime(clock);
% number of cores (i.e. parfor workers)
p = gcp;
ncores = p.NumWorkers;
%git commit version
gitcommit = get_gitcommit();

%% package into struct
%general model variables 
mdlgen = var_names(brkQ,method,projtol,zeroQ,usv,starttime,ncores,...
    gitcommit,uuid,predinput,queryinput,projQ);
%general parameters
mdlparsgen = var_names(brkQ,method,projtol,zeroQ,starttime,nmeshpts,...
    ndatapts,ncores,gitcommit,uuid,predinput,queryinput,projQ);

%% helper functions
%function to concatenate structures with all different fields (no common)
structcat = @(S1,S2) table2struct([struct2table(S1,'AsArray',true),struct2table(S2,'AsArray',true)]);

%% method-specific interpolation
tic
switch method
    case {'sphbary','pbary'}
        if isempty(NV.databary)
            %% get triangulation
            [pptstmp,usvtri] = proj_down(o,projtol,'zeroQ',true);
            K = sphconvhulln(pptstmp);
            mesh.pts = proj_up(pptstmp,usvtri);
            
            %% compute intersecting facet IDs
            nnMax = 10;
            inttol = 0.01;
            disp('intersect_facet')
            intfacetIDs = intersect_facet(ppts,K,ppts2,inttol,'inttype','planar','nnMax',nnMax);
            
            %% mesh triangulation and filename
            mesh.sphK = K;
            
            %% interpolation
            disp('interpolation')
            %method-specific parameters
            switch method
                case 'sphbary'
                    barytol = 0.2;
                    barytype = 'spherical';
                    mesh.ppts = normr(mesh.ppts);
                    data.ppts = normr(data.ppts);
                    
                case 'pbary'
                    barytol = 0.1;
                    barytype = 'planar';
            end
            %interpolation
            [ypred,databary,facetprops,facetIDs,barypars] = get_interp(mesh,data,intfacetIDs,barytype,barytol);
            
            %model command and interpolation function
            mdlcmd = @(mesh,data,intfacetIDs,barytype,barytol) get_interp(mesh,data,intfacetIDs,barytype,barytol);
            interpfn = @(qm2,nA2) interp_bary(mesh,[],qm2,nA2,usv,zeroQ,barytype,barytol,projtol,nnMax,brkQ);
            
            %unpack intersection metrics
            nints = barypars.nints;
            numnonints = barypars.numnonints;
            int_fraction = barypars.int_fraction;
            
            %unpack NN extrapolation RMSE and MAE values
            nn_rmse = barypars.nn_errmetrics.rmse;
            nn_mae = barypars.nn_errmetrics.mae;        
            
            %model-specific variables
            mdlspec = var_names(databary,facetprops,barytol,barytype,inttol,...
                intfacetIDs,nnMax,facetIDs,barypars,nn_rmse,nn_mae);
                        
            %model-specific parameters
            mdlparsspec = var_names(barytol,barytype,inttol,nnMax,...
                nn_rmse,nn_mae,nints,numnonints,int_fraction,barypars);
            
        else
            %unpack
            databary = NV.databary;
            facetIDs = NV.facetIDs;
            
            %interpolate using supplied barycentric coordinates
            [ypred,facetprops,NNextrapID,nnList] = interp_bary_fast(ppts,ppts2,meshprops,databary,facetIDs);
            
            %model command and interp function
            mdlcmd = @(databary,facetprops) dot(databary,facetprops,2);
            interpfn = @(propList) interp_bary_fast(ppts,ppts2,meshprops,databary,facetIDs);
            
            %model-specific variables
            mdlspec = var_names(databary,facetprops,NNextrapID,nnList,facetIDs);
            
            %model-specific parameters
            mdlparsspec = NV.modelparsspec;
            
        end
        
    case {'sphgpr','gpr'}
        if projQ
            X = ppts;
            X2 = ppts2;
        else
            X = o;
            X2 = o2;
        end
        
        %gpr options
        if isempty(NV.gpropts)
            %% interp5DOF's default gpr options
            if nmeshpts <= 50000
                PredictMethod = 'exact';
                gpropts = {};
            else
                PredictMethod = 'bcd';
                gpropts = {'BlockSize',10000};
            end
            gpropts = [gpropts {'PredictMethod',PredictMethod}];
%             gpropts = {'PredictMethod',PredictMethod};
           
            if strcmp(method,'sphgpr')
                %squared exponential kernel function with octonion distance
                kfcn = @(XN,XM,theta) (exp(theta(2))^2)*exp(-(pdist2(XN,XM,@get_omega).^2)/(2*exp(theta(1))^2));
                theta0 = [deg2rad(10), std(propList)/sqrt(2)]; %initial length scale and noise
                gpropts = [gpropts,{'KernelFunction',kfcn,'KernelParameters',theta0}];
            end
            
        else
            % user-supplied gpr options
            gpropts = NV.gpropts;
            gproptnames = gpropts{1:2:end};
            gproptvals = gpropts{2:2:end};
            gproptstruct = cell2struct(gproptvals,gproptnames,2);
            
            %extract parameters (for table)
            G = gproptstruct;
            if isempty(fieldnames(G))
                error('user-specified gpropts is empty')
            end
            gproptshort = struct();
            if isfield(G,'HyperparameterOptimizationOptions')
                G1 = G.HyperparameterOptimizationOptions;
                if isfield(G1,'UseParallel')
                    gprParallelQ = G1.UseParallel;
                    gproptshort.gprParallelQ = gprParallelQ;
                end
                if isfield(G1,'Optimizer')
                    hyperoptimizer = G1.Optimizer;
                    gproptshort.hyperoptimizer = hyperoptimizer;
                end
                if isfield(G1,'MaxObjectiveEvaluations')
                    maxhyperobj = G1.MaxObjectiveEvaluations;
                    gproptshort.maxhyperobj = maxhyperobj;
                end
            end
            if isfield(G,'PredictMethod')
                PredictMethod = G.PredictMethod;
                gproptshort.PredictMethod = PredictMethod;
            end
            if isfield(G,'ActiveSetMethod')
                ActiveSetMethod = G.ActiveSetMethod;
                gproptshort.ActiveSetMethod = ActiveSetMethod;
            end
            if isfield(G,'FitMethod')
                FitMethod = G.FitMethod;
                gproptshort.FitMethod = FitMethod;
            end
        end
        
        %Gaussian process regression        
        if ~isempty(gpropts)
            gprMdl = fitrgp(X,propList,gpropts{:}) %#ok<NOPRT>
        else
            gprMdl = fitrgp(X,propList);
        end
        %compact the model
        cgprMdl = compact(gprMdl);
        
        %predictions ("interpolated" points)
        if ~strcmp(PredictMethod,'bcd')
            [ypred,ysd,yint] = predict(cgprMdl,X2);
        else
            ypred = predict(cgprMdl,X2);
        end
        
        mdlcmd = @(gprMdl,X2) predict(cgprMdl,X2);
        interpfn = @(qm2,nA2) interp_gpr(cgprMdl,qm2,nA2,projtol,usv);
        
        %model-specific variables
        if ~strcmp(PredictMethod,'bcd')
            mdlspec = var_names(cgprMdl,gpropts,ysd,yint);
        else
            mdlspec = var_names(cgprMdl,gpropts);
        end
        %model-specific parameters
        if exist('gproptshort','var') == 1
            if ~isempty(fieldnames(gproptshort))
                mdlparsspec = var_names(maxhyperobj,gproptshort);
            else
                evalc([(gproptstruct{end}) ' = G.(gproptstruct{end}']);
                warning(['no tracked options were contained in user-specified gpropts. Consider tracking all. Tracking added via evalc() for ' ...
                    gproptstruct{end-1}])
                if isstruct(gproptshort.(gproptstruct{end}))
                    error(['Tracking automatically added for ' gproptstruct{end-1} ' but the added tracking option cannot be a struct'])
                end
                mdlparsspec = struct(gproptstruct{end-1},gproptshort.gproptstruct{end});
            end
        else
            mdlparsspec = var_names(PredictMethod);
        end
        
    case 'nn'
        %nearest neighbors (NN)
        nnList = dsearchn(ppts,ppts2);
        
        mdlcmd = @(ppts,ppts2,propList) propList(dsearchn(ppts,ppts2));
        interpfn = @(qm2,nA2) interp_nn(ppts,qm2,nA2,projtol,usv,propList);
        
        %assign NN property values
        ypred = propList(nnList);
        
        %model-specific variables
        mdlspec.nnList = nnList;
        %model-specific parameters
        mdlparsspec = struct();
        
    case 'avg'
        % "interpolation" (just constant model)
        [ypred,yavg] = interp_avg(propList,ndatapts);
        
        mdlcmd = @(propList,ndatapts) interp_avg(propList,ndatapts);
        interpfn = @(qm2,nA2) repelem(yavg,ndatapts,1); %any new point gets assigned yavg
        
        %model-specific variables
        mdlspec.yavg = yavg;
        %model-specific parameters
        mdlparsspec = struct();
        
    case 'insertnamehere'
        %implement your own interpolation scheme here
        % specify a custom name next to case, e.g. case 'insertnamehere' --> case 'myinterpscheme'
        %
        % variables to use:
        %   ppts - predictors
        %   propList - reponses
        %
        % variables to define:
        %   propOut - interpolated values
        %   mdl - struct containing model variables
        %   mdlpars - struct containing model parameters
        %
        % update the documentation
        %
        % consider suggesting as an addition in the GitHub repository via
        % GitHub issue or preferably pull request
end
runtime = toc; %time elapsed to do the interpolation (method-specific portion)

%% append extra general variables
%parity variables
parity.ypred = ypred;
parity.ytrue = ytrue;

%model
mdlgen.ypred = ypred;
mdlgen.mdlcmd = mdlcmd;
mdlgen.interpfn = interpfn;
mdlgen.runtime = runtime;
mdlgen.mesh = mesh;
mdlgen.data = data;
mdlgen.parity = parity;
%parameters
mdlparsgen.runtime = runtime;
mdlparsgen.parity = parity;

%% concatenate structs
%model variables 
mdl = structcat(mdlgen,mdlspec);
%model parameters
mdlpars = structcat(mdlparsgen,mdlparsspec);

end

%%
%---------------------------FUNCTION DEPENDENCIES--------------------------
%  addpathdir.m **
%
%  GB5DOF_setup.m
%   -GB5DOF.m
%   -constructGBMatrices.m
%    --gmat2q.m
%    --q2gmat.m
%    --qinv_johnson.m (renamed to prevent conflict with CMU qinv.m)
%     ---qconj.m
%     ---qnorm.m
%    --qmultiply.m
%
%  GBfive2oct.m
%   -ax2qu.m
%   -qinv.m
%
%  get_interp.m
%   -get_omega.m
%   -projray2hypersphere.m
%    --numStabBary.m (optional)
%   -sphbary.m
%    --projfacet2hyperplane.m
%     ---projray2hyperplane.m
%   -sqrt2norm.m
%
%  get_octpairs.m
%   -misFZfeatures.mat
%   -GBdist4.m
%    --osymsets.m
%     ---get_sympairs.m
%      ----allcomb.m (FEX)
%      ----PGsymops.mat
%      ----PGnames.mat (optional)
%     ---osymset.m
%     ---proj_up.m
%	 --get_omega.m
%    --sqrt2norm.m
%    --zeta_min2.m (distinct from 'zeta_min' to prevent conflicts with
%    GBdist.m)
%   -q2rod.m
%   -disorientation.m
%   -addpathdir.m
%   -GBfive2oct.m (see above)
%   -GBoct2five.m
%    --disorientation.m (optional)
%    --findgeometry.m
%     ---inmisFZ.m
%      ----misFZcon.m (optional because misFZcon.m is called in findgeometry.m)
%       -----vert2con.m (FEX)
%     ---misFZcon.m
%      ----vert2con.m (FEX)
%     ---q2rod.m
%    --q2rod.m
%    --normr.m *
%    --qinv.m
%    --qu2ax.m
%   -symaxis.m
%   -mustBeSqrt2Norm.m (argument validation function)
%   -plotFZrodriguez_vtx.m (optional)
%    --plotFZrodriguez.m
%
%  intersect_facet.m
%   -projray2hypersphere.m
%    --numStabBary.m (optional)
%   -sphbary_setup.m
%    --sphbary.m
%     ---projfacet2hyperplane.m
%      ----projray2hyperplane.m
%
%  normr.m *
%
%  proj_down.m
%
%  sphconvhulln.m
%   -normr.m
%   -projfacet2hyperplane.m
%    --projray2hyperplane.m
%
%  var_names.m (convenient for packaging structures)
%%
%------------------------------CODE GRAVEYARD------------------------------
%{
            % project points to hyperplane (introduce deg dimension)
            pptstmp = projfacet2hyperplane(mean(ppts),ppts); %valid for max arc length < pi
            % rotate to remove newly introduced deg dimension
            pptstmp = proj_down(pptstmp,1e-4);
            K = convhulln(pptstmp);

        %         opts = getopts(6);


        %regression
%         hyperopts = struct('UseParallel',true);
%         opts = { ...
%             'OptimizeHyperparameters',{'KernelScale','Sigma'},...
%             'HyperparameterOptimizationOptions',hyperopts,...
%             'PredictMethod','bcd',...
%             'ActiveSetMethod','entropy',...
%             'FitMethod','fic'};

            if strcmp(method,'sphbary')
%                 ppts = normr(ppts);
%                 ppts2 = normr(ppts2);
            end


            %             mdlvars = [mdlvars {databary,facetprops}];
            %             propOut(nanID) = propList(nnList);



            %             mdl.databary
                        %             mdlvars = [mdlvars {databary,facetprops,fname,intfacetIDs}];
            %             varargout{1} = databary;
            %             varargout{2} = facetprops;
            %             varargout{3} = fname;


        %         varargout{1} = @(ppts2) dsearchn(ppts,ppts2);

        %         mdlvars = [mdlvars {gprMdl,ysd,yint}];
        %         varargout = {gprMdl,ysd,yint};

% mdl = table2struct(tblvertcat(struct2table(mdlgen),struct2table(mdlspec)));
% mdlpars = table2struct(tblvertcat(struct2table(mdlparsgen),struct2table(mdlparsspec)));

%         mdlcmd = @(ppts,propList,gpropts) predict(fitrgp(ppts,propList,gpropts{:}),ppts2);

% 
%             gpropts = struct(...
%                 'OptimizeHyperparameters',{'KernelScale','Sigma'},...
%                 'HyperparameterOptimizationOptions',hyperopts,...
%                 'PredictMethod','exact',...
%                 'ActiveSetMethod','entropy',...
%                 'FitMethod','fic');
%             %namedargs2cell(gpropts)


            structcat(rmfield(data,{'pts','ppts'}),...
                struct('pts',get_pts(qm2,nA2),'ppts',get_ppts(qm2,nA2))),...

            %maxhyperobj = ncores*2;
            %gprParallelQ = true;
            %hyperoptimizer = 'bayesopt';
            %%if ndatapts > 10000
            %%    PredictMethod = 'bcd'; %'exact', 'bcd'
            %%else
            %PredictMethod = 'exact';
            %%end
            %ActiveSetMethod = 'entropy';
            %FitMethod = 'fic';
            %hyperopts = struct('UseParallel',gprParallelQ,'Optimizer',hyperoptimizer,'MaxObjectiveEvaluations',maxhyperobj);
            %gpropts = { ...
            %'OptimizeHyperparameters',{'KernelScale','Sigma'},...
            %'HyperparameterOptimizationOptions',hyperopts,...
            %'PredictMethod',PredictMethod,...
            %'ActiveSetMethod',ActiveSetMethod,...
            %'FitMethod',FitMethod};

%             mesh.fname = 'meshtemp.mat';
            
            %assigning query point properties
%             data.fname = 'datatemp.mat';



%             mdlspec.intfacetcell = intfacetcell;


%             facetprops = NV.facetprops;
%             
%             %find NaN values & replace with NN values (NN extrapolation)
%             [NNextrapID,~] = isnan(databary);
%             nnList = dsearchn(ppts2(NNextrapID),ppts);
%             d = size(databary,2);
%             
%             % e.g. databary == [NaN NaN NaN], facetprops == [NaN NaN NaN]
%             % --> [1 0 0], [1.213 0 0], such that dot([1 0 0],[1.213 0 0])
%             % == 1.213, where 1.213 is the extrapolated NN value
%             databary(NNextrapID,1) = 1;
%             databary(NNextrapID,2:d) = 0;
%             facetprops(NNextrapID,1) = propList(nnList);
%             facetprops(NNextrapID,2:d) = 0;
            
%             switch NV.modelparsspec.method
%                 case 'sphbary'
%                     getinterpmethod = 'spherical';
%                 case 'pbary'
%                     getinterpmethod = 'planar';
%             end


%             propOut = dot(databary,facetprops,2);
            %             interpfn = @(propList) dot(databary,facetprops,2);

%             intfacetcell = {intfacetIDs};



            % additional input checking
%             if sum([isempty(NV.databary) isempty(NV.facetprops)]) == 1
%                 error('both databary and facetprops should be supplied simultaneously')
%             end
%             
%             if ~all(size(NV.databary) == size(NV.facetprops))
%                 error('databary and facetprops must have the same dimensions')
%             end
%             
%             if sum([isempty(NV.databary) isempty(NV.facetprops)])==1
%                 error('both databary and facetprops should be defined or neither')
%             end


%             gpropts = { ...
%             'OptimizeHyperparameters',{'KernelScale','Sigma'},...
%             'HyperparameterOptimizationOptions',hyperopts};



%             maxhyperobj = 30; %default
%             gprParallelQ = true;
%             [hyperoptimizer,PredictMethod,ActiveSetMethod,FitMethod]=deal('default');
%             hyperopts = struct('UseParallel',gprParallelQ);

            mdlparsspec = var_names(maxhyperobj,gprParallelQ,hyperoptimizer,PredictMethod,ActiveSetMethod,FitMethod);


%             [ActiveSetMethod,FitMethod]=deal('default');

% %projected points
% ppts = a(1:nmeshpts,:);
% ppts2 = a(nmeshpts+1:end,:);

%}
