function [propOut,interpfn,mdl,mdlpars] = interp5DOF_setup(ndatapts,npredpts,method,uuid,inputtype)
arguments
   ndatapts
   npredpts
   method = 'gpr'
   uuid = get_uuid()
   inputtype char {mustBeMember(inputtype,{'5dof','octonion'})} = '5dof'
end
%INTERP5DOF_SETUP  setup for interpolating five-degree-of-freedom property
%data using random octochorically sampled octonions

%% setup
addpathdir({'cu2qu.m','q2rod.m','qmult.m','get_ocubo.m'})

%random octonions
ocubotype = 'random';
ocuboseed1 = 8;
ocuboseed2 = 15;
%ocuboseed1 = 'shuffle';
%ocuboseed2 = 'shuffle';

five = get_five(ndatapts);
five2 = get_five(npredpts);

o = GBfive2oct(five);

% o = get_ocubo(ndatapts,ocubotype,[],ocuboseed1);
% o2 = get_ocubo(npredpts,ocubotype,[],ocuboseed2);

o = get_octpairs(o);
o2 = get_octpairs(o2);

genseed = 10;
rng(genseed);

%convert to 5DOF
% oct2fivemethod = 'reverse'; %'reverse', 'simple'
% five = GBoct2five(o,true,oct2fivemethod);
% five2 = GBoct2five(o2,true,oct2fivemethod);

%get BRK function values
propList = GB5DOF_setup(five);
dataprops = GB5DOF_setup(five2);

switch inputtype
    case '5dof'
        %unpack
        qm = vertcat(five.q);
        nA = vertcat(five.nA);
        qm2 = vertcat(five2.q);
        nA2 = vertcat(five2.nA);
    case 'octonion'
        [qm,nA,qm2,nA2]=deal([]);
end

%% interpolation
[propOut,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,propList,qm2,nA2,method,...
    'uuid',uuid,'o',o,'o2',o2,'dataprops',dataprops);

%% error values
proptrue = mdl.data.props;

e = proptrue-propOut; %error
ae = abs(e); %absolute error
mae = mean(ae); %mean absolute error
se = e.^2; %square error
mse = mean(se); %mean square error
rmse = sqrt(mse); %root mean square error
disp(['RMSE = ' num2str(rmse) ' J/m^2'])
disp(['MAE = ' num2str(mae) ' J/m^2'])
errmetrics = var_names(e,ae,mae,se,mse,rmse);

%% repackage model and parameters
mdlparsextra = var_names(ocubotype,ocuboseed1,ocuboseed2,genseed,errmetrics,rmse,mae,oct2fivemethod);
mdlextra = var_names(ocubotype,ocuboseed1,ocuboseed2,genseed,errmetrics,rmse,mae,oct2fivemethod);

%function to concatenate structures with all different fields
structcat = @(S1,S2) table2struct([struct2table(S1,'AsArray',true),struct2table(S2,'AsArray',true)]);

%concatenation
mdlpars = structcat(mdlpars,mdlparsextra);
mdl = structcat(mdl,mdlextra);
end


%---------------------------------CODE GRAVEYARD---------------------------
%{
switch method
    case 'gpr'
        [propOut,mdl,ysd,yint] = interp5DOF(qm,nA,propList,qm2,nA2,method);
        mdl.ysd = ysd;
        mdl.yint = yint;
    case {'sphbary','pbary'}
        [propOut,mdl.databary,mdl.facetprops,mdl.fname] = interp5DOF(qm,nA,propList,qm2,nA2,method);
    case 'nn'
        [propOut,mdl] = interp5DOF(qm,nA,propList,qm2,nA2,method);
end

% mse = immse(propList2,propOut);


% ocuboOpts = var_names(ocubotype,ocuboseed1,ocuboseed2,genseed);


%o = get_ocubo(ndatapts,'random');


%o2 = get_ocubo(npredpts,'random');
%}
