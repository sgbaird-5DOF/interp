function [ypred,interpfn,mdl,mdlpars] = interp5DOF_setup(ndatapts,npredpts,method,datatype,NV)
arguments
   ndatapts
   npredpts
   method char = 'gpr'
   datatype char {mustBeMember(datatype,{'brk','kim','rohrer-Ni'})} = 'brk'
   NV.pgnum(1,1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
   NV.uuid = get_uuid()
end
%INTERP5DOF_SETUP  setup for interpolating five-degree-of-freedom property
%data using random octochorically sampled octonions

%% setup
addpathdir({'cu2qu.m','q2rod.m','qmult.m','get_ocubo.m','get_uuid.m'})

%unpack
pgnum = NV.pgnum;
uuid = NV.uuid;

%seed
seedstruct = rng;
genseed = seedstruct.Seed;

disp(['datatype = ' datatype])

switch datatype
    case 'kim'
        %load 5dof
        fname = 'Kim2011_Fe_oct_GBE.mat'; %produced via Kim2oct.m
        addpathdir({fname})
        S = load(fname,'five','propList');
        
        %unpack
        fivetmp = S.five;
        q = vertcat(fivetmp.q);
        nA = vertcat(fivetmp.nA);
        proptmp = S.propList;
        
        npts = ndatapts+npredpts;
        c = cvpartition(npts,'Holdout',npredpts);
        id1 = training(c);
        id2 = test(c);
                
        %split 5dof parameters
        five = struct('q',q(id1,:),'nA',nA(id1,:));
        five2 = struct('q',q(id2,:),'nA',nA(id2,:));

        %split properties
        propList = proptmp(id1);
        dataprops = proptmp(id2);
        
    case 'brk'
        %random 5dof parameters
        five = get_five(ndatapts);
        five2 = get_five(npredpts);
        
        %get BRK function values
        propList = GB5DOF_setup(five);
        dataprops = GB5DOF_setup(five2);
        
    case 'rohrer-Ni'
        load('../../TJ2GBE/output/Ni_0131_21520_Cub.mat','EAs','norms','resE')
        
        %convert
        [q,nA] = TJ2five(EAs,norms);
        %unpack
        assert(isvector(resE),'resE should be a vector');
        proptmp = resE(:);
        
        npts = ndatapts+npredpts;
        c = cvpartition(npts,'Holdout',npredpts);
        id1 = training(c);
        id2 = test(c);
        
        %split 5dof parameters
        five = struct('q',q(id1,:),'nA',nA(id1,:));
        five2 = struct('q',q(id2,:),'nA',nA(id2,:));
        
        %split properties
        propList = proptmp(id1);
        dataprops = proptmp(id2);
        
end

%unpack
qm = vertcat(five.q);
nA = vertcat(five.nA);
qm2 = vertcat(five2.q);
nA2 = vertcat(five2.nA);

%% interpolation
[ypred,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,propList,qm2,nA2,method,...
    'pgnum',pgnum,'uuid',uuid,'dataprops',dataprops);

%% error values
proptrue = mdl.data.props;

errmetrics = get_errmetrics(ypred,proptrue);

rmse = errmetrics.rmse;
mae = errmetrics.mae;
disp(['RMSE = ' num2str(rmse) ' J/m^2'])
disp(['MAE = ' num2str(mae) ' J/m^2'])

%% repackage model and parameters
% variables and parameters to prepend
mdlpre = var_names(datatype);
mdlparspre = var_names(datatype);

% variables and parameters to append
mdlparsextra = var_names(genseed,errmetrics,rmse,mae);
mdlextra = var_names(genseed,errmetrics,rmse,mae,seedstruct);

%concatenation
mdlpars = structhorzcat(mdlparspre,mdlpars,mdlparsextra);
mdl = structhorzcat(mdlpre,mdl,mdlextra);
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


e = proptrue-propOut; %error
ae = abs(e); %absolute error
mae = mean(ae); %mean absolute error
se = e.^2; %square error
mse = mean(se); %mean square error
rmse = sqrt(mse); %root mean square error
errmetrics = var_names(e,ae,mae,se,mse,rmse);



%convert to 5DOF
% oct2fivemethod = 'reverse'; %'reverse', 'simple'
% five = GBoct2five(o,true,oct2fivemethod);
% five2 = GBoct2five(o2,true,oct2fivemethod);


% ocubotype = 'random';
% ocuboseed1 = 8;
% ocuboseed2 = 15;
%ocuboseed1 = 'shuffle';
%ocuboseed2 = 'shuffle';


% o = get_ocubo(ndatapts,ocubotype,[],ocuboseed1);
% o2 = get_ocubo(npredpts,ocubotype,[],ocuboseed2);
% genseed = 10;
% rng(genseed);

mdlparsextra = var_names(ocubotype,ocuboseed1,ocuboseed2,genseed,errmetrics,rmse,mae,oct2fivemethod);
mdlextra = var_names(ocubotype,ocuboseed1,ocuboseed2,genseed,errmetrics,rmse,mae,oct2fivemethod);

% o = get_octpairs(o);
% o2 = get_octpairs(o2);

% rng('shuffle'); %to prevent getting '10' as the seed
% pause(5)



        %load 5dof
        fname = 'Kim2011_Fe_oct_GBE.mat'; %produced via Kim2oct.m
        addpathdir({fname})
        S = load(fname);
        
        %extract
        fivetmp = S.five;
        qtmp = vertcat(fivetmp.q);
        nAtmp = vertcat(fivetmp.nA);
        
        %get random IDs
        npts = size(qtmp,1);
        id1 = randi(npts,ndatapts,1);
        id2 = randi(npts,npredpts,1);
        
        %5dof parameters (randomly sampled from Kim Fe GBs)
        five.q = qtmp(id1,:);
        five.nA = nAtmp(id1,:);
        
        %five2
        five2.q = qtmp(id2,:);
        five2.nA = nAtmp(id2,:);  



   NV.inputtype char {mustBeMember(NV.inputtype,{'5dof','octonion'})} = '5dof'


inputtype = NV.inputtype;
'o',o,'o2',o2,


        switch inputtype
            case 'octonion'
                %octonions
                o = otmp(id1,:);
                o2 = otmp(id2,:);
                
            case '5dof'
                o = [];
                o2 = [];
        end


        switch inputtype
            case 'octonion'
                %convert to octonions
                o = GBfive2oct(five);
                o2 = GBfive2oct(five2);
                
            case '5dof'
                o = [];
                o2 = [];
        end

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


%function to concatenate structures with all different fields
structcat = @(S1,S2,S3) ...
    table2struct([...
    struct2table(S1,'AsArray',true),...
    struct2table(S2,'AsArray',true),...
    struct2table(S3,'AsArray',true)]); %distinct from structvertcat.m


        %get random IDs to split up the data
        npts = size(qtmp,1);
        idfull = 1:npts;
        id1 = randi(npts,ndatapts,1);
        idsub = setdiff(idfull,id1);
        id2 = randi(npts,npredpts,1);



%}
