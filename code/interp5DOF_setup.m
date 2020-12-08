function [ypred,interpfn,mdl,mdlpars] = interp5DOF_setup(ndatapts,npredpts,method,datatype,epsijk,NV)
arguments
    ndatapts
    npredpts
    method char = 'gpr'
    datatype char {mustBeMember(datatype,{'brk','kim','rohrer-Ni','rohrer-test','rohrer-brk-test'})} = 'brk'
    epsijk(1,1) double = 1
    NV.pgnum(1,1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
    NV.uuid = get_uuid()
    NV.K(1,1) double = 1 %# of VFZO ensembles
    NV.sigma(1,1) double = 0 %mJ/m^2, standard deviation to add to y
end
%INTERP5DOF_SETUP  setup for interpolating five-degree-of-freedom property
%data using random octochorically sampled octonions

%% setup
addpathdir({'cu2qu.m','q2rod.m','qmult.m','get_ocubo.m','get_uuid.m'})

%unpack
pgnum = NV.pgnum;
uuid = NV.uuid;
K = NV.K;
sigma = NV.sigma;

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
        ytmp = S.propList;
        
        npts = ndatapts+npredpts;
        c = cvpartition(npts,'Holdout',npredpts);
        id1 = training(c);
        id2 = test(c);
        
        %split 5dof parameters
        five = struct('q',q(id1,:),'nA',nA(id1,:));
        five2 = struct('q',q(id2,:),'nA',nA(id2,:));
        
        %split properties
        y = ytmp(id1);
        ytrue = ytmp(id2);
        
    case 'brk'
        %random 5dof parameters
        five = get_five(ndatapts);
        five2 = get_five(npredpts);
        
        %get BRK function values
        y = GB5DOF_setup([],five.q,five.nA,epsijk);
        ytrue = GB5DOF_setup([],five2.q,five2.nA,epsijk);
        
    case {'rohrer-Ni','rohrer-test','rohrer-brk-test'}
        switch datatype
            case 'rohrer-Ni'
%                 datfpath = '../../TJ2GBE/TJdata/triples_Ni_0131_21520.txt';
                datfpath = '../../TJ2GBE/TJdata/triples_Ni_cat_88092.dat';
                resEfpath = '../../TJ2GBE/output/Ni_cat_88092_Cub_lamb300mresE.txt';
                
%                 load('../../TJ2GBE/output/Ni_0131_21520_Cub.mat','EAs','norms','resE')
%                 %convert
%                 [q,nA] = TJ2five(EAs,norms,epsijk);
%                 oct = TJ2oct(EAs,norms,epsijk);
            case {'rohrer-test','rohrer-brk-test'}
                datfpath = '../../TJ2GBE/TJdata/triples_30000.dat';
%                 [~,EAs,norms] = read_dat(datfpath);
                resEfpath = '../../TJ2GBE/TJdata/trueE_30000.txt';
        end
        [q,nA] = datfile2five(datfpath,0,epsijk);
        
        % assign GBE values
        switch datatype
            case {'rohrer-test','rohrer-Ni'}
                resE = importdata(resEfpath);
                resE(resE < 0) = 0;
            case {'rohrer-brk-test'}
                resE = GB5DOF_setup([],q,nA,epsijk);
        end
        
        % average repeats, remove repeats except 1
%         qnA = [q,nA];
%         [qnA,resE] = avgrepeats(qnA,resE);
%         q = qnA(:,1:4);
%         nA = qnA(:,5:8);

        %unpack
        assert(isvector(resE),'resE should be a vector');
        ytmp = resE(:);
        
        oct = five2oct(q,nA,epsijk);
        
        npts = ndatapts+npredpts;
        c = cvpartition(npts,'Holdout',npredpts);
        id1 = training(c);
        id2 = test(c);
        
        %split 5dof parameters
        five = struct('q',q(id1,:),'nA',nA(id1,:));
        five2 = struct('q',q(id2,:),'nA',nA(id2,:));
        
        o = oct(id1,:);
        o2 = oct(id2,:);
        
        %split properties
        y = ytmp(id1);
        ytrue = ytmp(id2);
        
    case 'olmsted-Ni'
        
        
end

noisetype = 'normal';
switch noisetype
    case 'normal'
        y = normrnd(y,sigma);
    case 'uniform'
        y = y + sigma*2*(rand(size(y))-0.5); %uniform
end

%unpack
qm = vertcat(five.q);
nA = vertcat(five.nA);
qm2 = vertcat(five2.q);
nA2 = vertcat(five2.nA);

%% interpolation
[ypredlist,interpfnlist,mdllist,mdlparslist] = deal(cell(K,1));
for k = 1:K
    switch datatype
        case 'rohrer-Ni'
            % use octonions obtained via GBlab2oct
            [qm,nA,qm2,nA2]=deal([]);
            [ypredlist{k},interpfnlist{k},mdllist{k},mdlparslist{k}] = interp5DOF(qm,nA,y,qm2,nA2,method,...
                'pgnum',pgnum,'uuid',uuid,'ytrue',ytrue,'o',o,'o2',o2);
        otherwise
            [ypredlist{k},interpfnlist{k},mdllist{k},mdlparslist{k}] = interp5DOF(qm,nA,y,qm2,nA2,method,...
                'pgnum',pgnum,'uuid',uuid,'ytrue',ytrue);
    end
end
ypredtmp = [ypredlist{:}];
ypred = median(ypredtmp,2);
interpfn = interpfnlist{1};

%% error values
proptrue = mdllist{1}.data.props;

errmetrics = get_errmetrics(ypred,proptrue);

rmse = errmetrics.rmse;
mae = errmetrics.mae;
disp(['input sigma = ',num2str(sigma),' J/m^2'])
disp(['RMSE = ' num2str(rmse) ' J/m^2'])
disp(['MAE = ' num2str(mae) ' J/m^2'])

%% repackage model and parameters
% variables and parameters to prepend
mdlpre = var_names(datatype,sigma);
mdlparspre = var_names(datatype,sigma);

% variables and parameters to append
mdlparsextra = var_names(genseed,errmetrics,rmse,mae);
mdlextra = var_names(genseed,errmetrics,rmse,mae,seedstruct);

%concatenation
mdlpars = structhorzcat(mdlparspre,mdlparslist{1},mdlparsextra);
mdl = structhorzcat(mdlpre,mdllist{1},mdlextra);
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

%         [q,nA] = TJ2five(EAs,norms,'francis');

%         [pA,pB,mA] = get_qmA(ndatapts);
%         [pA2,pB2,mA2] = get_qmA(npredpts);
%         y = GB5DOF_setup(pA,pB,mA,epsijk);
%         ytrue = GB5DOF_setup(pA,pB,mA,epsijk);

%                 oct2 = TJ2oct(EAs,norms,epsijk);

                [~,e1,e2,e3,m1,m2,m3] = datfile2em(fpath,nheaderlines);
                [q,nA] = em2five(e1,e2,e3,m1,m2,m3,epsijk);


%}
