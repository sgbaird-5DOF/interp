%INTERP5DOF_TEST  simple test case for five degree-of-freedom GB property interpolation
clear; close all

testnum = 2;

addpathdir({'cu2qu.m','q2rod.m','qmult.m','get_uuid.m','qmA2nA.m'})

%% Test Parameters
switch testnum
    case 1
        ninputpts = 388;
        npredpts = 1000;
        
    case 2
        ninputpts = 10000;
        npredpts = 10000;
end

%% random 5DOF parameters
[five,qm,nA] = get_five(ninputpts);
[five2,qm2,nA2] = get_five(npredpts);

%% get BRK function values
tic
y = GB5DOF_setup(five);
ytrue = GB5DOF_setup(five2);
toc

%% Interpolation
tstart = tic;
[ypred,interpfn,mdl,mdlpars] = interp5DOF(qm,nA,y,qm2,nA2,'gpr');
toc(tstart)

%% Results
errmetrics = get_errmetrics(ypred,ytrue,'dispQ',true);

%% Plotting
paperfigure();
parityplot(ytrue,ypred)


%-----------------------------CODE GRAVEYARD-------------------------------
%{
%         rmse = sqrt(immse(propList2,propOut));
%         mae = mean(abs(propList2-propOut));


%         o = get_ocubo(ndatapts);
%         o2 = get_ocubo(npredpts);
         
        %convert to 5DOF
%         five = GBoct2five(o);
%         five2 = GBoct2five(o2);

%         %unpack
%         qm = vertcat(five.q);
%         nA = vertcat(five.nA);
%         qm2 = vertcat(five2.q);
%         nA2 = vertcat(five2.nA);

%         [propOut,databary] = interp5DOF(qm,nA,propList,qm2,nA2,'pbary');

%         o = get_ocubo(ndatapts,'random',[],8);
%         o2 = get_ocubo(npredpts,'random',[],15);

% %% convert to octonions
% o = GBfive2oct(qm,nA);
% o2 = GBfive2oct(qm2,nA2);


%}