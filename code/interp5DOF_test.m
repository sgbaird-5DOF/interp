%interp5DOF test
clear; close all

testnum = 1;

addpathdir({'cu2qu.m','q2rod.m','qmult.m','get_uuid.m','qmA2nA.m'})

switch testnum
    case 1
        %random octonions
        ndatapts = 388;
        npredpts = 10000;
%         o = get_ocubo(ndatapts,'random',[],8);
%         o2 = get_ocubo(npredpts,'random',[],15);
        five = get_five(ndatapts);
        five2 = get_five(npredpts);
        
        %unpack
        qm = vertcat(five.q);
        nA = vertcat(five.nA);
        qm2 = vertcat(five2.q);
        nA2 = vertcat(five2.nA);
        
        %convert to octonions
        o = GBfive2oct(qm,nA);
        o2 = GBfive2oct(qm2,nA2);
        
        %get BRK function values
        propList = GB5DOF_setup(five);
        propList2 = GB5DOF_setup(five2);
        
        [propOut,gprMdl] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr');
%         [propOut,databary] = interp5DOF(qm,nA,propList,qm2,nA2,'pbary');
        
        errmetrics = get_errmetrics(propOut,propList2);
        rmse = errmetrics.rmse;
        mae = errmetrics.mae;
        disp(['RMSE = ' num2str(rmse) ' J/m^2'])
        disp(['MAE = ' num2str(mae) ' J/m^2'])
end

figure
parityplot(propOut,propList2)


%-----------------------------CODE GRAVEYARD-------------------------------
%{
%         rmse = sqrt(immse(propList2,propOut));
%         mae = mean(abs(propList2-propOut));


%         o = get_ocubo(ndatapts);
%         o2 = get_ocubo(npredpts);
         
        %convert to 5DOF
%         five = GBoct2five(o);
%         five2 = GBoct2five(o2);
%}