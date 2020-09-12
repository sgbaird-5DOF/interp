%interp5DOF test
clear; close all

testnum = 1;

addpathdir({'cu2qu.m','q2rod.m','qmult.m'})

switch testnum
    case 1
        %random octonions
        ndatapts = 388;
        npredpts = 388;
        o = get_ocubo(ndatapts,'random',[],8);
        o2 = get_ocubo(npredpts,'random',[],15);
%         o = get_ocubo(ndatapts);
%         o2 = get_ocubo(npredpts);
         
        %convert to 5DOF
        five = GBoct2five(o);
        five2 = GBoct2five(o2);
        
        %unpack
        qm = vertcat(five.q);
        nA = vertcat(five.nA);
        qm2 = vertcat(five2.q);
        nA2 = vertcat(five2.nA);
        
        %get BRK function values
        propList = GB5DOF_setup(five);
        propList2 = GB5DOF_setup(five2);
        
        [propOut,gprMdl] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr');
%         [propOut,databary] = interp5DOF(qm,nA,propList,qm2,nA2,'pbary');
        
        rmse = sqrt(immse(propList2,propOut));
        mae = mean(abs(propList2-propOut));
        disp(['RMSE = ' num2str(rmse) ' J/m^2'])
        disp(['MAE = ' num2str(mae) ' J/m^2'])
end

figure
parityplot(propOut,propList2)
