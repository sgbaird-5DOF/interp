%verify that id sequence leads to expected energies for [100] STGBs

property_datapre = importdata('../olm_properties.txt'); pdata = property_datapre.data; %properties for 388 Ni GB's
dat100pre = importdata('stgb100_en.txt'); dat100 = dat100pre.data; %energy / geodesic angle plotted for subset of GB's


id100 = dat100(:,1);
ang100 = dat100(:,2);

% do the global properties for the dataset agree with the tabulation given
% by dat100? 

NiGBE_100_test1 = [0; pdata(id100(2:end-1),2); 0]; %NB's are added
% NiGBE_100_test2 = dat100(:,3); %already includes NB's

%% Plot GBE energy

figure
plot(ang100,NiGBE_100_test1)
xlabel('tilt angle')
ylabel('GB energy')
